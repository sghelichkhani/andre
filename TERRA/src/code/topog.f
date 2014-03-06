*DK geoid
      subroutine geoid(g,h,w,gh,plm,rr,lmax)
 
c...  This routine computes the geoid g from the temperature temp
c...  and boundary deformation h.
 
      include 'size.h'
      include 'pcom.h'
      real g((nt+1)**2,nd,2), h((nt+1)**2,nd,2)
      real w((nt+1)**2,nd,*), gh(2,(lmax+1)*(lmax+2)/2), tavg(nr+1)
      common /temp/ tpb(nt+1), temp((nt+1)**2,nd,nr+1), tpe(nt+1)
      common /radl/ rshl(nr+1), ird
      common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
      common /eos2/ rhocv(nr+1), gamma(nr+1), alpha(nr+1), prf(nr+1),
     &               cond(nr+1)
 
      call layrav(temp,tavg)
 
c...  Load array w with the nodal density fluctuation times
c...  4*pi*G*r*dr/g.
 
      a0 = 4.*3.14159265*6.6732e-11/grv(1)
 
      do ir=2,nr+1
 
         rt = 0.5*(rshl(ir-1) + rshl(ir))
         rb = rshl(nr+1)
         if(ir .le. nr) rb = 0.5*(rshl(ir+1) + rshl(ir))
 
         a1 = a0*(rt**3 - rb**3)/(3.*rshl(ir))*rhorf(ir)*alpha(ir)
         a1 = 0.
 
         do id=1,nd
            do ii=1,(nt+1)**2
               w(ii,id,ir) = a1*(tavg(ir) - temp(ii,id,ir))
            end do
         end do
 
      end do
 
c...  Include effect of topography on the top and bottom boundaries.
 
      a1 = a0*rhorf(1)*rshl(1)
      a2 = a0*(9910. - rhorf(nr+1))*rshl(nr+1)
 
      do id=1,nd
         do ii=1,(nt+1)**2
            w(ii,id,1)    = a1*h(ii,id,1)
            w(ii,id,nr+1) = a2*h(ii,id,2) + w(ii,id,nr+1)
         end do
      end do
 
c...  Compute the spherical harmonic geoid coefficients gh from the
c...  density variations within the shell and from boundary deflections.
 
      call geoidhar(gh,w,plm,rr,lmax)
 
c...  Convert from spherical harmonic representation of geoid to
c...  finite element nodal representation of geoid and gravity fields.
 
      call ndgeoid(g,gh,plm,lmax)
 
      if(mynum .eq. 0) then
         write(6,'(/"   Harmonic Coefficients for Geoid Field:"/)')
         ihar = 1
         do l=1,16
            do m=0,l
               ihar = ihar + 1
               write(6,'(2i5,1p2e15.7)') l, m, gh(1,ihar), gh(2,ihar)
            end do
         end do
      end if
 
      end
*dk geoidhar
      subroutine geoidhar(gh,d,plm,rr,lmax)
 
c...  This routine generates spherical harmonic geoid coefficients gh
c...  resulting from the nodal finite element density field d.
 
      include 'size.h'
      include 'pcom.h'
      real d((nt+1)**2,nd,nr+1), gh(2,(lmax+1)*(lmax+2)/2)
      real plm((lmax+1)*(lmax+2)), csm(0:128,2), rr(nr+1,lmax)
      common /radl/ rshl(nr+1), ird
      common /mesh/ xn((nt+1)**2,nd,3)
      common /ndar/  arn((nt+1)**2),  arne((nt+1)**2),
     &              rarn((nt+1)**2), rarne((nt+1)**2)
 
      r4pi = 0.125/asin(1.)
 
      do l=1,lmax
         do ir=1,nr+1
            rr(ir,l) = (rshl(ir)/rshl(1))**(l+1)
         end do
      end do
 
      call nulvec(gh, (lmax+1)*(lmax+2))
 
      do ii=1,(nt+1)**2
 
         a0 = r4pi*arn(ii)
 
         do id=1,nd
 
            phi = atan2(xn(ii,id,2) + 1.e-30, xn(ii,id,1))
 
            do m=0,lmax
               csm(m,1) = cos(m*phi)
               csm(m,2) = sin(m*phi)
            end do
 
            if(mod(id, 5) .eq. 1)
     &         call plmbar(plm,plm,lmax,xn(ii,id,3),0)
 
            k = 1
 
            do l=1,lmax
 
               aa = a0/(2*l + 1)
 
               do m=0,l
 
                  k = k + 1
 
                  do ir=1,nr+1
 
                     bb = aa*rr(ir,l)*plm(k)*d(ii,id,ir)
 
                     gh(1,k) = gh(1,k) + bb*csm(m,1)
                     gh(2,k) = gh(2,k) + bb*csm(m,2)
 
                  end do
 
               end do
 
            end do
 
         end do
 
      end do
 
      if(nproc .gt. 1) call psumlong(gh, plm, (lmax+1)*(lmax+2))
 
      end
*dk ndgeoid
      subroutine ndgeoid(g,gh,plm,lmax)
 
c...  This routine uses the spherical harmonics gh for the geoid
c...  to generate the nodal geoid and gravity fields.
 
      include 'size.h'
      include 'pcom.h'
      real g((nt+1)**2,nd,2), gh(2,(lmax+1)*(lmax+2)/2)
      real plm((lmax+1)*(lmax+2)), csm(0:128,2)
      common /radl/ rshl(nr+1), ird
      common /mesh/ xn((nt+1)**2,nd,3)
      common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
 
      call nulvec(g, (nt+1)**2*nd*2)
 
      a0 = grv(1)/rshl(1)
 
      do ii=1,(nt+1)**2
 
         do id=1,nd
 
            phi = atan2(xn(ii,id,2) + 1.e-30, xn(ii,id,1))
 
            do m=0,lmax
               csm(m,1) = cos(m*phi)
               csm(m,2) = sin(m*phi)
            end do
 
            if(mod(id, 5) .eq. 1)
     &         call plmbar(plm,plm,lmax,xn(ii,id,3),0)
 
            k = 3
 
            do l=2,lmax
 
               aa = a0*(l - 1)
 
               do m=0,l
 
                  k  = k + 1
                  bb = plm(k)*csm(m,1)*gh(1,k)
     &               + plm(k)*csm(m,2)*gh(2,k)
 
                  g(ii,id,1) = g(ii,id,1) + bb
                  g(ii,id,2) = g(ii,id,2) + bb*aa
 
               end do
 
            end do
 
         end do
 
      end do
 
      end


*dk sealevel
      subroutine sealevel(h,g,map2)
 
c...  This routine computes the surface topography relative to sea
c...  level, including contributions from dynamic topography, geoid,
c...  and continental crust.
 
      include 'size.h'
      include 'pcom.h'
      real h((nt+1)**2,nd,3), g((nt+1)**2,nd)

      common /ndar/  arn((nt+1)**2),  arne((nt+1)**2),
     &              rarn((nt+1)**2), rarne((nt+1)**2)
 
c...  Compute the surface topography relative to the geoid including
c...  the presence of continental crust.  Use 4500 m as the height
c...  associated with continental crust.

	integer map2(0:nt,nt+1,nd)
	integer i,j

      hmin = 1.e10
 
      do id=1,nd
         do ii=1,(nt+1)**2
				i=mod(ii-1,nt+1)
				j=ceiling(real(ii/(nt+1.0)))
            cont = 0.
            if(map2(i,j,id) .ne. 0) cont = 4500.
            h(ii,id,3) = h(ii,id,1) - g(ii,id) + cont
            hmin = min(hmin, h(ii,id,3))
         end do
      end do
 
      if(nproc .gt. 1) call pmin(hmin)
 
c...  Reference the surface topography such that zero occurs as its
c...  minimum value.
 
      hmax = 0.
 
      do id=1,nd
         do ii=1,(nt+1)**2
            h(ii,id,3) = h(ii,id,3) - hmin
            hmax = max(hmax, h(ii,id,3))
         end do
      end do
 
      if(nproc .gt. 1) call pmax(hmax)
 
c...  Use an iterative method for finding sea level relative to the
c...  surface topography field.
 
      radsq   = 6370.e3**2
      seavolm = 1.4e18
      sealevl = 0.7*hmax
 
      do it=1,5
 
         vol1 = 0.
         vol2 = 0.
 
         do id=1,nd
            do ii=1,(nt+1)**2
               aa = radsq*arn(ii)
               v1 = aa*(sealevl - h(ii,id,3))
               v2 = v1 + 10.*aa
               if(v1 .gt. 0.) vol1 = vol1 + v1
               if(v2 .gt. 0.) vol2 = vol2 + v2
            end do
         end do
 
         if(nproc .gt. 1) call psum(vol1, 1)
         if(nproc .gt. 1) call psum(vol2, 1)
 
         dv = vol2 - vol1
         ds = 0.
         if(dv .ne. 0.) ds = 10.*(seavolm - vol1)/dv
         sealevl = sealevl + ds
 
c...     write(*,10) ds, sealevl
c10      format(' ds = ',1pe10.3,'  sealevl = ',e10.3)
 
      end do
 
c...  Move reference value for surface topography to sea level.
 
      do id=1,nd
         do ii=1,(nt+1)**2
            h(ii,id,1) = h(ii,id,3) - sealevl
         end do
      end do
 
      end

*dk topogrph
      subroutine topogrph
 
c...  This routine computes the variations in topographic height h
c...  at the shell boundaries due to computed normal forces.  It then
c...  calls routines to compute the geoid, gravity, and sea level.
 
      include 'size.h'
      common /flds/ f((nt+1)**2*nd*(nr+1),8)
      common /topo/ h((nt+1)**2,nd,3), g((nt+1)**2,nd,2)
      common /velo/ upb(nt+1), u((nt+1)**2,nd,3,nr+1),  upe(nt+1)
      common /pres/ ppb(nt+1), pres((nt+1)**2,nd,nr+1), ppe(nt+1)
      common /mesh/ xn((nt+1)**2,nd,3)
      common /radl/ rshl(nr+1), ird
      common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
      common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
 
c...  Convert normal boundary velocity u and pressure pres to
c...  boundary topographic height assuming reference densities.
 
      k  = nr + 1
      b1 = 1./(rhorf(1)*grv(1))
      b2 = 1./((9910. - rhorf(k))*grv(k))
      a1 = 2.*b1*visc/(rshl(1)  - rshl(2))
      a2 = 2.*b2*visc/(rshl(nr) - rshl(k))
 
      do id=1,nd
         do ii=1,(nt+1)**2
            h(ii,id,1) = a1*((u(ii,id,1, 2)*xn(ii,id,1)
     &                      + u(ii,id,2, 2)*xn(ii,id,2))
     &                      + u(ii,id,3, 2)*xn(ii,id,3))
     &                      + b1*pres(ii,id,1)
            h(ii,id,2) = a2*((u(ii,id,1,nr)*xn(ii,id,1)
     &                      + u(ii,id,2,nr)*xn(ii,id,2))
     &                      + u(ii,id,3,nr)*xn(ii,id,3))
     &                      - b2*pres(ii,id,k)
         end do
      end do
 
c...  Compute geoid and surface gravity.
 
      call geoid(g,h,f,f(1,2),f(1,3),f(1,4),32)
 
c...  Compute surface topography and sea level.
 
c     call sealevel(h,g,map2)
 
c     call uharmonic(u,f,f(1,2),mt/2)
c     call uharmonic(u(1,1,1,nr+1),f,f(1,2),mt/2)
 
      end
