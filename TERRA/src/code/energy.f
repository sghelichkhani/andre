*dk energy
	subroutine energy(dhdt)
!	This routine updates the temperature field.
 
	include 'size.h'
	include 'pcom.h'

	common /flds/ f((nt+1)**2,nd,nr+1,8)
	common /temp/ tpb(nt+1), temp((nt+1)**2,nd,nr+1), tpe(nt+1)
	common /shht/ shb(nt+1),  shr((nt+1)**2,nd,nr+1), she(nt+1)
	common /volm/ vol((nt+1)**2,nr+1,2)
	common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
	common /heat/ htop, hbot, hrad, heat, hshr, hvol, hnet, tnrm,
     &              tav(nr+1), qc(nr)
	common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
	common /eos2/ rhocv(nr+1), gamma(nr+1), alpha(nr+1), prf(nr+1),
     &               cond(nr+1)
	common /volc/ tmpmlt
	common /call/ ncall
	common /clck/ itmng, sec(50)

	real dhdt((nt+1)**2,nd,nr+1)

	if(itmng==1) call mytime(tin)
	if(ncall==1) call facegen(f)
 
!	Compute the rate of temperature change due to advection and
!	conduction.
	call advect(dhdt)
 
!	Add volume heating contribution (radiogenic heating)
	if(hgen>0.0) then
		do ir=1,nr+1
			do id=1,nd
				do ii=1,(nt+1)**2
					dhdt(ii,id,ir) = dhdt(ii,id,ir) + hgen*rhorf(ir)
				enddo
			enddo
		enddo
	endif
 
!	Add shear heating contribution.
	if(ieos>=10) then
		do ir=1,nr+1
			do id=1,nd
				do ii=1,(nt+1)**2
					dhdt(ii,id,ir) = dhdt(ii,id,ir) + shr(ii,id,ir)
				enddo
			enddo
		enddo
	endif

	htop = 0.0
	hbot = 0.0
 
!	Add cooling due to melting and magma migration.
	if(tmpmlt.ne.0.0) call volcano(dhdt,temp)
 
!	Use heating rate of the nodal boundary layers as heat flux
!	leaving or entering the shell.
	do id=1,nd
 
		do ii=1,(nt+1)**2
			htop = htop + dhdt(ii,id,1)*vol(ii,1,1)
			dhdt(ii,id,1) = 0.0
		enddo
 
		if(tb(2)>0.0) then
			do ii=1,(nt+1)**2
				hbot = hbot - dhdt(ii,id,nr+1)*vol(ii,nr+1,1)
				dhdt(ii,id,nr+1) = 0.
			enddo
		endif

	enddo
 
	if(nproc>1) then
		call psum(htop,1)
		call psum(hbot,1)
	endif
 
!	Convert from heating rate to rate of change in temperature.
	do ir=2,nr+1
		a1 = 1.0/rhocv(ir)
		do id=1,nd
			do ii=1,(nt+1)**2
				dhdt(ii,id,ir) = a1*dhdt(ii,id,ir)
			enddo
		enddo
	enddo

	if(imoon==1.and.tb(2).ne.0.0) then
		do id=1,nd
			do ii=1,(nt+1)**2
				temp(ii,id,nr+1) = tb(2)
			enddo
		enddo
	endif
 
	if(itmng==1) then
		call mytime(tout)
		sec(20) = sec(20) + tout - tin
	endif

	end subroutine


*dk advect
	subroutine advect(dhdt)
!	This routine computes the heat advected and conducted per unit
!	time and also accounts for phase transitions if they are present.
 
	include 'size.h'

	common /velo/ upb(nt+1), u((nt+1)**2,nd,3,nr+1),  upe(nt+1)
	common /temp/ tpb(nt+1), temp((nt+1)**2,nd,nr+1), tpe(nt+1)
	common /mesh/ xn((nt+1)**2,nd,3)
	common /radl/ rshl(nr+1), ird
	common /volm/ vol((nt+1)**2,nr+1,2)
	common /face/ psi((nt+1)**2,3), rdxn((nt+1)**2,3),
     &              facet((nt+1)**2,nd,3,3)
	common /ndar/ arn((nt+1)**2),  arne((nt+1)**2),
     &              rarn((nt+1)**2), rarne((nt+1)**2)
	common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
	common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
	common /eos2/ rhocv(nr+1), gamma(nr+1), alpha(nr+1), prf(nr+1),
     &               cond(nr+1)
	common /qdfc/ j1f(4), j2f(4), jef(4), i1n(3), i2n(3)
	common /clck/ itmng, sec(50)

	real dhdt((nt+1)**2,nd,nr+1), dh((nt+1)**2,2)

	if(itmng==1) call mytime(tin)
 
	call nulvec(dhdt,(nt+1)**2*nd*(nr+1))

	! parameters for phase transitions 
	i410 = 0
	i660 = 0

	do ir=1,nr
 
		jr = ir+1
		rm = 0.5*(rshl(ir) + rshl(jr))
		r2 = rm**2
		
		if(jr<=nr) ra = r2 - 0.25*(rshl(jr) + rshl(jr+1))**2
		if(jr>nr) ra = r2 - rshl(nr+1)**2

		cr = 0.5*(cond(ir) + cond(jr))*rshl(ir)*rshl(jr)
     &          /(rshl(ir) - rshl(jr))
		ct = ra*cond(jr)/rshl(jr)
		ai = (gamma(ir) - 1.)*rhocv(ir)
		aj = (gamma(jr) - 1.)*rhocv(jr)

		if(i410==0.and.rshl(1)-rshl(jr)>410.e3) i410 = ir
		if(i660==0.and.rshl(1)-rshl(jr)>660.e3) i660 = ir

		do id=1,nd
!	Add contributions from radial advection and conduction,
!	adiabatic heating, and phase changes.
			ii1 = 2
			if(mod(id,5)==1) ii1 = 1

			do ii=ii1,(nt+1)**2
				uf1  = 0.5*(u(ii,id,1,ir) + u(ii,id,1,jr))
				uf2  = 0.5*(u(ii,id,2,ir) + u(ii,id,2,jr))
				uf3  = 0.5*(u(ii,id,3,ir) + u(ii,id,3,jr))
				vfl  = ((uf1*xn(ii,id,1) + uf2*xn(ii,id,2))
     &				+ uf3*xn(ii,id,3))*r2*arn(ii)
				tadv = 0.5*(temp(ii,id,ir) + temp(ii,id,jr))*vfl
				cnd  = (temp(ii,id,jr) - temp(ii,id,ir))*cr*arn(ii)
				dh(ii,1) = vfl
				dhdt(ii,id,ir) = (dhdt(ii,id,ir) + tadv*rhocv(ir)) + cnd
     &							+ temp(ii,id,ir)*vfl*ai
				dhdt(ii,id,jr) = (dhdt(ii,id,jr) - tadv*rhocv(jr)) - cnd
     &							- temp(ii,id,jr)*vfl*aj
			enddo
 
!	Add latent heat from phase transition if appropriate.
 			if(ieos>=10.and.(ir==i410.or.ir==i660)) then

				if(ir==i410) then
					cl =  cl410*230./3610.
					bb = 0.60*1.4e10*230./3610.
					if(ieos.ne.10) bb = 0.0
					at = (rshl(1) - 410.e3 - rshl(i410+1))
     &					/(rshl(i410) - rshl(i410+1))
				else
					cl =  cl660*380./4185.
					bb = 0.60*2.4e10*380./4185.
					if(ieos.ne.10) bb = 0.0
					at = (rshl(1) - 660.e3 - rshl(i660+1))
     &					/(rshl(i660) - rshl(i660+1))
				endif

				ab = 1. - at
				do ii=ii1,(nt+1)**2
					aa = (cl*(at*temp(ii,id,ir) + ab*temp(ii,id,jr))
     &				- bb)*dh(ii,1)
					dhdt(ii,id,ir) = dhdt(ii,id,ir) - aa*at
					dhdt(ii,id,jr) = dhdt(ii,id,jr) - aa*ab
				enddo

			endif
 
!	Add contributions from lateral advection and conduction.
			do ifc=1,3
		
				jj = i1n(ifc) + (nt+1)*i2n(ifc)
				ib = max(1, 1-jj)
				ie = min((nt+1)**2, (nt+1)**2-jj)

				do ii=ib,ie
					uf1  = 0.5*(u(ii,id,1,jr) + u(ii+jj,id,1,jr))
					uf2  = 0.5*(u(ii,id,2,jr) + u(ii+jj,id,2,jr))
					uf3  = 0.5*(u(ii,id,3,jr) + u(ii+jj,id,3,jr))
					vfl  = ((uf1*facet(ii,id,ifc,1)
     &					+ uf2*facet(ii,id,ifc,2))
     &					+ uf3*facet(ii,id,ifc,3))*psi(ii,ifc)*ra
					hadv = 0.5*(temp(ii,id,jr) + temp(ii+jj,id,jr))
     &					*vfl*rhocv(jr)
					cnd  = ct*(temp(ii,id,jr) - temp(ii+jj,id,jr))
     &					*psi(ii,ifc)*rdxn(ii,ifc)
					dh(ii,1) = hadv + cnd
					dh(ii,2) = vfl*aj
					dhdt(ii,id,jr) = dhdt(ii,id,jr) - dh(ii,1)
     &								- temp(ii,id,jr) * dh(ii,2)
				enddo

				do ii=ib,ie
					dhdt(ii+jj,id,jr) = dhdt(ii+jj,id,jr) + dh(ii,1)
     &									+ temp(ii+jj,id,jr) * dh(ii,2)
				enddo
			enddo
		enddo
	enddo
	 
	call comm3s(dhdt,nr,nt,1)
 
!	Divide by nodal volume.
	do ir=1,nr+1
		do id=1,nd
			do ii=1,(nt+1)**2
				dhdt(ii,id,ir) = dhdt(ii,id,ir)*vol(ii,ir,2)
			enddo
		enddo
	enddo
 
	if(itmng==1) then
		call mytime(tout)
		sec(21) = sec(21) + tout - tin
	endif

	end subroutine


*dk facegen
      subroutine facegen(xc)
 
c...  This routine computes the arrays facet, psi, and rdxn.  Array
c...  facet contains the unit normals to the quadrilateral faces; psi
c...  contains the half-angle subtended by the faces with respect to the
c...  center of the sphere, and rdxn contains the reciprocal distances
c...  between the nodes on opposite sides of the faces.  Three faces
c...  are associated with each node, and these are numbered in a
c...  counter-clockwise order for diamonds linked with the north pole.
c...  Since these quantities are identical for all faces in a given
c...  radial stack, they are computed only for one layer.
 
      include 'size.h'
      include 'pcom.h'
      real xc(nt+1,0:nt,2,3)
      common /xnex/ xne(-1:nt+1,0:nt+2,nd,3)
      common /face/   psi(0:nt,nt+1,3), rdxn(0:nt,nt+1,3),
     &              facet(0:nt,nt+1,nd,3,3)
      common /qdfc/ j1f(4), j2f(4), jef(4), i1n(3), i2n(3)
 
      call nulvec(facet,(nt+1)**2*nd*9)
      call nulvec(psi  ,(nt+1)**2*3)
      call nulvec(rdxn ,(nt+1)**2*3)
 
      id0 = (mynum/mproc)*nd
 
      do id=1,nd
 
      jd  = id + id0
 
      call tricntr(xc,xne,id)
 
      do ifc=1,3
 
         je = jef(ifc)
         ke = jef(ifc+1)
 
         do i2=1,nt
 
            j2 = i2 + j2f(ifc)
            k2 = i2 + j2f(ifc+1)
            m2 = i2 + i2n(ifc)
 
            do i1=1,nt
 
               j1 = i1 + j1f(ifc)
               k1 = i1 + j1f(ifc+1)
               m1 = i1 + i1n(ifc)
 
               if(id .eq. 1) then
 
c...              Compute the reciprocal distance between nodes.
 
                  aa = (xne(i1,i2,1,1) - xne(m1,m2,1,1))**2
     &               + (xne(i1,i2,1,2) - xne(m1,m2,1,2))**2
     &               + (xne(i1,i2,1,3) - xne(m1,m2,1,3))**2
                  rdxn(i1,i2,ifc) = 1./sqrt(aa)
 
c...              Compute the half-angle psi between the two points
c...              defining the face with respect to the center of
c...              the sphere.
 
                  aa = (xc(k1,k2,ke,1) - xc(j1,j2,je,1))**2
     &               + (xc(k1,k2,ke,2) - xc(j1,j2,je,2))**2
     &               + (xc(k1,k2,ke,3) - xc(j1,j2,je,3))**2
                  psi(i1,i2,ifc) = asin(0.5*sqrt(aa))
 
               endif
 
c...           The array facet is given by the cross product of vectors
c...           from the origin to the two points defining the face,
c...           normalized to unity and pointing outward from the center
c...           of the cell.
 
               aa = 1.
               if(jd .ge. 6) aa = -1.
               aa = aa/sin(2.*psi(i1,i2,ifc))
 
               facet(i1,i2,id,ifc,1) =
     &                             aa*(xc(k1,k2,ke,2)*xc(j1,j2,je,3)
     &                               - xc(k1,k2,ke,3)*xc(j1,j2,je,2))
               facet(i1,i2,id,ifc,2) =
     &                             aa*(xc(k1,k2,ke,3)*xc(j1,j2,je,1)
     &                               - xc(k1,k2,ke,1)*xc(j1,j2,je,3))
               facet(i1,i2,id,ifc,3) =
     &                             aa*(xc(k1,k2,ke,1)*xc(j1,j2,je,2)
     &                               - xc(k1,k2,ke,2)*xc(j1,j2,je,1))
 
            end do
 
         end do
 
      end do
 
      end do
 
      end subroutine
      
      
*dk heating
      subroutine heating(heat,hrad,tdot)
 
c...  This routine computes the total heating rate from the time rate
c...  of change of temperature tdot as well as the radiogenic heating
c...  rate.
 
      include 'size.h'
      include 'pcom.h'
      real tdot((nt+1)**2,nd,nr+1)
      common /volm/ vol((nt+1)**2,nr+1,2)
      common /radl/ rshl(nr+1), ird
      common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
      common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
      common /eos2/ rhocv(nr+1), gamma(nr+1), alpha(nr+1), prf(nr+1),
     &               cond(nr+1)
 
      heat = 0.
      davg = 0.
      svol = 0.
 
      do ir=1,nr+1
 
         do id=1,nd
            do ii=1,(nt+1)**2
               heat = heat + rhocv(ir)*tdot(ii,id,ir)*vol(ii,ir,1)
            end do
         end do
 
         davg = davg + vol(2,ir,1)*rhorf(ir)
         svol = svol + vol(2,ir,1)
 
      end do
 
      hrad = 4.18879*hgen*davg/svol*(rshl(1)**3 - rshl(nr+1)**3)
 
      if(nproc .gt. 1) call psum(heat,1)
 
      end
            
*dk shrheat
	subroutine shrheat(dhdt,u,z,g)
!	This routine computes the shear heating rate.
 
      include 'size.h'
      include 'pcom.h'
      include 'para.h'
      
      common /grad/ rg(3,2,nr+1), grd(7,(nt+1)**2,3,2)
      common /ndar/  arn((nt+1)**2),  arne((nt+1)**2),
     &              rarn((nt+1)**2), rarne((nt+1)**2)
      common /volm/ vol((nt+1)**2,nr+1,2)
      common /ofst/ j1n(7), j2n(7), md(7)
      common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
      common /heat/ htop, hbot, hrad, heat, hshr, hvol, hnet, tnrm,
     &              tav(nr+1), qc(nr)
      common /eos2/ rhocv(nr+1), gamma(nr+1), alpha(nr+1), prf(nr+1),
     &               cond(nr+1)
      common /vis0/ vb(nt+1), vv((nt+1)**2*nd*(nr+1),2), ve(nt+1)
      common /vis1/ vscmax, rvscscl, tvscscl, pwrlawn, pwrlawsr, yldstrs
      common /shrh/ shrlayr(nr+1)
	  common /radl/ rshl(nr+1), ird
      common /call/ ncall
      common /clck/ itmng, sec(50)

	integer shearcut
	real dhdt((nt+1)**2,nd,*), u((nt+1)**2,nd,3,*)
	real g((nt+1)**2,nd,3), z((nt+1)**2,nd,3,3), rr(3,2,nr+1)
      
	if(itmng==1) call mytime(tin)
 
	twothd = 2.0/3.0
 	call nulvec(dhdt,(nt+1)**2*nd*(nr+1))
 
	aa = sqrt(visc)
	do k=1,3
		do ir=1,nr+1
			rr(k,1,ir) = aa*rg(k,1,ir)
			rr(k,2,ir) = aa*rg(k,2,ir)
		enddo
	enddo
	
!	compute the shear heating
	do ir=1,nr+1
		if((rshl(1)-rshl(ir))/1000.0<=dpth_shcut) shearcut=ir
	enddo
	shearcut=min(shearcut+1,nr+1)
	do ir=shearcut,nr+1

 !	Compute the node-to-node gradient of velocity.
		call nulvec(z,(nt+1)**2*nd*9)
 
		k1 = 1
		k2 = 3
		if(ir==1) k1 = 2
		if(ir==nr+1) k2 = 2
 
		do k=k1,k2
			do m=1,7
 				jj  = j1n(m) + (nt+1)*j2n(m)
				jr = ir + k - 2
 
			      do ii=1,(nt+1)**2
					g(ii,1,1) = vol(ii,ir,2)*(grd(m,ii,1,1)*rr(k,1,ir)
     &                        + grd(m,ii,1,2)*rr(k,2,ir))
					g(ii,1,2) = vol(ii,ir,2)*(grd(m,ii,2,1)*rr(k,1,ir)
     &                        + grd(m,ii,2,2)*rr(k,2,ir))
					g(ii,1,3) = vol(ii,ir,2)*(grd(m,ii,3,1)*rr(k,1,ir)
     &                        + grd(m,ii,3,2)*rr(k,2,ir))
				enddo
				
				do id=2,nd
					do ii=1,(nt+1)**2
						g(ii,id,1) = g(ii,1,1)
						g(ii,id,2) = g(ii,1,2)
						g(ii,id,3) = g(ii,1,3)
					enddo
				enddo

				call rotate(g,nd,0,nt,-1)
 
				do id=1,nd
 					do ii=1,(nt+1)**2
						z(ii,id,1,1) = z(ii,id,1,1) + g(ii,id,1)*u(ii+jj,id,1,jr)
						z(ii,id,2,1) = z(ii,id,2,1) + g(ii,id,2)*u(ii+jj,id,1,jr)
						z(ii,id,3,1) = z(ii,id,3,1) + g(ii,id,3)*u(ii+jj,id,1,jr)
						z(ii,id,1,2) = z(ii,id,1,2) + g(ii,id,1)*u(ii+jj,id,2,jr)
						z(ii,id,2,2) = z(ii,id,2,2) + g(ii,id,2)*u(ii+jj,id,2,jr)
						z(ii,id,3,2) = z(ii,id,3,2) + g(ii,id,3)*u(ii+jj,id,2,jr)
			 		enddo
 					
 					do ii=1,(nt+1)**2
						z(ii,id,1,3) = z(ii,id,1,3) + g(ii,id,1)*u(ii+jj,id,3,jr)
						z(ii,id,2,3) = z(ii,id,2,3) + g(ii,id,2)*u(ii+jj,id,3,jr)
						z(ii,id,3,3) = z(ii,id,3,3) + g(ii,id,3)*u(ii+jj,id,3,jr)
					enddo
				enddo
			enddo
		enddo

		call comm3s(z,2,nt,3)
 
		do id=1,nd
			do ii=1,(nt+1)**2
				dhdt(ii,id,ir) = (((2.*((z(ii,id,1,1)**2
     &			   +  z(ii,id,2,2)**2) + z(ii,id,3,3)**2)
     &			   + (z(ii,id,1,2) + z(ii,id,2,1))**2)
     &			   + (z(ii,id,2,3) + z(ii,id,3,2))**2)
     &			   + (z(ii,id,3,1) + z(ii,id,1,3))**2)
     & 			   - (z(ii,id,1,1) + z(ii,id,2,2) + z(ii,id,3,3))**2*twothd
			enddo
		enddo
	enddo
 
	if(vscmax>0.0) then
		do ir=1,nr+1
			do id=1,nd
				jj = (nd*(ir - 1) + (id - 1))*(nt + 1)**2
				do ii=1,(nt+1)**2
					dhdt(ii,id,ir) = dhdt(ii,id,ir)*vv(ii+jj,1)**2
				enddo
			enddo
		enddo
	endif
 
	hshr = 0.0
	do ir=1,nr+1
		shrlayr(ir) = 0.
		do id=1,nd
			do ii=1,(nt+1)**2
				dh = dhdt(ii,id,ir)*vol(ii,ir,1)
				hshr = hshr + dh
				shrlayr(ir) = shrlayr(ir) + dh
			enddo
		enddo
		shrlayr(ir) = shrlayr(ir)*arn(2)/(4.*3.14159265*vol(2,ir,1))
	enddo
	 
      if(nproc>1) then
      	call psum(hshr,1)
      	call psum(shrlayr,nr+1)
	endif
	
      if(itmng==1) then
      	call mytime(tout)
      	sec(22) = sec(22) + 0.25*(tout - tin)
      endif
      
	end subroutine
	
	
*dk volcano
      subroutine volcano(dhdt,temp)

! VORSICHT!!! TSTEP NICHT DEFINIERT!!!!!!!!!

c...  This routine adds cooling to reduce the temperature to the
c...  melting temperature where this temperature is exceeded.  This
c...  is intended to model melting and rapid heat removal by magma
c...  migration to the surface.
 
      include 'size.h'
      include 'pcom.h'
      real dhdt((nt+1)**2,nd,nr+1), temp((nt+1)**2,nd,nr+1)
      common /volm/ vol((nt+1)**2,nr+1,2)
      common /radl/ rshl(nr+1), ird
      common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
      common /heat/ htop, hbot, hrad, heat, hshr, hvol, hnet, tnrm,
     &              tav(nr+1), qc(nr)
      common /eos2/ rhocv(nr+1), gamma(nr+1), alpha(nr+1), prf(nr+1),
     &               cond(nr+1)
      common /volc/ tmpmlt
 
      if(tmpmlt .eq. 0.) return
 
      hvol = 0.
 
      irmx = nr/3
      if(tb(2).eq.0. .and. ieos.eq.1) irmx = nr
 
      do ir=2,irmx
 
         aa   = 0.
         if(tstep .gt. 0.) aa = rhocv(ir)/tstep
         tmlt = tmpmlt + 1.8e-3*(rshl(1) - rshl(ir))
         if(ir .gt. nr/4) tmlt = min(tmlt, 3000.)
 
         do id=1,nd
            do ii=1,(nt+1)**2
               cool = max(0., aa*(temp(ii,id,ir) - tmlt))
               hvol = hvol + cool*vol(ii,ir,1)
               dhdt(ii,id,ir) = dhdt(ii,id,ir) - cool
            end do
         end do
 
      end do
 
      htop = htop + hvol
 
      if(nproc .gt. 1) call psum(hvol, 1)
 
      end
*dk bdenergy
      blockdata bdenergy
 
      common /qdfc/ j1f(4), j2f(4), jef(4), i1n(3), i2n(3)
      data j1f/1,0,0,0/, j2f/0,0,0,-1/, jef/1,2,1,2/
      data i1n/0,-1,-1/, i2n/1,1,0/
 
      end blockdata
