*dk usolve
	subroutine usolve(map,urot,ireport)
!	This routine solves for the new velocity field u.

	include 'size.h'
	include 'pcom.h'
	include 'para.h'

	common /flds/ f((nt+1)**2*nd,3,nr+1), r((nt+1)**2*nd,nr+1),
     &              v((nt+1)**2*nd,3,nr+1), s((nt+1)**2*nd,nr+1)
	common /mgwk/ w((nt+1)**2*nd,nr+1,5)
	common /velo/ upb(nt+1), u((nt+1)**2*nd,3,nr+1),  upe(nt+1)
	common /pres/ ppb(nt+1), pres((nt+1)**2*nd,nr+1), ppe(nt+1)
	common /mesh/ xn((nt+1)**2*nd,3)
	common /radl/ rshl(nr+1), ird
	common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
	common /vis1/ vscmax, rvscscl, tvscscl, pwrlawn, pwrlawsr, yldstrs
	common /vis2/ rdvsc(nr+1), tactv(nr+1), vscl(nr+1)
	common /nrms/ fnrm, rnrm, unrm, ekin
	common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
	common /mgrd/ itlimit, convtol, itsolve
	common /opwt/ wd, wm
	common /solv/ itreq, convrg, divnorm, cutdiv
	common /io01/ casenum, gpath, lpath, opath, tpath
	common /clck/ itmng, sec(50)
	common /call/ ncall

	integer itr, ibc0, izerou, ireport
	integer map(0:nt,nt+1,nd)
	
	real divnorm, cutdiv
	real urot(3,pl_size)
	character char1*4, char2*4, cname*14, casenum*3
	character lpath*80, opath*80
	
	if(itmng==1) call mytime(tin)
 
!	Compute the nodal force field f.
!	and substitude the first layer of the velocity field u
!	by the given plate velocity
	call forces(map,urot)
 
!	Solve for an estimate of the new velocity field.
	wd = 1.0
	wm	= 0.0
 
	if(vscmax>0.0) call uscale(u,1,nr)

	! calculates u from f (Au=f)
	! input u here serves as a first guess
	call multigrid(u,f,3,map,urot,ibc,0,ireport)

	if(vscmax>0.0) call uscale(u,2,nr)
 
	velnrm  = unrm
	fnrm0   = fnrm
	rnrm0   = rnrm
	convrg  = rnrm/fnrm
	cutdiv  = 8.0*(rdvsc(1)**0.25)*convtol*unrm/rshl(1)
	! itsolve: number of multigrid iterations needed
	! has influence on the next time step
	itreq   = itsolve

	itr=1
	izerou=1
	divnorm=cutdiv
	
	!if(mynum==0) write(6,*) "cutdiv: ", cutdiv
	! pressure correction
	do while(itr<=npres.and.divnorm>=cutdiv)

		if(itr==1) then
 
			call divergence(r,w,u,w(1,1,2),w(1,1,5),w,1)

			call norm3s(r,divnorm,1,nd,nr,nt)
 
			if(divnorm>=cutdiv) then
				do ir=1,nr+1
					do ii=1,(nt+1)**2*nd
						s(ii,ir) = r(ii,ir)
					enddo
				enddo
			endif
			
		else
		
			delta = (divnorm/divnorm0)**2.
			do ir=1,nr+1
				do ii=1,(nt+1)**2*nd
					s(ii,ir) = r(ii,ir) + delta*s(ii,ir)
				enddo
			enddo
 
		endif
 
		if(divnorm>=cutdiv) then

			call gradient(f,s)
 
			if(vscmax>0.0) call uscale(f,2,nr)
 
			do ir=1,nr+1
				a1 = rhorf(ir)/rho0
				do ii=1,(nt+1)**2*nd
					f(ii,1,ir) = a1*f(ii,1,ir)
					f(ii,2,ir) = a1*f(ii,2,ir)
					f(ii,3,ir) = a1*f(ii,3,ir)
				enddo
			enddo
 
			call nulvec(v, (nt+1)**2*3*nd*(nr+1))
 
			if(ibc==6) then
				ibc0=5
			else
				ibc0=ibc
			endif

			call multigrid(v,f,3,map,urot,ibc0,izerou,ireport)
			izerou=0

			if(vscmax>0.0) call uscale(v,2,nr)			

			call divergence(w(1,1,2),w,v,w(1,1,2),w(1,1,5),w,1)
 
			call dotprod(a0,w(1,1,2),s)
 
			if(a0.ne.0.0) a0 = -divnorm**2/a0
			a1 = a0*visc
 
			do ir=1,nr+1
				do ii=1,(nt+1)**2*nd
					r(ii,ir)    = r(ii,ir)    + a0*w(ii,ir,2)
					pres(ii,ir) = pres(ii,ir) + a1*s(ii,ir)
					u(ii,1,ir)  = u(ii,1,ir)  + a0*v(ii,1,ir)
					u(ii,2,ir)  = u(ii,2,ir)  + a0*v(ii,2,ir)
					u(ii,3,ir)  = u(ii,3,ir)  + a0*v(ii,3,ir)
				enddo
			enddo

			divnorm0 = divnorm
 
			call norm3s(r,divnorm,1,nd,nr,nt)
 
		endif
		
		itr=itr+1
	enddo
 
 	!if(mod(ncall,4)==1) then
 	!	write(char1,'(I4.4)') mynum
 	!	write(char2,'(I4.4)') ncall
 	!	cname=	'g'//casenum//'.'//char1//'.'//char2
	!	open(62, file=lpath(1:nlpath)//cname, status='unknown')
	!	call vecout(r,rshl,1,nr,nt,nr,62,0)
	!	close(62)
	!endif
			
	!if(mynum==0) write(6,*) "pres: ", itr-1, divnorm
	unrm = velnrm
	fnrm = fnrm0
	rnrm = rnrm0
 
	if(itmng==1) then
		call mytime(tout)
		sec(5) = sec(5) + tout - tin
	endif
      
	end subroutine
      

*dk forces
	subroutine forces(map,urot)
 
	include 'size.h'
	include 'para.h'

	common /flds/ f((nt+1)**2*nd,3,nr+1), b((nt+1)**2*nd,nr+1),
     &              w((nt+1)**2*nd,4,nr+1)
	common /velo/ upb(nt+1), u((nt+1)**2*nd,3,nr+1),  upe(nt+1)
	common /pres/ ppb(nt+1), pres((nt+1)**2*nd,nr+1), ppe(nt+1)
	common /mesh/ xn
	common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
	common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
	common /vis1/ vscmax, rvscscl, tvscscl, pwrlawn, pwrlawsr, yldstrs               
	common /clck/ itmng, sec(50)

	integer map(0:nt,nt+1,nd)
	real xn((nt+1)**2*nd,3)
	real urot(3,pl_size)

	if(itmng==1) call mytime(tin)

	! calculates b (weak form.) AND w (nodal)
	call buoyancy(b,w)
	
	! calculates f from pres 
	call gradient(f,pres)
 
	do ir=1,nr+1
		aa = rhorf(ir)/rho0
		! total force = buoyancy force + pressure gradient
		do ii=1,(nt+1)**2*nd
			f(ii,1,ir) = aa*f(ii,1,ir) - b(ii,ir)*xn(ii,1)
			f(ii,2,ir) = aa*f(ii,2,ir) - b(ii,ir)*xn(ii,2)
			f(ii,3,ir) = aa*f(ii,3,ir) - b(ii,ir)*xn(ii,3)
		enddo
	enddo

	! f is changed here due to phase changes
	call fphase
 
	if(vscmax>0.0) call uscale(f,2,nr)
 
	aa = 1.0/visc
	do ir=1,nr+1
		do ii=1,(nt+1)**2*nd
			f(ii,1,ir) = aa*f(ii,1,ir)
			f(ii,2,ir) = aa*f(ii,2,ir)
			f(ii,3,ir) = aa*f(ii,3,ir)
		enddo
	enddo

!	at this point, we have f, w and b (just calculated)
!	and u from last time-step (time 1: u=0)
!	now, we replace the first layer of the velocity field
!	by the given plate velocity
	if(ibc==6) call platevelreplace(u,urot,xn,map)
 
	if(itmng==1) then
		call mytime(tout)
		sec(15) = sec(15) + tout - tin
	endif
      
	end subroutine


*dk buoyancy
	subroutine buoyancy(buoy,w)
!	This routine computes the scalar buoyancy field arising from
!	temperature and pressure variations.

	include 'size.h'

	real buoy(*), w((nt+1)**2*nd,nr+1)
	common /temp/ tpb(nt+1), temp((nt+1)**2*nd,nr+1), tpe(nt+1)
	common /pres/ ppb(nt+1), pres((nt+1)**2*nd,nr+1), ppe(nt+1)
	common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
	common /heat/ htop, hbot, hrad, heat, hshr, hvol, hnet, tnrm,
     &              tav(nr+1), qc(nr)
	common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
	common /eos2/ rhocv(nr+1), gamma(nr+1), alpha(nr+1), prf(nr+1),
     &               cond(nr+1)
 
	do ir=1,nr+1
 
		a1 = alpha(ir)*rhorf(ir)*grv(ir)
		a2 = -drhdp(ir)*grv(ir)*rhorf(ir)/rho0
		tref = tmprf(ir)
		if(ieos/10==0) tref = tnrm

		do ii=1,(nt+1)**2*nd
			w(ii,ir) = a1*(temp(ii,ir) - tref) + a2*pres(ii,ir)
		enddo
 
	enddo
 	
 	! buoy = M*w
	call massmtrx(buoy,w)
 
	end subroutine


*dk divergence
	subroutine divergence(d,t,u,v,q,r,iflag)
!	This routine computes the divergence d of various possible
c...  scalar products with the velocity field u, depending on the
c...  value of the flag iflag as follows:
c...     When iflag = 0, the routine simply computes div(u).
c...     When iflag = 1, the routine computes
c...                     div[(rhorf/rho)*u]/(rhorf/rho)**2
c...     When iflag = 2, the routine computes div(t*u).
c...  The second case involves weighting the velocity by the normalized
c...  radial density variation needed in the anelastic approximation.
c...  The last case corresponds to advection of the scalar field t by
c...  the velocity field u.
 
c...  The routine applies the transposed of the finite element gradient
c...  operator to obtain the divergence.  Since the finite element
c...  operator contains a factor of volume, it is necessary to solve
c...  the linear system  M*d = D*(s*u) for d, where M is the finite
c...  element mass matrix, D is the transpose of the gradient, and s
c...  is the appropriate field multiplying the velocity u.  Jacobi
c...  iteration is used to solve this well-conditioned system.  The
c...  arrays v, q, and r are work arrays.
 
      include 'size.h'
	include 'para.h'

      real d(0:nt,nt+1,nd,nr+1), v(0:nt,nt+1,nd,3,nr+1)
      real t(0:nt,nt+1,nd,nr+1), r(0:nt,nt+1,nd,  nr+1)
      real q(0:nt,nt+1,nd,nr+1), w(2,0:nt,nr+1), u(*)
      common /grad/ rg(3,2,nr+1), grd(7,0:nt,nt+1,3,2)
      common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
      common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
      common /volm/ vol(0:nt,nt+1,nr+1,2)
      common /radl/ rshl(nr+1), ird
      common /clck/ itmng, sec(50)
      if(itmng.eq.1) call mytime(tin)
 
      call scopy((nt+1)**2*(nr+1)*3*nd,u,1,v,1)
      call rotate(v,nd,nr,nt,1)
 
      if(iflag==1) then
 
         do ir=1,nr+1
            a1 = rhorf(ir)/rho0
            do id=1,nd
               do i2=1,nt+1
                  do i1=0,nt
                     v(i1,i2,id,1,ir) = a1*v(i1,i2,id,1,ir)
                     v(i1,i2,id,2,ir) = a1*v(i1,i2,id,2,ir)
                     v(i1,i2,id,3,ir) = a1*v(i1,i2,id,3,ir)
                  end do
               end do
            end do
         end do
 
      elseif(iflag==2) then
 
         do ir=1,nr+1
         	do id=1,nd
         	   do i2=1,nt+1
         	      do i1=0,nt
         	         v(i1,i2,id,1,ir) = t(i1,i2,id,ir)*v(i1,i2,id,1,ir)
         	         v(i1,i2,id,2,ir) = t(i1,i2,id,ir)*v(i1,i2,id,2,ir)
         	         v(i1,i2,id,3,ir) = t(i1,i2,id,ir)*v(i1,i2,id,3,ir)
         			enddo
               enddo
            enddo
         enddo
 
      elseif(iflag .ne. 0) then
 
         print *, 'STOP--invalid value for iflag in routine divergence.'
         stop
 
      end if
 
      do id=1,nd
 
         do i2=1,nt+1
 
            do ir=1,nr+1
               do i1=0,nt
                  w(1,i1,ir) = 0.
                  w(2,i1,ir) = 0.
               end do
            end do
 
            do j=1,3
               do ir=1,nr,2
                  do i1=0,nt
                     w(1,i1,ir) = w(1,i1,ir) + (((((((
     &                 grd(1,i1,i2,j,1)*v(i1  ,i2  ,id,j,ir))
     &               + grd(2,i1,i2,j,1)*v(i1+1,i2  ,id,j,ir))
     &               + grd(3,i1,i2,j,1)*v(i1  ,i2+1,id,j,ir))
     &               + grd(4,i1,i2,j,1)*v(i1-1,i2+1,id,j,ir))
     &               + grd(5,i1,i2,j,1)*v(i1-1,i2  ,id,j,ir))
     &               + grd(6,i1,i2,j,1)*v(i1  ,i2-1,id,j,ir))
     &               + grd(7,i1,i2,j,1)*v(i1+1,i2-1,id,j,ir))
                     w(1,i1,ir+1) = w(1,i1,ir+1) + (((((((
     &                 grd(1,i1,i2,j,1)*v(i1  ,i2  ,id,j,ir+1))
     &               + grd(2,i1,i2,j,1)*v(i1+1,i2  ,id,j,ir+1))
     &               + grd(3,i1,i2,j,1)*v(i1  ,i2+1,id,j,ir+1))
     &               + grd(4,i1,i2,j,1)*v(i1-1,i2+1,id,j,ir+1))
     &               + grd(5,i1,i2,j,1)*v(i1-1,i2  ,id,j,ir+1))
     &               + grd(6,i1,i2,j,1)*v(i1  ,i2-1,id,j,ir+1))
     &               + grd(7,i1,i2,j,1)*v(i1+1,i2-1,id,j,ir+1))
                     w(2,i1,ir) = w(2,i1,ir) + (((((((
     &                 grd(1,i1,i2,j,2)*v(i1  ,i2  ,id,j,ir))
     &               + grd(2,i1,i2,j,2)*v(i1+1,i2  ,id,j,ir))
     &               + grd(3,i1,i2,j,2)*v(i1  ,i2+1,id,j,ir))
     &               + grd(4,i1,i2,j,2)*v(i1-1,i2+1,id,j,ir))
     &               + grd(5,i1,i2,j,2)*v(i1-1,i2  ,id,j,ir))
     &               + grd(6,i1,i2,j,2)*v(i1  ,i2-1,id,j,ir))
     &               + grd(7,i1,i2,j,2)*v(i1+1,i2-1,id,j,ir))
                     w(2,i1,ir+1) = w(2,i1,ir+1) + (((((((
     &                 grd(1,i1,i2,j,2)*v(i1  ,i2  ,id,j,ir+1))
     &               + grd(2,i1,i2,j,2)*v(i1+1,i2  ,id,j,ir+1))
     &               + grd(3,i1,i2,j,2)*v(i1  ,i2+1,id,j,ir+1))
     &               + grd(4,i1,i2,j,2)*v(i1-1,i2+1,id,j,ir+1))
     &               + grd(5,i1,i2,j,2)*v(i1-1,i2  ,id,j,ir+1))
     &               + grd(6,i1,i2,j,2)*v(i1  ,i2-1,id,j,ir+1))
     &               + grd(7,i1,i2,j,2)*v(i1+1,i2-1,id,j,ir+1))
                  end do
               end do
                  do i1=0,nt
                     w(1,i1,nr+1) = w(1,i1,nr+1) + (((((((
     &                 grd(1,i1,i2,j,1)*v(i1  ,i2  ,id,j,nr+1))
     &               + grd(2,i1,i2,j,1)*v(i1+1,i2  ,id,j,nr+1))
     &               + grd(3,i1,i2,j,1)*v(i1  ,i2+1,id,j,nr+1))
     &               + grd(4,i1,i2,j,1)*v(i1-1,i2+1,id,j,nr+1))
     &               + grd(5,i1,i2,j,1)*v(i1-1,i2  ,id,j,nr+1))
     &               + grd(6,i1,i2,j,1)*v(i1  ,i2-1,id,j,nr+1))
     &               + grd(7,i1,i2,j,1)*v(i1+1,i2-1,id,j,nr+1))
                     w(2,i1,nr+1) = w(2,i1,nr+1) + (((((((
     &                 grd(1,i1,i2,j,2)*v(i1  ,i2  ,id,j,nr+1))
     &               + grd(2,i1,i2,j,2)*v(i1+1,i2  ,id,j,nr+1))
     &               + grd(3,i1,i2,j,2)*v(i1  ,i2+1,id,j,nr+1))
     &               + grd(4,i1,i2,j,2)*v(i1-1,i2+1,id,j,nr+1))
     &               + grd(5,i1,i2,j,2)*v(i1-1,i2  ,id,j,nr+1))
     &               + grd(6,i1,i2,j,2)*v(i1  ,i2-1,id,j,nr+1))
     &               + grd(7,i1,i2,j,2)*v(i1+1,i2-1,id,j,nr+1))
                  end do
            end do
 
            do i1=0,nt
               q(i1,i2,id,1) =
     &          (((rg(2,1,1)*w(1,i1,1)
     &           + rg(2,2,1)*w(2,i1,1))
     &           + rg(3,1,1)*w(1,i1,2))
     &           + rg(3,2,1)*w(2,i1,2))
            end do
 
            do ir=2,nr
               do i1=0,nt
                  q(i1,i2,id,ir) =
     &             (((((rg(1,1,ir)*w(1,i1,ir-1)
     &                + rg(1,2,ir)*w(2,i1,ir-1))
     &                + rg(2,1,ir)*w(1,i1,ir  ))
     &                + rg(2,2,ir)*w(2,i1,ir  ))
     &                + rg(3,1,ir)*w(1,i1,ir+1))
     &                + rg(3,2,ir)*w(2,i1,ir+1))
               end do
            end do
 
            do i1=0,nt
               q(i1,i2,id,nr+1) =
     &          (((rg(1,1,nr+1)*w(1,i1,nr)
     &           + rg(1,2,nr+1)*w(2,i1,nr))
     &           + rg(2,1,nr+1)*w(1,i1,nr+1))
     &           + rg(2,2,nr+1)*w(2,i1,nr+1))
            end do
 
         end do
 
      end do
 
      call comm3s(q,nr,nt,1)
 
      do ir=1,nr+1
         do id=1,nd
            do i2=1,nt+1
               do i1=0,nt
                  d(i1,i2,id,ir) = q(i1,i2,id,ir)*vol(i1,i2,ir,2)
               end do
            end do
         end do
      end do
 
      itermax = 3
      if(iflag==1.and.ibc==6) itermax = 0
 
      do it=1,itermax
 
         call massmtrx(r,d)
 
         do ir=1,nr+1
            do id=1,nd
               do i2=1,nt+1
                  do i1=0,nt
                     r(i1,i2,id,ir) = q(i1,i2,id,ir) - r(i1,i2,id,ir)
                     d(i1,i2,id,ir) = d(i1,i2,id,ir) + r(i1,i2,id,ir)
     &                                          *1.3*vol(i1,i2,ir,2)
                  end do
               end do
            end do
         end do
 
      end do
 
      if(iflag .eq. 1) then
 
         do ir=1,nr+1
            a1 = (rho0/rhorf(ir))**2
            do id=1,nd
               do i2=1,nt+1
                  do i1=0,nt
                     d(i1,i2,id,ir) = a1*d(i1,i2,id,ir)
                  end do
               end do
            end do
         end do
 
      end if
 
      if(itmng.eq.1) call mytime(tout)
      if(itmng.eq.1) sec(17) = sec(17) + tout - tin
      
	end subroutine


*dk dotprod
	subroutine dotprod(dotprd,u,v)
!	This routine computes the inner product of the fields u and v.
 
      include 'size.h'
      include 'pcom.h'
      real u((nt+1)**2,nd,nr+1), v((nt+1)**2,nd,nr+1)
 
      sum = 0.
 
      do ir=1,nr+1
 
         if(mynum .eq. 0) sum = sum + u(1,1,ir)*v(1,1,ir)
         if(mynum .eq. 0 .and. nd .eq. 10)
     &                    sum = sum + u(1,6,ir)*v(1,6,ir)
         if(mynum .eq. mproc .and. nd .le. 5)
     &                    sum = sum + u(1,1,ir)*v(1,1,ir)
 
         do id=1,nd
 
            do ii=2,nt*(nt+1)
               sum = sum + u(ii,id,ir)*v(ii,id,ir)
            end do
 
            do ii=nt+2,nt*(nt+1),nt+1
               sum = sum - u(ii,id,ir)*v(ii,id,ir)
            end do
 
         end do
 
      end do
 
      if(nproc .gt. 1) call psum(sum,1)
 
      dotprd = sum/((10*mt*mt+2)*(nr+1))

	end subroutine


*dk fphase
	subroutine fphase
!	This routine adds a radial force at nodes near the phase
!	boundaries to account for the depth perturbation of the 410
!	and 660 km phase transitions.
 
      include 'size.h'
      common /flds/ f((nt+1)**2,nd,3,nr+1), w((nt+1)**2,nd,nr+1,5)
      common /temp/ tpb(nt+1), temp((nt+1)**2,nd,nr+1), tpe(nt+1)
      common /mesh/ xn((nt+1)**2,nd,3)
      common /radl/ rshl(nr+1), ird
      common /ndar/  arn((nt+1)**2),  arne((nt+1)**2),
     &              rarn((nt+1)**2), rarne((nt+1)**2)
      common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
      common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
 
      if(cl410.eq.0. .and. cl660.eq.0.) return
 
      ir410 = 0
      ir660 = 0
      do ir=2,nr
         if(ir410.eq.0 .and. rshl(1)-rshl(ir).gt.410.e3) ir410 = ir - 1
         if(ir660.eq.0 .and. rshl(1)-rshl(ir).gt.660.e3) ir660 = ir - 1
      end do
 
      at410 = (rshl(1) - 410.e3 - rshl(ir410+1))
     &            /(rshl(ir410) - rshl(ir410+1))
      at660 = (rshl(1) - 660.e3 - rshl(ir660+1))
     &            /(rshl(ir660) - rshl(ir660+1))
      ab410 = 1. - at410
      ab660 = 1. - at660
      cr410 = cl410*230./3610.*(rshl(1) - 410.e3)**2
      cr660 = cl660*380./4185.*(rshl(1) - 660.e3)**2
 
      do id=1,nd
         do ii=1,(nt+1)**2
            aa = cr410*arne(ii)
     &         *(at410*(temp(ii,id,ir410)   - tmprf(ir410))
     &         + ab410*(temp(ii,id,ir410+1) - tmprf(ir410+1)))
            at = at410*aa
            ab = ab410*aa
            f(ii,id,1,ir410)   = f(ii,id,1,ir410)   - at*xn(ii,id,1)
            f(ii,id,2,ir410)   = f(ii,id,2,ir410)   - at*xn(ii,id,2)
            f(ii,id,3,ir410)   = f(ii,id,3,ir410)   - at*xn(ii,id,3)
            f(ii,id,1,ir410+1) = f(ii,id,1,ir410+1) - ab*xn(ii,id,1)
            f(ii,id,2,ir410+1) = f(ii,id,2,ir410+1) - ab*xn(ii,id,2)
            f(ii,id,3,ir410+1) = f(ii,id,3,ir410+1) - ab*xn(ii,id,3)
         end do
      end do
 
      do id=1,nd
         do ii=1,(nt+1)**2
            aa = cr660*arne(ii)
     &         *(at660*(temp(ii,id,ir660)   - tmprf(ir660))
     &         + ab660*(temp(ii,id,ir660+1) - tmprf(ir660+1)))
            at = at660*aa
            ab = ab660*aa
            f(ii,id,1,ir660)   = f(ii,id,1,ir660)   - at*xn(ii,id,1)
            f(ii,id,2,ir660)   = f(ii,id,2,ir660)   - at*xn(ii,id,2)
            f(ii,id,3,ir660)   = f(ii,id,3,ir660)   - at*xn(ii,id,3)
            f(ii,id,1,ir660+1) = f(ii,id,1,ir660+1) - ab*xn(ii,id,1)
            f(ii,id,2,ir660+1) = f(ii,id,2,ir660+1) - ab*xn(ii,id,2)
            f(ii,id,3,ir660+1) = f(ii,id,3,ir660+1) - ab*xn(ii,id,3)
         end do
      end do

	end subroutine
      
      
*dk gradient
	subroutine gradient(g,s)
!	This routine computes the gradient g of the nodal scalar field s.
 
      include 'size.h'
      real g(0:nt,nt+1,nd,3,nr+1), s(0:nt,nt+1,nd,nr+1)
      real w(2,0:nt,nr+1)
      common /grad/ rg(3,2,nr+1), grd(7,0:nt,nt+1,3,2)
      common /clck/ itmng, sec(50)
      if(itmng.eq.1) call mytime(tin)
 
      do id=1,nd
 
         do i2=1,nt+1
 
            do j=1,3
 
               do ir=1,nr,2
                  do i1=0,nt
                     w(1,i1,ir) =
     &                ((((((grd(1,i1,i2,j,1)*s(i1  ,i2  ,id,ir)
     &                    + grd(2,i1,i2,j,1)*s(i1+1,i2  ,id,ir))
     &                    + grd(3,i1,i2,j,1)*s(i1  ,i2+1,id,ir))
     &                    + grd(4,i1,i2,j,1)*s(i1-1,i2+1,id,ir))
     &                    + grd(5,i1,i2,j,1)*s(i1-1,i2  ,id,ir))
     &                    + grd(6,i1,i2,j,1)*s(i1  ,i2-1,id,ir))
     &                    + grd(7,i1,i2,j,1)*s(i1+1,i2-1,id,ir))
                     w(1,i1,ir+1) =
     &                ((((((grd(1,i1,i2,j,1)*s(i1  ,i2  ,id,ir+1)
     &                    + grd(2,i1,i2,j,1)*s(i1+1,i2  ,id,ir+1))
     &                    + grd(3,i1,i2,j,1)*s(i1  ,i2+1,id,ir+1))
     &                    + grd(4,i1,i2,j,1)*s(i1-1,i2+1,id,ir+1))
     &                    + grd(5,i1,i2,j,1)*s(i1-1,i2  ,id,ir+1))
     &                    + grd(6,i1,i2,j,1)*s(i1  ,i2-1,id,ir+1))
     &                    + grd(7,i1,i2,j,1)*s(i1+1,i2-1,id,ir+1))
                     w(2,i1,ir) =
     &                ((((((grd(1,i1,i2,j,2)*s(i1  ,i2  ,id,ir)
     &                    + grd(2,i1,i2,j,2)*s(i1+1,i2  ,id,ir))
     &                    + grd(3,i1,i2,j,2)*s(i1  ,i2+1,id,ir))
     &                    + grd(4,i1,i2,j,2)*s(i1-1,i2+1,id,ir))
     &                    + grd(5,i1,i2,j,2)*s(i1-1,i2  ,id,ir))
     &                    + grd(6,i1,i2,j,2)*s(i1  ,i2-1,id,ir))
     &                    + grd(7,i1,i2,j,2)*s(i1+1,i2-1,id,ir))
                     w(2,i1,ir+1) =
     &                ((((((grd(1,i1,i2,j,2)*s(i1  ,i2  ,id,ir+1)
     &                    + grd(2,i1,i2,j,2)*s(i1+1,i2  ,id,ir+1))
     &                    + grd(3,i1,i2,j,2)*s(i1  ,i2+1,id,ir+1))
     &                    + grd(4,i1,i2,j,2)*s(i1-1,i2+1,id,ir+1))
     &                    + grd(5,i1,i2,j,2)*s(i1-1,i2  ,id,ir+1))
     &                    + grd(6,i1,i2,j,2)*s(i1  ,i2-1,id,ir+1))
     &                    + grd(7,i1,i2,j,2)*s(i1+1,i2-1,id,ir+1))
                  end do
               end do
 
               do i1=0,nt
                  w(1,i1,nr+1) =
     &             ((((((grd(1,i1,i2,j,1)*s(i1  ,i2  ,id,nr+1)
     &                 + grd(2,i1,i2,j,1)*s(i1+1,i2  ,id,nr+1))
     &                 + grd(3,i1,i2,j,1)*s(i1  ,i2+1,id,nr+1))
     &                 + grd(4,i1,i2,j,1)*s(i1-1,i2+1,id,nr+1))
     &                 + grd(5,i1,i2,j,1)*s(i1-1,i2  ,id,nr+1))
     &                 + grd(6,i1,i2,j,1)*s(i1  ,i2-1,id,nr+1))
     &                 + grd(7,i1,i2,j,1)*s(i1+1,i2-1,id,nr+1))
                  w(2,i1,nr+1) =
     &             ((((((grd(1,i1,i2,j,2)*s(i1  ,i2  ,id,nr+1)
     &                 + grd(2,i1,i2,j,2)*s(i1+1,i2  ,id,nr+1))
     &                 + grd(3,i1,i2,j,2)*s(i1  ,i2+1,id,nr+1))
     &                 + grd(4,i1,i2,j,2)*s(i1-1,i2+1,id,nr+1))
     &                 + grd(5,i1,i2,j,2)*s(i1-1,i2  ,id,nr+1))
     &                 + grd(6,i1,i2,j,2)*s(i1  ,i2-1,id,nr+1))
     &                 + grd(7,i1,i2,j,2)*s(i1+1,i2-1,id,nr+1))
               end do
 
               do i1=0,nt
                  g(i1,i2,id,j,1) =
     &             (((rg(2,1,1)*w(1,i1,1)
     &              + rg(2,2,1)*w(2,i1,1))
     &              + rg(3,1,1)*w(1,i1,2))
     &              + rg(3,2,1)*w(2,i1,2))
               end do
 
               do ir=2,nr
                  do i1=0,nt
                     g(i1,i2,id,j,ir) =
     &                (((((rg(1,1,ir)*w(1,i1,ir-1)
     &                   + rg(1,2,ir)*w(2,i1,ir-1))
     &                   + rg(2,1,ir)*w(1,i1,ir  ))
     &                   + rg(2,2,ir)*w(2,i1,ir  ))
     &                   + rg(3,1,ir)*w(1,i1,ir+1))
     &                   + rg(3,2,ir)*w(2,i1,ir+1))
                  end do
               end do
 
               do i1=0,nt
                  g(i1,i2,id,j,nr+1) =
     &             (((rg(1,1,nr+1)*w(1,i1,nr)
     &              + rg(1,2,nr+1)*w(2,i1,nr))
     &              + rg(2,1,nr+1)*w(1,i1,nr+1))
     &              + rg(2,2,nr+1)*w(2,i1,nr+1))
               enddo
            enddo
         enddo
      enddo
 
      call rotate(g,nd,nr,nt,-1)
      call comm3s(g,nr,nt,3)
 
      if(itmng==1) then
      	call mytime(tout)
      	sec(19) = sec(19) + tout - tin
      endif
      
	end subroutine

