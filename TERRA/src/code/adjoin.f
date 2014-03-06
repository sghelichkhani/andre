*dk convectadjoint
	subroutine convectadjoint(tstepadj,iadj,itforw)

	include 'size.h'
	include 'pcom.h'
	include 'para.h'
	
	common /flds/ f((nt+1)**2*nd*(nr+1),8)
	common /velo/ upb(nt+1), u((nt+1)**2,nd,3,nr+1),  upe(nt+1)
	common /temp/ tpb(nt+1), temp((nt+1)**2,nd,nr+1), tpe(nt+1)
	common /pres/ ppb(nt+1), pres((nt+1)**2,nd,nr+1), ppe(nt+1)
	common /mesh/ xn((nt+1)**2,nd,3)
	common /radl/ rshl(nr+1), ird
	common /ndar/ arn((nt+1)**2),  arne((nt+1)**2),
     &              rarn((nt+1)**2), rarne((nt+1)**2)
	common /work/ wk((nt+1)**2,(nr+1)*2)
	common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
	common /conv/ step0, stepmin, stepmax
	common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
	common /eos2/ rhocv(nr+1), gamma(nr+1), alpha(nr+1), prf(nr+1),
     &               cond(nr+1)
	common /mgrd/ itlimit, convtol, itsolve
	common /vis1/ vscmax, rvscscl, tvscscl, pwrlawn, pwrlawsr, yldstrs
	common /io01/ casenum, gpath,  lpath, opath, tpath
	common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output
	common /solv/ itreq, convrg, divnorm, cutdiv
	common /call/ ncall
	common /prty/ propr(20)
	common /name/ titl(4,4)
	common /fnam/ vname, sname, cname
	common /clck/ itmng, sec(50)

	integer iter, iadj, itforw, ncall, ibc0
	real tb01, tb02, a, tstepadj(itmax+1)
	character*8 titl, vname, sname, cname
	character gpath*80, lpath*80, opath*80, tpath*80
	
	call mytime(tin)
	call nulvec(sec,50)

	ncall = 0
	a = randomnum(-1)
	time=tstepadj(itforw+1)

!	The adjoint velocity boundary conditions are no-slip on the
!	outer surface and free-slip on the inner boundary for the adjoint
!	backward integration. The adjoint temperature boundary conditions
!	are always zero on both boundaries. While we store the BC for the 
!	forward run, here we define new "adjoint" BC for the adjoint run.
	ibc0=ibc
	if(ibc==6) ibc0=5	
	tb01  = tb(1)
	tb02  = tb(2)
	tb(1) = 0.0
	tb(2) = 0.0

	call opinit
	call eosset

!	Peter: Use the residual temperature field at the end of the forward
!	run as the initial condition of the adjoint backward integration.
	call diffin

!	Begin main time step loop.
!	We need to define the number of backward timesteps based on the 
!	actual number of timesteps we took in the forward run. This info
!	from the forward run is stored in the 'itforw' variable.
	do iter=itforw,1,-1

		tstep=tstepadj(iter)
		if(mynum==0) then
			write(6,*) iter, time/1000000.0,
     &			3.1688e-8*tstep/1000000.0
   			write(8,*) iter, time/1000000.0,
     &			3.1688e-8*tstep/1000000.0
     		endif

		call advanceadjoint(tstep,ibc0,iter,itforw)

!	Peter: Reverse time calculation for backward integration.
		time  = time  - 3.1688e-8*tstep
		propr(2) = time			
	enddo

!	Peter: Write the adjoint temperature, i.e. the gradient of
!	the cost function at the end of the adjoint integration
!	and update the initial temperature field
	call adjointout(iadj)
	call perturbinitadjoint(iadj,iter)
	
	call mytime(tout)
	sec(3) = sec(3) + tout - tin

!	Peter: Set velocity and temperature boundary conditions back to 
!	their original value for use in the next iteration forward run.
	tb(1)= tb01
	tb(2)= tb02
	
	end subroutine


*dk advanceadjoint
	subroutine advanceadjoint(tstep,ibc0,iter,itforw)
!	This routine advances the adjoint Navier-Stokes 
!	equations one time step back in time.

	include 'size.h'
	include 'pcom.h'
	include 'para.h'

	common /flds/ f((nt+1)**2*nd*(nr+1),8)
	common /temp/ tpb(nt+1), temp((nt+1)**2*nd*(nr+1)), tpe(nt+1)
	common /shht/ shb(nt+1),  shr((nt+1)**2*nd*(nr+1)), she(nt+1)
	common /tdot/ tdb(nt+1), tdot((nt+1)**2*nd*(nr+1),2), tde(nt+1)
	common /velo/ upb(nt+1), u((nt+1)**2*nd*3*(nr+1)),  upe(nt+1)
	common /conv/ step0, stepmin, stepmax
	common /nrms/ fnrm, rnrm, unrm, ekin
	common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
	common /heat/ htop, hbot, hrad, heat, hshr, hvol, hnet, tnrm,
     &              tav(nr+1), qc(nr)
	common /eos2/ rhocv(nr+1), gamma(nr+1), alpha(nr+1), prf(nr+1),
     &               cond(nr+1)
	common /io01/ casenum, gpath,  lpath, opath, tpath
	common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output
	common /radl/ rshl(nr+1), ird
	common /solv/ itreq, convrg, divnorm, cutdiv
	common /mgrd/ itlimit, convtol, lreport, itsolve, izerou
	common /call/ ncall
	common /clck/ itmng, sec(50)

	integer iter, ibc0, itforw, iplot
	real fu(3*nv), ftemp(nv), tstep
	
	character opath*80

	if(itmng==1) call mytime(tin)

	msgnum = 1
	ncall=ncall+1

!	Peter: Read forward temperature and velocity field.
	call forwardin(fu,ftemp,iter,itforw)

!	Update the viscosity variation field.
	call oprupdate(ibc0)

!	Update the shear heating rate if appropriate.
	if(ieos>=10.and.mod(ncall,4)==1) call shrheat(shr,u,f,f(1,3))

	if(iter==itforw-3) then
		iplot=1
	else
		iplot=0
	endif
!	Use a second-order Runge-Kutta scheme to update the temperature.
!	Peter: Use modified momentum and energy routines for backward run.
	if((mod(ncall,time_skp)==0.or.ncall==1).and.adj_u==1) then
		call usolveadjoint(ftemp,ibc0,0)
	endif
	call energyadjoint(fu,tdot(1,1),0)

!	Peter: Reverse the time-stepping for the backward integration.
	do ii=1,nv
		temp(ii) = temp(ii) - 0.5*tstep*tdot(ii,1)
	enddo

	if((mod(ncall,time_skp)==0.or.ncall==1).and.adj_u==1) then
		call usolveadjoint(ftemp,ibc0,iplot)
	endif
	call energyadjoint(fu,tdot(1,2),iplot)

	do ii=1,nv
		temp(ii) = temp(ii) - tstep*(tdot(ii,2) - 0.5*tdot(ii,1))
	enddo

!	Compute norms and check energy conservation.
	if(mod(ncall,10)==0) then
		call norm3s(u,unrm,3,nd,nr,nt)
		call kenergy(ekin,u)
		call layrav(temp,tav)
		call heating(heat,hrad,tdot(1,2))
		tnrm = smean(temp)
		hnet = htop - hbot + heat - hrad
		do ir=1,nr
			qc(ir) = 12.56637*(tav(ir+1) - tav(ir))*rshl(ir)*rshl(ir+1)
     &			*sqrt(cond(ir)*cond(ir+1))/(rshl(ir) - rshl(ir+1))
		enddo
	endif

	if(ncall==1) call norm3s(u,unrm,3,nd,nr,nt)
      
	if(itmng==1) then
		call mytime(tout)
		sec(4) = sec(4) + tout - tin
	endif
	     
	end subroutine


*dk usolveadjoint
	subroutine usolveadjoint(ftemp,ibc0,inorm)
!	This routine solves for the new velocity field u.

	include 'size.h'
	include 'pcom.h'
	include 'para.h'
      
	common /flds/ f((nt+1)**2*nd,3,nr+1), r((nt+1)**2*nd,nr+1),
     &              v((nt+1)**2*nd,3,nr+1), s((nt+1)**2*nd,nr+1)
	common /mgwk/ w((nt+1)**2*nd,nr+1,5)
	common /temp/ tpb(nt+1), temp((nt+1)**2*nd*(nr+1)), tpe(nt+1)
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
	common /clck/ itmng, sec(50)

	integer ibc0, map(0:nt,nt+1,nd), inorm
	integer itr, izerou
	real ftemp(nv), urot(3,pl_size)
	real divnorm, cutdiv

	map=0
	urot=0.0
	if(itmng==1) call mytime(tin)

!	Compute the nodal force field f.
	call forcesadjoint(ftemp,inorm)

!	Solve for an estimate of the new velocity field.
	wd = 1.0
	wm	= 0.0
 
	if(vscmax>0.0) call uscale(u,1,nr)

	! calculates u from f (Au=f)
	! input u here serves as a first guess
	call multigrid(u,f,3,map,urot,ibc0,0,0)

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

	! pressure correction
	do while(itr<=npres_adj.and.divnorm>=cutdiv)

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

			call multigrid(v,f,3,map,urot,ibc0,izerou,0)
			izerou=0

			if(vscmax>0.0) call uscale(v,2,nr)			

			call divergence(w(1,1,2),w,v,w(1,1,2),w(1,1,5),w,1)
			call dotprod(a0,w(1,1,2),s)
 
			if(a0/=0.0) a0 = -divnorm**2/a0
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

	unrm = velnrm
	fnrm = fnrm0
	rnrm = rnrm0

	if(itmng==1) then
		call mytime(tout)
		sec(5) = sec(5) + tout - tin
	endif
	
	end subroutine
            

*dk forcesadjoint
	subroutine forcesadjoint(ftemp,inorm)

	include 'size.h'
	include 'pcom.h'
	
	common /flds/ f((nt+1)**2*nd,3,nr+1), b((nt+1)**2*nd,nr+1),
     &              w((nt+1)**2*nd,4,nr+1)
	common /velo/ upb(nt+1), u((nt+1)**2*nd,3,nr+1),  upe(nt+1)
	common /pres/ ppb(nt+1), pres((nt+1)**2*nd,nr+1), ppe(nt+1)
	common /temp/ tpb(nt+1), temp((nt+1)**2*nd,(nr+1)), tpe(nt+1)
	common /mesh/ xn((nt+1)**2*nd,3)
	common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
	common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
	common /vis1/ vscmax, rvscscl, tvscscl, pwrlawn, pwrlawsr, yldstrs
	common /clck/ itmng, sec(50)
	common /radl/ rshl(nr+1), ird	

	integer ir_beg(3), ir_end(3), inorm
	real lbeg, lend, tgrad((nt+1)**2*nd,3,nr+1), sum_s(3)
	real ftemp(nv)

	if(itmng==1) call mytime(tin)

	ir_beg=1
	ir_end=1
	do i=1,3
	if(i==1) then
		lbeg=0.0
		lend=300.0
	elseif(i==2) then
		lbeg=350.0
		lend=2400.0
	else
		lbeg=2500.0
		lend=3000.0
	endif

	do ir=1,nr+1
		if(lbeg-(6370.0-rshl(ir)/1000.0)>=0) ir_beg(i)=ir
		if(lend-(6370.0-rshl(ir)/1000.0)>=0) ir_end(i)=ir
	enddo
	enddo
	
!	Compute gradient 'w' of the forward temperature field 'ftemp'
!	This will be the 'adjoint-buoyancy' term.
	call gradient(w,ftemp)

!	Compute the gradient of the 'adjoint-pressure'.
	call gradient(f,pres)

	do ir=1,nr+1
		aa = rhorf(ir)/rho0
!	Peter: Add the adjoint momentum-forcing term (which comes from the
!	adjoint temperature times the gradient of forward temperature) to
!	the gradient of the 'adjoint-pressure'. (Previously this loop 
!	added the buoyancy term to the pressure gradient)
		do ii=1,(nt+1)**2*nd
			tgrad(ii,1,ir) = w(ii,1,ir)*temp(ii,ir)
			tgrad(ii,2,ir) = w(ii,2,ir)*temp(ii,ir)
			tgrad(ii,3,ir) = w(ii,3,ir)*temp(ii,ir)
			f(ii,1,ir) = aa*f(ii,1,ir) + tgrad(ii,1,ir)
			f(ii,2,ir) = aa*f(ii,2,ir) + tgrad(ii,2,ir)			
			f(ii,3,ir) = aa*f(ii,3,ir) + tgrad(ii,3,ir)
		enddo
	enddo
	
	if(inorm==1) then
	do i=1,3
		 call sum_sqrt(tgrad,sum_s(i),ir_beg(i),ir_end(i))
		 sum_s(i)=sqrt(sum_s(i)/visc)
	enddo
	if(mynum==0) then
		write(223,*) sum_s(1), sum_s(2), sum_s(3)
	endif
	call gmtout2(tgrad/visc,3,0,0,0,'g')
	endif

	if(vscmax>0.0) call uscale(f,2,nr)
	aa = 1.0/visc
	do ir=1,nr+1
		do ii=1,(nt+1)**2*nd
			f(ii,1,ir) = aa*f(ii,1,ir)
			f(ii,2,ir) = aa*f(ii,2,ir)
			f(ii,3,ir) = aa*f(ii,3,ir)
		enddo
	enddo

	if(itmng==1) then
		call mytime(tout)
		sec(15) = sec(15) + tout - tin
	endif

	end subroutine
      
      
*dk energyadjoint
	subroutine energyadjoint(fu,dhdt,inorm)
!	This routine updates the temperature field.

	include 'size.h'
	include 'pcom.h'
	include 'para.h'

	common /flds/ f((nt+1)**2,nd,nr+1,8)
	common /temp/ tpb(nt+1), temp((nt+1)**2,nd,nr+1), tpe(nt+1)
	common /velo/ upb(nt+1), u((nt+1)**2,nd,3,nr+1),  upe(nt+1)
	common /shht/ shb(nt+1),  shr((nt+1)**2,nd,nr+1), she(nt+1)
	common /volm/ vol((nt+1)**2,nr+1,2)
	common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
	common /heat/ htop, hbot, hrad, heat, hshr, hvol, hnet, tnrm,
     &              tav(nr+1), qc(nr)
	common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
	common /eos2/ rhocv(nr+1), gamma(nr+1), alpha(nr+1), prf(nr+1),
     &               cond(nr+1)
	common /conv/ step0, stepmin, stepmax
	common /mesh/ xn((nt+1)**2,nd,3)
	common /radl/ rshl(nr+1), ird
	common /volc/ tmpmlt
	common /call/ ncall
	common /clck/ itmng, sec(50)

	integer iter, inorm, ir_beg(3), ir_end(3)
	real dhdt((nt+1)**2,nd,nr+1), dhdt_adju((nt+1)**2,nd,nr+1)
	real fu(3*nv), norm1, norm2, norm3
	real lbeg, lend, sum_s(3)
	
	if(itmng==1) call mytime(tin)

	if(ncall==1) call facegen(f)

	ir_beg=1
	ir_end=1
	do i=1,3
	if(i==1) then
		lbeg=0.0
		lend=300.0
	elseif(i==2) then
		lbeg=350.0
		lend=2400.0
	else
		lbeg=2500.0
		lend=3000.0
	endif

	do ir=1,nr+1
		if(lbeg-(6370.0-rshl(ir)/1000.0)>=0) ir_beg(i)=ir
		if(lend-(6370.0-rshl(ir)/1000.0)>=0) ir_end(i)=ir
	enddo
	enddo

!	Compute the rate of temperature change due to advection and
!	conduction.
	call advectadjoint(fu,dhdt)	
	if(inorm==1) then
		call norm3s(dhdt,norm1,1,nd,nr,nt)
		call gmtout2(dhdt,1,0,0,0,'h')
	endif
	
!	Peter: Include the adjoint-energy source term, which consists of the
!	product of 'adjoint-velocities' and a buoyancy term without temperature.
!	In other words, this term is derived from the forward-buoyancy term.
	if(adj_u==1) then
	do ir=1,nr+1
		do id=1,nd
			do ii=1,(nt+1)**2 
				dhdt_adju(ii,id,ir) = alpha(ir)*grv(ir)*rhorf(ir)
     &			* (u(ii,id,1,ir)*xn(ii,id,1) + u(ii,id,2,ir)*xn(ii,id,2)
     &				+ u(ii,id,3,ir)*xn(ii,id,3))
     				dhdt(ii,id,ir)=dhdt(ii,id,ir)+dhdt_adju(ii,id,ir)
			enddo
		enddo
	enddo
	endif
	
	if(inorm==1) then
	do i=1,3
		 call sum_sqrt(dhdt_adju,sum_s(i),ir_beg(i),ir_end(i))
		 sum_s(i)=sqrt(sum_s(i))
	enddo
	call norm3s(dhdt_adju,norm2,1,nd,nr,nt)
	if(mynum==0) then
		write(223,*) norm1, norm2
		write(223,*) sum_s(1), sum_s(2), sum_s(3)
	endif
	call gmtout2(dhdt_adju,1,0,0,0,'d')
	endif
	
	htop = 0.0
	hbot = 0.0

!	Use heating rate of the nodal boundary layers as proxy for
!	the heat flux leaving or entering the shell.
	do id=1,nd

		do ii=1,(nt+1)**2
			htop = htop + dhdt(ii,id,1)*vol(ii,1,1)
			dhdt(ii,id,1) = 0.
		enddo

!	Peter: Always enforce adjoint temperature boundary conditions (rather
!	than adjoint temperature flux boundary conditions) on the adjoint
!	temperature both at the top AND! the bottom (at the CMB). This is,
!	because we will assume that the forward temperature is fixed both
!	at the top and the bottom (CMB).
		do ii=1,(nt+1)**2
			hbot = hbot - dhdt(ii,id,nr+1)*vol(ii,nr+1,1)
			dhdt(ii,id,nr+1) = 0.
		enddo

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

	if(itmng==1) then
		call mytime(tout)
		sec(20) = sec(20) + tout - tin
	endif

	end subroutine


*dk advectadjoint
	subroutine advectadjoint(fu,dhdt)
!	This routine computes the heat advected and conducted per unit
!	time and also accounts for phase transitions if they are present.
!	Note that "fu" corresponds to velocities of the forward run.

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
	real fu((nt+1)**2,nd,3,nr+1)

	if(itmng==1) call mytime(tin)

	call nulvec(dhdt,(nt+1)**2*nd*(nr+1))

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
			!	Peter: Use original velocities from forward run (called fu).
				uf1  = 0.5*(fu(ii,id,1,ir) + fu(ii,id,1,jr))
				uf2  = 0.5*(fu(ii,id,2,ir) + fu(ii,id,2,jr))
				uf3  = 0.5*(fu(ii,id,3,ir) + fu(ii,id,3,jr))
				vfl  = ((uf1*xn(ii,id,1) + uf2*xn(ii,id,2))
     &                + uf3*xn(ii,id,3))*r2*arn(ii)
			!	Peter: Advect with true velocities from forward run.
				tadv = 0.5*(temp(ii,id,ir) + temp(ii,id,jr))*vfl
			!	Peter: Reverse the sign of the adjoint heat conduction term.
				cnd  = ((temp(ii,id,jr)-temp(ii,id,ir))*cr*arn(ii)) *-1.0
				dh(ii,1) = vfl
				dhdt(ii,id,ir) = (dhdt(ii,id,ir) + tadv*rhocv(ir)) + cnd
     &                        +  temp(ii,id,ir)*vfl*ai
				dhdt(ii,id,jr) = (dhdt(ii,id,jr) - tadv*rhocv(jr)) - cnd
     &                        -  temp(ii,id,jr)*vfl*aj
			enddo

!	Add contributions from lateral advection and conduction.
			do ifc=1,3

				jj = i1n(ifc) + (nt+1)*i2n(ifc)
				ib = max(1, 1-jj)
				ie = min((nt+1)**2, (nt+1)**2-jj)

				do ii=ib,ie
				!	Peter: Use original velocities from forward run (called fu).
					uf1  = 0.5*(fu(ii,id,1,jr) + fu(ii+jj,id,1,jr))
					uf2  = 0.5*(fu(ii,id,2,jr) + fu(ii+jj,id,2,jr))
					uf3  = 0.5*(fu(ii,id,3,jr) + fu(ii+jj,id,3,jr))
					vfl  = ((uf1*facet(ii,id,ifc,1)
     &                   + uf2*facet(ii,id,ifc,2))
     &                   + uf3*facet(ii,id,ifc,3))*psi(ii,ifc)*ra
				!	Peter: Advect with true velocities from forward run.
					hadv = 0.5*(temp(ii,id,jr) + temp(ii+jj,id,jr))
     &                      *vfl*rhocv(jr)
				!	Peter: Reverse the sign of the adjoint heat conduction term.
					cnd  = (ct*(temp(ii,id,jr) - temp(ii+jj,id,jr))
     &                  *psi(ii,ifc)*rdxn(ii,ifc)) * -1.0
					dh(ii,1) = hadv + cnd
					dh(ii,2) = vfl*aj
					dhdt(ii,id,jr) = dhdt(ii,id,jr) - dh(ii,1)
     &                           - temp(ii,id,jr) * dh(ii,2)
				enddo

				do ii=ib,ie
					dhdt(ii+jj,id,jr) = dhdt(ii+jj,id,jr) + dh(ii,1)
     &                              + temp(ii+jj,id,jr) * dh(ii,2)
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


*dk perturbinitadjoint
	subroutine perturbinitadjoint(iadj,iter)
!	This routine updates the initial condition with a correction field
!	obtained from the adjoint-backward-integration and writes the
!	new field back to tape

	include 'size.h'
	include 'pcom.h'
	include 'para.h'	
	
	common /flds/ f((nt+1)**2*nd*(nr+1),8)
	common /temp/ tpb(nt+1), temp((nt+1)**2*nd*(nr+1)), tpe(nt+1)
	common /velo/ upb(nt+1), u((nt+1)**2,nd,3,nr+1),  upe(nt+1)
	common /pres/ ppb(nt+1), pres((nt+1)**2,nd,nr+1), ppe(nt+1)
	common /conv/ step0, stepmin, stepmax
	common /heat/ htop, hbot, hrad, heat, hshr, hvol, hnet, tnrm,
     &              tav(nr+1), qc(nr)
	common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
	common /radl/ rshl(nr+1), ird
	common /io01/ casenum, gpath,  lpath, opath, tpath
	common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output

	integer irlimit, iter, iadj, ii, ik, lay, layskp_appl, lay_damp
	real orshl(nr+1), diffmax, dfac
	real temp_tomo((nt+1)**2*nd,nr+1)
	real min_tomo(nr+1), max_tomo(nr+1), mean_tomo(nr+1)
	real min_glob(nr+1), max_glob(nr+1), mean_glob(nr+1)
	real tol_fac, max_border, min_border
	
	character char1*4, cname*15, casenum*3, char2*2
	character*8 titl, otitl(4,4)
	character gpath*80, lpath*80, opath*80

	write(char1,'(I4.4)') mynum
	write(char2,'(I2.2)') iadj

	cname= 'b'//casenum//'.'//char1//'.00_'//char2
	open(20, file=trim(opath)//'c-files/'//cname,status='unknown')
	call vecin(f(1,1),otitl,orshl,1,nr,20,1)
	close(20)

!	read in old inital temperature field
	cname=	'c'//casenum//'.'//char1//'.00'
	open(61, file=trim(gpath)//cname, status='unknown')
	call vecin(temp,otitl,orshl,1,nr,61,1)
	close(61)

!	Add a fraction of the gradient of the cost-function. The
!	damping factor is set to a rather restrictive value of 0.8
	dfac=0.8
!	Skip the first 'res_skp_appl' km.
	do ir=1,nr+1
		if((rshl(1)-rshl(ir))/1000.0<=res_skp_appl) layskp_appl=ir
	enddo
	layskp_apply=min(layskp_appl,nr)
	do ii=(nt+1)**2*nd*layskp_appl+1,(nt+1)**2*nd*(nr+1)
		temp(ii) = temp(ii) + dfac*f(ii,1)
	enddo

! ################################################
! ########### DAMPING *NEW* *NEW* ####################
! ################################################
!	Read in the reference model to compute the min/max temperatures
	cname=	'c511.'//char1//'.01'
	open(62, file=trim(gpath)//cname, status='unknown')
	call vecin(temp_tomo,otitl,orshl,1,nr,62,1)
	close(62)
	
	do ii=1,nr+1
		max_tomo(ii)=0.0
		min_tomo(ii)=5000.0
		mean_tomo(ii)=0.0
		do ik=1,(nt+1)**2*nd
			if(temp_tomo(ik,ii)>max_tomo(ii)) then
				max_tomo(ii)=temp_tomo(ik,ii)
			endif
			if(temp_tomo(ik,ii)<min_tomo(ii)) then
				min_tomo(ii)=temp_tomo(ik,ii)
			endif
			mean_tomo(ii)=mean_tomo(ii)+temp_tomo(ik,ii)
		enddo
		mean_tomo(ii)=mean_tomo(ii)/real((nt+1)**2*nd)/real(nproc)
	enddo
	
!	Computation of the min/max temp values in each layer of the reference model (tomography)
	call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	call MPI_REDUCE(max_tomo,max_glob,nr+1,MPI_REAL8,MPI_MAX,0,
     &	MPI_COMM_WORLD,IERROR)
	call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	call MPI_BCAST(max_glob,nr+1,MPI_REAL8,0,MPI_COMM_WORLD,IERROR)
	
	call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	call MPI_REDUCE(min_tomo,min_glob,nr+1,MPI_REAL8,MPI_MIN,0,
     &	MPI_COMM_WORLD,IERROR)
	call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	call MPI_BCAST(min_glob,nr+1,MPI_REAL8,0,MPI_COMM_WORLD,IERROR)
	call MPI_BARRIER(MPI_COMM_WORLD,IERROR)

	call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	call MPI_REDUCE(mean_tomo,mean_glob,nr+1,MPI_REAL8,MPI_SUM,0,
     &	MPI_COMM_WORLD,IERROR)
	call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	call MPI_BCAST(mean_glob,nr+1,MPI_REAL8,0,MPI_COMM_WORLD,IERROR)
	call MPI_BARRIER(MPI_COMM_WORLD,IERROR)

!	Restrict temp. values in the initial condition in the first dpth_damp km to be
!	max. 500 K greater (plumes) or 1000 K smaller (slabs) than the reference mean value
	do ir=1,nr+1
		max_glob(ir)=max_glob(ir)/mean_glob(ir)
		min_glob(ir)=min_glob(ir)/mean_glob(ir)
	enddo

	! average temperature in each layer of the initial model
	call layrav(temp,tav)
	
	tol_fac=0.1
	! update temperatures to dpth_damp km.
	do ir=1,nr+1
		if((rshl(1)-rshl(ir))/1000.0<=dpth_damp) lay_damp=ir
	enddo
	do ii=1,(nt+1)**2*nd*lay_damp
		lay=ceiling(real(ii)/real((nt+1)**2*nd))
		
		!max_border=(1.0+tol_fac)*max_glob(lay)
		!min_border=(1.0-tol_fac)*min_glob(lay)
		
		!max_border=mean_glob(lay)+500.0
		!min_border=mean_glob(lay)-1000.0
		
		!if(temp(ii)>max_border) temp(ii) = max_border
		!if(temp(ii)<min_border) temp(ii) = min_border
		
		
		! ########### 11.10.11 #################
		! NEW DAMPING ! 
		! T deviations in initial state must not be larger than those in final state
		! ##################################
		if(lay>1) then
		if(temp(ii)/tav(lay)>max_glob(lay)) then
			!if(mynum==0) write(6,*) temp(ii), tav(lay), max_glob(lay)
			temp(ii)=tav(lay)*max_glob(lay)
			!if(mynum==0) then
			!	write(6,*) temp(ii)
			!	write(6,*)
			!endif
		endif
		if(temp(ii)/tav(lay)<min_glob(lay)) then
			!if(mynum==0) write(6,*) temp(ii), tav(lay), min_glob(lay)
			temp(ii)=tav(lay)*min_glob(lay)
			!if(mynum==0) then
			!	write(6,*) temp(ii)
			!	write(6,*)
			!endif
		endif
		endif
	enddo
! #################################################

!	Compute difference between 'true' and 'estimated' initial condition.
	cname= 'c511.'//char1//'.00'
	open (90, file=trim(gpath)//cname, status='unknown')
	call vecin(f(1,1),otitl,orshl,1,nr,90,1)
	call vecdiff(f(1,2),f(1,1),temp,((nt+1)**2)*nd*(nr+1))
	close(90)
	
	call proprty(iter)
	cname= 'd'//casenum//'.'//char1//'.00'
	open (91,file=trim(gpath)//cname, status='unknown')
	call vecout(f(1,2),rshl,1,nr,nt,nr,91,1)
	close(91)

!	Peter: Check norm and maxium value of the temperature residual (misfit)
!	between estimated and 'true' temp at the beginning of the forward run.
	call norm3s(f(1,2),diffnrm,1,nd,nr,nt)
	call normnew3s(f(1,2),diffnrm2,1)

	diffmax=0.0
	do ii=1,(nt+1)**2*nd*(nr+1)
		diffmax=max(diffmax,abs(f(ii,2)))
	enddo
	if(nproc>1) call pmax(diffmax)
      
	if(mynum==0) write(7,15) iadj, diffnrm, diffnrm2, diffmax
 15	format(' Difference in the initial state of the forward run after'/
     &' ADJT-ITER=', i4,' DIFFNRM=',1pe10.3,' DIFFNRM2=',1pe10.3,
     &' DIFFMAX=',1pe10.3)

!	Overwrite the temp/velocity initial condition field (first c-file)
	call fldsout(0,iter)
	call fldsout2(0,iadj+1)

	end subroutine
	
	
*dk normnew3s
	subroutine normnew3s(r,rnorm,nj)
	implicit none
!	This routine computes the l2 norm of r.

	include 'size.h'
	include 'pcom.h'

	common /clck/ itmng, sec

	integer nj, itmng, ic,id,i1,i2

	real sum1, rnorm, sec(50), tin, tout
	real r(0:nt,nt+1,nd,nj*(nr+1))

	if(itmng==1) call mytime(tin)
 
	if(mynum<mproc*10/nd) then
 
		sum1=0.0

		if(mynum==0.or.(mynum==mproc.and.nd<=5)) then
			do ic=4,nj*(nr+1)
				sum1 = sum1 + r(0,1,1,ic)**2
			enddo
		endif

		if(mynum==0.and.nd==10) then
			do ic=4,nj*(nr+1)
				sum1 = sum1 + r(0,1,6,ic)**2
			enddo
		endif

		do ic=4,nj*(nr+1)
			do id=1,nd
				do i2=1,nt
					do i1=1,nt
						sum1 = sum1 + r(i1,i2,id,ic)**2
					enddo
				enddo
			enddo
		enddo
 
		if(mproc*10/nd>1) call psum(sum1, 1)
 
		rnorm = sqrt(sum1/((nt*nt*10*mproc+2)*(nr+1-4)))
 
		if(itmng==1) then
			call mytime(tout)
			sec(13) = sec(13) + tout - tin
		endif

	endif

	end subroutine
	
	
*dk vecdiff
	subroutine vecdiff(s1,s2,s3,nn)
	implicit none
	
	integer nn, ii
	real s1(nn), s2(nn), s3(nn)

	do ii=1,nn
		s1(ii) = s2(ii)-s3(ii)
	enddo

	end subroutine

	
!	##########################################
!	######## NEW I/O SUBROUTINES #############
!	##########################################
*dk forwardout
	subroutine forwardout(idump)
!	This routine writes out the temperature and velocity
!	fields of the forward run at each time-step for later
!	use in the adjoint backward integration.

	include 'size.h'
	include 'pcom.h'

	common /temp/ tpb(nt+1), temp((nt+1)**2,nd,nr+1), tpe(nt+1)
	common /velo/ upb(nt+1), u((nt+1)**2,nd,3,nr+1),  upe(nt+1)
	common /radl/ rshl(nr+1), ird
	common /io01/ casenum, gpath,  lpath, opath, tpath
	common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output
      
	character char1*4, char2*4, cname*14, casenum*3
	character gpath*80, lpath*80, opath*80

	write(char1,'(I4.4)') mynum
	write(char2,'(I4.4)') idump
	
	cname='a'//casenum//'.'//char1//'.'//char2

	open(111, file=trim(gpath)//cname, form='unformatted',
     &status='unknown')
	call proprty(idump)
	call vecoutunform2(temp,1,111)
	call vecoutunform2(u,3,111)
	close(111)

	end subroutine  
	
	
*dk forwardout2
	subroutine forwardout2(idump)
!	This routine writes out the temperature and velocity
!	fields of the forward run at each time-step for later
!	use in the adjoint backward integration.

	include 'size.h'
	include 'pcom.h'

	common /temp/ tpb(nt+1), temp((nt+1)**2,nd,nr+1), tpe(nt+1)
	common /velo/ upb(nt+1), u((nt+1)**2,nd,3,nr+1),  upe(nt+1)
	common /radl/ rshl(nr+1), ird
	common /io01/ casenum, gpath,  lpath, opath, tpath
	common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output
      
	character char1*4, char2*4, cname*14, casenum*3
	character gpath*80, lpath*80, opath*80

	write(char1,'(I4.4)') mynum
	write(char2,'(I4.4)') idump
	
	cname='a'//casenum//'.'//char1//'.'//char2

	open(111, file=trim(opath)//'c-files/'//cname,
     &		 form='unformatted', status='unknown')
     
	call proprty(idump)
	call vecoutunform2(temp,1,111)
	call vecoutunform2(u,3,111)
	close(111)

	end subroutine  


*dk forwardin
	subroutine forwardin(fu,ftemp,idump,itforw)
!	This routine reads in the forward fields of temperature 
!	and velocity for use in the adjoint backward integration.

	include 'size.h'
	include 'para.h'
	include 'pcom.h'

	common /io01/ casenum, gpath,  lpath, opath, tpath
	common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output

	integer idump, nlpath, itforw, idump_tmp

	real fu(3*nv), ftemp(nv)

	character char1*4, char2*4, cname*14, casenum*3
	character gpath*80, lpath*80, opath*80

	write(char1,'(I4.4)') mynum
	
	idump_tmp=idump
	do while(mod(idump_tmp,adj_skp)>0)
		idump_tmp=idump_tmp+1
	enddo
	
	if(idump_tmp>itforw) idump_tmp=itforw

	write(char2,'(I4.4)') idump_tmp		
	cname='a'//casenum//'.'//char1//'.'//char2

	open(37,file=trim(gpath)//cname, form='unformatted',
     & status='unknown')
	
	call vecinunform2(ftemp,1,37)
	call vecinunform2(fu,3,37)
	close(37)

	end subroutine
	
      
*dk diffout
	subroutine diffout(iadj,iter)
!	This routine writes the misfit between the 'observed/true'
!	temperatures and the temperatures computed at the end of
!	the forward model run.

	include 'size.h'
	include 'pcom.h'
	include 'para.h'

	common /flds/ f((nt+1)**2*nd*(nr+1),8)
	common /temp/ tpb(nt+1), temp((nt+1)**2,nd,nr+1), tpe(nt+1)
	common /velo/ upb(nt+1), u((nt+1)**2,nd,3,nr+1),  upe(nt+1)
	common /radl/ rshl(nr+1), ird
	common /io01/ casenum, gpath,  lpath, opath, tpath
	common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output

	integer iadj, iter, layskp_calc
	real orshl(nr+1), diffnrm
	
	character*8 otitl(4,4)
	character char1*4, cname1*12, cname2*12      
	character gpath*80,lpath*80, casenum*3, opath*80

	write(char1,'(I4.4)') mynum
	cname1='c511.'//char1//'.01'
	cname2='d'//casenum//'.'//char1//'.01'

!	Peter: Read in the 'true' (observed) field at the final-stage.
	open (99, file=trim(gpath)//cname1, status='unknown')
	call vecin(f(:,1),otitl,orshl,1,nr,99,1)
	close(99)

!	Peter: Compute difference w.r.t. to 'estimated' final-stage.
!	(f2=f1-temp=TRUE-CALC)
	call vecdiff(f(:,2),f(:,1),temp,((nt+1)**2)*nd*(nr+1))

!	Peter: Check norm and maximum value of the temperature residual (misfit)
!	between estimated and 'true' temperatures at the end of forward run.
	call norm3s(f(:,2),diffnrm,1,nd,nr,nt)

	diffmax=0.0
	do ii=1,(nt+1)**2*nd*(nr+1)
		diffmax=max(diffmax,abs(f(ii,2)))
	enddo
	if(nproc>1) call pmax(diffmax)
	
	if(mynum==0) write(7,15) iadj, diffnrm, diffmax
 15   format(' Difference in the final state of the forward run'/
     &' ITERATION=', i4,'  DIFFNRM=',1pe10.3,'  DIFFMAX=',1pe10.3/)

!	Set the residual to zero for the first 'res_skp_calc' km.
	do ir=1,nr+1
		if((rshl(1)-rshl(ir))/1000.0<=res_skp_calc) layskp_calc=ir
	enddo
	do ii=1,(nt+1)**2*nd*layskp_calc
		f(ii,2)=0.0
	enddo

	call proprty(iter)
	open (79,file=trim(gpath)//cname2,status='unknown',
     & form='formatted')
	call vecout(f(:,2),rshl,1,nr,nt,nr,79,1)
	close(79)

	end subroutine


*dk diffin
	subroutine diffin
!	This routine reads in the temperature residual between
!	the "observed" temperature and the temperature that was
!	calculated at the end of the forward integration.

	include 'size.h'
	include 'pcom.h'

	common /temp/ tpb(nt+1), temp((nt+1)**2,nd,nr+1), tpe(nt+1)
	common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
	common /heat/ htop, hbot, hrad, heat, hshr, hvol, hnet, tnrm,
     &              tav(nr+1), qc(nr) 
	common /io01/ casenum, gpath,  lpath, opath, tpath
	common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output

	real orshl(nr+1)

	character*8 otitl(4,4)
	character char1*4, cname*12      
	character gpath*80,lpath*80, casenum*3, opath*80

	write(char1,'(I4.4)') mynum
	cname= 'd'//casenum//'.'//char1//'.01'

	open (77,file=trim(gpath)//cname, status='unknown')
	call vecin(temp,otitl,orshl,1,nr,77,1)
	close(77)

!	Impose adjoint temperature boundary conditions.
	do id=1,nd
		do ii=1,(nt+1)**2
			temp(ii,id,1)=tb(1)
			temp(ii,id,nr+1)=tb(2)
		enddo
	enddo

	call layrav(temp,tav)
	tnrm = smean(temp)

	end subroutine


*dk adjointout
	subroutine adjointout(iadj)
!	This routine writes out the adjoint solution at the end
!	of the adjoint-backward integration. A simple conjugate
!	gradient method is applied.

	include 'size.h'
	include 'pcom.h'

	common /temp/ tpb(nt+1), temp((nt+1)**2*nd,(nr+1)), tpe(nt+1)
	common /radl/ rshl(nr+1), ird
	common /io01/ casenum, gpath,  lpath, opath, tpath
	common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output

	integer iadj, istat
	real resid((nt+1)**2*nd,nr+1), orshl(nr+1)
	real delta, residnorm, tempnorm

	character*8 otitl(4,4)
	character char1*4, cname*15, opath*80, char2*2
	character gpath*80,lpath*80, casenum*3

	delta=0.0
	resid=0.0
	residnorm=0.0
	write(char1,'(I4.4)') mynum
	write(char2,'(I2.2)') iadj-1
	cname= 'b'//casenum//'.'//char1//'.00_'//char2

!	Apply a conjugate gradient method by adding a weighted
!	combination of previous search directions. The CG weight
!	is set to delta
	call norm3s(temp,tempnorm,1,nd,nr,nt)
	
	if(iadj>1) then
	open(99,file=trim(opath)//'c-files/'//cname,status='old',iostat=istat)
	if(istat==0) then
		! read residual of last step from file
		open(99,file=trim(opath)//'c-files/'//cname,status='unknown')
		call vecin(resid,otitl,orshl,1,nr,99,1)
		close(99)

		call norm3s(resid,residnorm,1,nd,nr,nt)
		if(residnorm>0) delta = (tempnorm/residnorm)**2
		do ir=1,(nr+1)
			do ii=1,(nt+1)**2*nd
				temp(ii,ir) = temp(ii,ir) + delta*resid(ii,ir)
			enddo
		enddo
	endif
	endif

	write(char2,'(I2.2)') iadj
	cname= 'b'//casenum//'.'//char1//'.00_'//char2
	call proprty(1)
	open(101,file=trim(opath)//'c-files/'//cname,status='unknown')
	call vecout(temp,rshl,1,nr,nt,nr,101,1)
	close(101)
	call gmtout2(temp,1,0,0,iadj,'b')

!	Write the step-length of the CG algorithm
	if(mynum==0) write(7,15) residnorm, tempnorm, delta
 15   format('  Conjugate Gradient Step Length in the adjoint-solution'/
     &'  RESID =',1pe10.3,'   TEMP =',1pe10.3,
     &'  DELTA = ',1pe10.3/)

	end subroutine
      
      
*dk vecinunform2
	subroutine vecinunform2(u,nj,nf)
	implicit none
!	This routine reads the nodal field u from logical unit nf
!	using unformated I/O and without any additional information.
 
	include 'size.h'
	include 'pcom.h'

	integer nf, nj, ii
	real u((nt+1)**2*nd*nj*(nr+1))

	read(nf) (u(ii),ii=1,(nt+1)**2*nd*nj*(nr+1))
 
	end subroutine


*dk vecoutunform2
	subroutine vecoutunform2(u,nj,nf)
	implicit none
!	This routine writes the nodal field u to logical unit nf
!	using unformated I/O and without any additional information.
 
	include 'size.h'
	include 'pcom.h'

	integer nf, nj, ii
	real u((nt+1)**2*nd*nj*(nr+1))
 
	write(nf) (u(ii),ii=1,(nt+1)**2*nd*nj*(nr+1))
 
	end subroutine

!##############################################

	subroutine sum_sqrt(s,sum_s,ir_beg,ir_end)
	implicit none

	include 'size.h'
	include 'pcom.h'

	integer k, ir, id, ii, ir_beg, ir_end
	real s((nt+1)**2*nd*(nr+1))
	real sum_s_pro, sum_s

	k=(ir_beg-1)*nd*(nt+1)**2+1
	sum_s=0.0
	sum_s_pro=0.0

	do ir=ir_beg,ir_end	
	do id=1,nd
		do ii=1,(nt+1)**2
			sum_s_pro=sum_s_pro+s(k)*s(k)
			k=k+1
		enddo
	enddo
	enddo
	
	sum_s_pro=sum_s_pro/((ir_end-ir_beg+1)*nproc*nd*(nt+1)**2)
	
	call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	call MPI_REDUCE(sum_s_pro,sum_s,1,MPI_REAL8,MPI_SUM,0,
     &	MPI_COMM_WORLD,IERROR)
	call MPI_BCAST(sum_s,1,MPI_REAL8,0,MPI_COMM_WORLD,IERROR)
	call MPI_BARRIER(MPI_COMM_WORLD,IERROR)

	end subroutine sum_sqrt


!##############################################

	subroutine gmtout2(vec,isize,imean,idump,iadj,char_file)
	implicit none
!	This routine writes out scalar fields and absolute values of
!	vector fields on a 'isamp' (see 'para.h') level coarser grid in long/lat/val format
 
	include 'size.h'
	include 'para.h'
	include 'pcom.h'

	common /mesh/ xn
	common /radl/ rshl
	common /io01/ casenum, gpath, lpath, opath, tpath

	integer isize, k, ir, id, ii, ipro, ifile, ihelp, istat, iadj, idump
	integer imean
	real xn((nt+1)**2,nd,3)
	real vec((nt+1)**2*nd*(nr+1)*isize), vec_fine((nt+1)**2*nd*(nr+1))
	real xn_coar((nt/2+1)**2,nd,3), xn_coar2((nt/4+1)**2,nd,3)
	real lat((nt/max(1,2*isamp)+1)**2,nd)
	real long((nt/max(1,2*isamp)+1)**2,nd)
	real vec_end((nt/max(1,2*isamp)+1)**2*nd*(nr/max(1,2*isamp)+1))
	real rsh(nr/max(1,2*isamp)+1), rshl(nr+1)
	real vec_coar((nt/2+1)**2*nd*(nr/2+1))
	real vec_coar2((nt/4+1)**2*nd*(nr/4+1))
	real vec_fine_av(nr+1), vec_av(nr/max(1,2*isamp)+1)
	real pi, convert, x, y, z, sgn, dot, norm

	character ch1*4, ch2*4, ch3*2, ch4*2, casenum*3, char_file, cname2*15
	character gpath*80, lpath*80, opath*80, tpath*80, cname*14

	pi = 3.141592653589793
	convert = 180.0/pi
	ihelp=max(1,2*isamp)
	write(ch3,'(I2.2)') idump
	write(ch4,'(I2.2)') iadj
	
	! in case of a vector field, calculate the amplitude (norm)
	if(isize==3) then
		k=1
		do ir=1,nr+1
		do id=1,nd
		do ii=1,(nt+1)**2			
			x=vec((ir-1)*3*nd*(nt+1)**2+(id-1)*(nt+1)**2+ii)
			y=vec((ir-1)*3*nd*(nt+1)**2+nd*(nt+1)**2+(id-1)*(nt+1)**2+ii)
			z=vec((ir-1)*3*nd*(nt+1)**2+2*nd*(nt+1)**2+(id-1)*(nt+1)**2+ii)
			! if vec points outside, use positive amplitude, else the negative one
			norm=sqrt(x**2.0+y**2.0+z**2.0)
			dot=aint((x*xn(ii,id,1)+y*xn(ii,id,2)+z*xn(ii,id,3))/norm
     &					*100.0)/100.0
			if(dot<0) then
				sgn=-1.0
			else
				sgn=1.0
			endif
			vec_fine(k)=sgn*norm
			k=k+1
		enddo
		enddo
		enddo
	else
		vec_fine=vec
	endif
	
	call layrav(vec_fine,vec_fine_av)

	do ir=1,nr/ihelp+1
		rsh(ir) = rshl(ihelp*(ir-1)+1)
		vec_av(ir)=vec_fine_av(ihelp*(ir-1)+1)
	enddo
	
	if(isamp==0) then
		k=1
		do ir=1,nr/ihelp+1
		do id=1,nd
		do ii=1,(nt/ihelp+1)**2
			if(ir==1) then
				lat(ii,id)=convert*asin(xn(ii,id,3))
				long(ii,id)=convert*atan2(xn(ii,id,2),xn(ii,id,1))
			endif
			if(imean==1) then
				vec_end(k)=(vec_fine(k)/vec_av(ir)-1.0)*100.0
			else
				vec_end(k)=vec_fine(k)
			endif
			k=k+1
		enddo
		enddo
		enddo
	endif
	if(isamp>0) then
		call sample(vec_coar,vec_fine,1,nr,nt)
		call sample(xn_coar,xn,3,0,nt)
		if(isamp==1) then
			k=1
			do ir=1,nr/ihelp+1
			do id=1,nd
			do ii=1,(nt/ihelp+1)**2
				if(ir==1) then
					lat(ii,id)=convert*asin(xn_coar(ii,id,3))
					long(ii,id)=convert*atan2(xn_coar(ii,id,2),xn_coar(ii,id,1))
				endif
				if(imean==1) then
					vec_end(k)=(vec_coar(k)/vec_av(ir)-1.0)*100.0
				else
					vec_end(k)=vec_coar(k)
				endif
				k=k+1
			enddo
			enddo
			enddo
		endif
	endif
	if(isamp>1) then
		call sample(vec_coar2,vec_coar,1,nr/2,nt/2)
		call sample(xn_coar2,xn_coar,3,0,nt/2)
		k=1
		do ir=1,nr/ihelp+1
		do id=1,nd
		do ii=1,(nt/ihelp+1)**2
			if(ir==1) then
				lat(ii,id)=convert*asin(xn_coar2(ii,id,3))
				long(ii,id)=convert*atan2(xn_coar2(ii,id,2),xn_coar2(ii,id,1))
			endif
			if(imean==1) then
				vec_end(k)=(vec_coar2(k)/vec_av(ir)-1.0)*100.0
			else
				vec_end(k)=vec_coar2(k)
			endif
			k=k+1
		enddo
		enddo
		enddo
	endif

	k=1
	write(ch2,'(I4.4)') mynum
	do ir=1,nr/ihelp+1		
		write(ch1,'(I4.4)') (nint((rsh(1)-rsh(ir))/10000.0))*10
		cname=char_file//casenum//'.'//ch1//'.'//ch2
		open(101,file=trim(opath)//'gmt/'//trim(cname), status='replace')

		do id=1,nd
		do ii=1,(nt/ihelp+1)**2
			write(101,*) long(ii,id), lat(ii,id), vec_end(k)
			k=k+1
		enddo
		enddo
		close(101)
	enddo

	call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	if(mynum==0) then
		do ir=1,nr/ihelp+1
			write(ch1,'(I4.4)') (nint((rsh(1)-rsh(ir))/10000.0))*10
			cname2=char_file//casenum//'.'//ch1//'.'//ch3//'_'//ch4
			open(101,file=trim(opath)//'gmt/'//trim(cname2), status='replace')
			
			do ipro=1,nproc
				write(ch2,'(I4.4)') ipro-1
				cname=char_file//casenum//'.'//ch1//'.'//ch2
				open(102,file=trim(opath)//'gmt/'//trim(cname), status='old',
     &					iostat=istat)
				if(istat==0) then
				do id=1,nd
				do ii=1,(nt/ihelp+1)**2
					read(102,*) x,y,z
					write(101,*) x,y,z
				enddo
				enddo
				close(102)
				endif
			enddo
			close(101)
		enddo
		call system('rm '//trim(opath)//'gmt/'//char_file//casenum//'.*.????')
	endif

	end subroutine gmtout2
	
