*dk advance
	subroutine advance(iadj,step,tstep,map,urot,ireport)
!	This routine advances the Navier-Stokes equations one time step.
 
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
	common /radl/ rshl(nr+1), ird
	common /solv/ itreq, convrg, divnorm, cutdiv
	common /mgrd/ itlimit, convtol, itsolve
	common /io01/ casenum, gpath, lpath, opath, tpath
	common /call/ ncall
	common /clck/ itmng, sec(50)

	integer map(0:nt,nt+1,nd), plcount, ireport
	real urot(3,pl_size), step, tstep
	character char1*4, char2*4, cname*14, casenum*3
	character lpath*80, opath*80

	if(itmng==1) call mytime(tin)
	msgnum = 1
	ncall=ncall+1

	if(iadj==1.and.init<0) then
		if(step<=step0.or.ncall==1) plcount=1
	else
		if((step<=step0.and.ncall>min(nr/2,64)).or.ncall==1) plcount=1
	endif
	
!	update operators
	call oprupdate(ibc)
 
!	update the shear heating rate if appropriate.
	if(ieos>=10.and.mod(ncall,4)==1) then
		call shrheat(shr,u,f,f(1,3))
	endif
	
	! convergence output multigrid
	if(ireport==1.and.mynum==0) then
		write(222,*) "CONVECTION STEP ", ncall
		write(222,*)
	endif
		
!	use a second-order Runge-Kutta scheme to update the temperature.
	if(nadj==0.or.real(1.0*ncall/min(nr/2,64))<1.0.or.plcount<=4.or.
     &	mod(ncall,time_skp)==0) then 
		call usolve(map,urot,ireport)
	endif

!	use 'u' coming from usolve and calculate 
!	the temperature change 'tdot'
	call energy(tdot(1,1))
 
	do ii=1,nv
		temp(ii) = temp(ii) + 0.5*tstep*tdot(ii,1)
	enddo
 
	if(nadj==0.or.real(1.0*ncall/min(nr/2,64))<1.0.or.plcount<=4.or.
     &	mod(ncall,time_skp)==0) then 
		call usolve(map,urot,ireport)
 	endif

	call energy(tdot(1,2))
 
	do ii=1,nv
		temp(ii) = temp(ii) + tstep*(tdot(ii,2) - 0.5*tdot(ii,1))
	enddo
 
!-------------------------------------------------------------
!	Compute the new time step!
	call fluxmax(fmax,tstep)
	
!	Small step sizes at the beginning (and at the beginning in case
!	of a restart (iadj==1.and.init<0) because of multigrid.
!	Then, step is fixed to 'stepmax' in the adjoint version, or once
!	set to 'stepmax' in the regular version, except for the begin
!	of a new plate stage (plcount).
	if(ncall>5.or.(iadj==1.and.init<0)) then
		if(itreq>3) then
			step = 0.70*step
		elseif(itreq==1.and.convrg<0.3*convtol) then
			step = min(1.25*step, stepmax)
		endif
		! set timestep ONCE to 'stepmax' after 5 steps of a new plate stage 
		if(real(1.0*ncall/min(nr/2,64))==1.0.or.plcount==5) then
			step=max(step,stepmax)
		endif
		! adj. version: set timestep constant to 'stepmax' (to reduce output)
		! after 5 steps of a new plate stage
		if(nadj>0.and.(real(1.0*ncall/min(nr/2,64))>=1.0
     &		.or.(init<0.and.iadj==1)).and.plcount>4) then
			step=max(step,stepmax)
		endif
		tstep = step*tstep/fmax
	endif
	plcount=plcount+1
!-------------------------------------------------------------

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
      
      
*dk convect
	subroutine convect(tstepadj,iadj,iter)
 
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
	common /heat/ htop, hbot, hrad, heat, hshr, hvol, hnet, tnrm,
     &              tav(nr+1), qc(nr)
	common /work/ wk((nt+1)**2,(nr+1)*2)
	common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
	common /conv/ step0, stepmin, stepmax
	common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
	common /eos2/ rhocv(nr+1), gamma(nr+1), alpha(nr+1), prf(nr+1),
     &               cond(nr+1)
	common /mgrd/ itlimit, convtol, itsolve
	common /vis1/ vscmax, rvscscl, tvscscl, pwrlawn, pwrlawsr, yldstrs
	common /io01/ casenum, gpath, lpath, opath, tpath
	common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output
	common /solv/ itreq, convrg, divnorm, cutdiv
	common /call/ ncall
	common /prty/ propr(20)
	common /name/ titl(4,4)
	common /clck/ itmng, sec(50)

	integer sta, iskp, nplate_tot, oflag, iadj, nout0
	integer map(0:nt,nt+1,nd), begstage
	integer iter, itmp, pout, nout_surf0, oflag_surf, oflag_gmt
	integer ireport, nout_gmt0

	real a, time_stage(0:pl_size,2)
	real urtn(3,pl_size,pl_size)

	real tstepadj(itmax+1), usurf
	real time, step, tstep, tmp, tstep0

	logical vtkw_output, cfile_output
	character*8 titl, ctmp, opath*80

	call mytime(tin)
	call nulvec(sec,50)

	ncall=0
	pout=0
	ireport=0
	tstep0=1.0e-5*visc/(rho0*grav*(rshl(1)-rshl(nr+1)))
	step = step0
	tstep = tstep0

	! set idump to 0 in case of a restart in the following
	! forward iteration
	if(iadj==1) then
		idump = idump0
	else
		idump = 0
	endif

	! time measurement parameter	
	itmng=0
	
	! output variables
	oflag=idump+1
	oflag_surf=oflag
	oflag_gmt=oflag
	nout0=nout
	nout_surf0=nout_surf
	nout_gmt0=nout_gmt
	if(nout0<2) nout0=2
	if(nout_surf0<2) nout_surf0=2
	if(nout_gmt0<2) nout_gmt0=2
	
	! plate variables
	map=0
	urtn=0.0
	begstage=1

	! initialization of function 'randomnum' (very strange!)
	a = randomnum(-1)
	
	! operator initialization
	call opinit

	! equation of state parameter initialization
	! (also ref. temp. is set)
	call eosset

	! field initialization (temp,vel,pres)
	call fldsinit(iadj,iter,time,tstepadj)

	! plate initialization (urtn,time_stage)
	if(ibc==6) then
		call plateinit(time,urtn,time_stage,begstage,nplate_tot)
		call platestage(map,int(time_stage(begstage,2)),nplate_tot)		
		pout=1
	endif
	
!	Begin time step loop
	sta=begstage
	do while(time<tsim*velfac.and.step>=stepmin.and.iter<itmax)
		
		iter=iter+1
		
		tstepadj(iter)=tstep
		time = time + 3.1688e-8*tstep		! seconds-> years
		propr(2) = time
		tp=tbeg*velfac-time	

		if(mynum==0) then
			write(6,*) iter, time/1000000.0,
     &			3.1688e-8*tstep/1000000.0
			write(8,*) iter, time/1000000.0,
     &			3.1688e-8*tstep/1000000.0
		endif
		
		call advance(iadj,step,tstep,map,urtn(:,:,sta),ireport)
			
		! write forward field to file in every adj_skp time step
		! and in the final time step
		! for use in adjoint iteration
		if(nadj>0.and.(mod(iter,adj_skp)==0.or.time>=tsim*velfac).and.
     &			(iadj<=nadj)) then
			call forwardout(iter)
		endif

		if(pout==1) then
			call platemap_out(map,int(time_stage(sta,2)),u,urtn(:,:,sta),
     &					nplate_tot)
			pout=0
		endif
		
		if(tp>0.and.tp<=time_stage(sta-1,1).and.sta>1) then			
			sta=sta-1
			if(mynum==0) then
				write(6,*)
				write(6,'(/24x,"PLATE STAGE ",i7/)') int(time_stage(sta,2))
			endif

			call platestage(map,int(time_stage(sta,2)),nplate_tot)
			
			pout=1
			step = step0
			tstep = tstep0
		endif

		ireport=0
		iskp = 25
		if(iter<100) iskp = 10
		if(mod(iter,iskp)==0) then
			call history2(step,tstep,time)
			ireport=1
		endif
		
		if(mynum==0.and.iter==1) then
			call raylnuss(time)
		endif

		! ###############################
		! ############ OUTPUT #############
		! we have 'nout0' outputs for temp/u/pres in equal time steps
		if(time>=(tsim*velfac)/(nout0-1)*oflag.and.oflag<nout0-1) then
			idump = idump+1
			oflag = oflag+1
			call fldsout2(idump,iadj)		 ! output in separate folder
		endif
		
		! we have 'nout_gmt0' outputs gmt-files in equal time steps
		if(time>=(tsim*velfac)/(nout_gmt0-1)*oflag_gmt.and.
     &			oflag_gmt<nout_gmt0-1) then
			call gmtout2(temp,1,1,oflag_gmt,iadj-1,'t')
			call gmtout2(u,3,1,oflag_gmt,iadj-1,'v')
			oflag_gmt = oflag_gmt+1
		endif
		
		! we have 'nout_surf0' outputs for surf. vel. and CMB heat flux in equal time steps
		if(time>=(tsim*velfac)/(nout_surf0-1)*oflag_surf
     &			.and.oflag_surf<nout_surf0-1) then
			oflag_surf = oflag_surf+1
			call surfvel(usurf,u)
			if(mynum==0) then
				write(456,*) time/velfac/1.0e+06, hbot/1.0e+12, 
     &						   usurf*3.1558e+09
			endif
		endif
		
		! write most recent a-files in separate folder for geoid trend calculations
		if(time>tsim*velfac*0.995.and.nadj>0) call forwardout2(iter)
		! ###############################
		! ###############################
	enddo

	tstepadj(iter+1)=time
	call history2(step,tstep,time)
	call fldsout2(idump+1,iadj)
	call gmtout2(temp,1,1,oflag_gmt,iadj-1,'t')
	call gmtout2(u,3,1,oflag_gmt,iadj-1,'v')
	call surfvel(usurf,u)
	if(mynum==0) then
		write(456,*) (tbeg-time/velfac)/1.0e+06, hbot/1.0e+12, 
     &				   usurf*3.1558e+09
	endif
	! radial averaged velocities
	call uradtng
	
	if(mynum==0) then
		!call radialout
		call raylnuss(time)
	endif
 
!	call wrapup(time)
	if(nadj>0) call diffout(iadj,iter)

	call mytime(tout)
	sec(3) = sec(3) + tout - tin
 
!	call profile
 
	end subroutine


*dk fldsinit
	subroutine fldsinit(iadj,iter,time,tstepadj)
!	This routine obtains initial temperature, velocity and pressure field

	include 'size.h'
	include 'pcom.h'
	include 'para.h'

	common /flds/ f((nt+1)**2*nd*(nr+1)*8)
	common /velo/ upb(nt+1), u((nt+1)**2,nd,3,nr+1),  upe(nt+1)
	common /temp/ tpb(nt+1), temp((nt+1)**2,nd,nr+1), tpe(nt+1)
	common /pres/ ppb(nt+1), pres((nt+1)**2,nd,nr+1), ppe(nt+1)
	common /work/ wk((nt+1)**2*(nr+1)*2)
	common /name/ titl(4,4)
	common /mesh/ xn((nt+1)**2,nd,3)
	common /radl/ rshl(nr+1), ird
	common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
	common /heat/ htop, hbot, hrad, heat, hshr, hvol, hnet, tnrm,
     &              tav(nr+1), qc(nr)
	common /conv/ step0, stepmin, stepmax                                     
	common /io01/ casenum, gpath, lpath, opath, tpath
	common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output
	common /prty/ propr(20)

	integer iadj, i1,i2,ir, nlpath, iter, itmp
	
	real orshl(nr+1), time, tmp, tstepadj(itmax+1)
	
	character*8 titl, otitl(4,4), ctmp
	character gpath*80, lpath*80, cname*12, casenum*3
	character char1*4, char2*3, char3*2, opath*80, tpath*80

	iter = 0
	time = 0.0
	call proprty(iter)

	!	Read in a tomography field in case of using the adjoint method
	!	or when using tomography as initial condition of a forward simulation
	if((iadj==1.and.(nadj>0.or.init==0)).or.
     &			(init<0.and.iadj==-init+1)) then

		write(char1,'(I4.4)') mynum
		write(char2,'(I3.3)') mt
				
		cname='tomo.'//char1//'.00'
		
		open(61, file=trim(tpath)//cname,status='old')
		call vecin(temp,otitl,orshl,1,nr,61,1)
		close(61)
		
		! write temp. field to file as a reference model (final stage)
		if(nadj>0) then
			cname='c511.'//char1//'.01'
			open(61, file=trim(lpath)//cname, status='unknown')
			call vecout(temp,rshl,1,nr,nt,nr,61,1)
			close(61)
		endif
	endif

!	declaration of the initial fields
	if(iadj==1.or.(init<0.and.iadj==-init+1)) then
	!	tomography
		if(init==0) then
			! temp was already read in above
			call nulvec(u,(nt+1)**2*nd*3*(nr+1))
			call nulvec(pres,(nt+1)**2*nd*(nr+1))
	
	! provided initial c-file (also contained in the folder ./tomo/)
		else if(init==2) then
			write(char1,'(I4.4)') mynum
			cname=	'c'//casenum//'.'//char1//'.00'
			open(61, file=trim(tpath)//cname, status='old')
			call vecin(temp,otitl,orshl,1,nr,61,1)
			close(61)
			call nulvec(u,   (nt+1)**2*nd*3*(nr+1))
			call nulvec(pres,(nt+1)**2*nd*(nr+1))
	
	!	specified initial field:
		else if(init>0) then
		
			call tempinit(temp,init)
			call nulvec(u,(nt+1)**2*nd*3*(nr+1))
			call nulvec(pres,(nt+1)**2*nd*(nr+1))

	!	restart:
		else if(init<0) then
		
			write(char1,'(I4.4)') mynum
			write(char3,'(I2.2)') -init
				
			cname='c'//casenum//'.'//char1//'.00'
		
			open(61, file=trim(opath)//'c-files/'//cname//
     &				'_'//char3,status='old')
			call vecin(temp,otitl,orshl,1,nr,61,1)
			close(61)
		
		else
			write(cname,'(''c'',A3,''.'',I4.4,''.'',I2.2)' )
     &			        casenum, mynum, idump0

			open(61, file=trim(lpath)//cname, status='unknown')
			call vecin(temp,otitl,orshl,1,nr,61,1)	! propr is read in here, # of start iter. and time
			call vecin(u   ,otitl,orshl,3,nr,61,0)	! in general, propr is read in by every "vecin"
			call vecin(pres,otitl,orshl,1,nr,61,0)
			close(61)
			
			iter = propr(1)
			time = propr(2)
			
			! Read old time steps from file in case of restart
			if(mynum==0) then
				open(666,file='outtime'//casenum//'_1',status='old')
				read(666,*) ctmp
				do i=1,iter+1
				read(666,*) itmp, tmp, tstepadj(i)
					tstepadj(i)=tstepadj(i)/3.1688e-8*1000000.0
					if(i==iter) time=tmp*1000000.0
				enddo
				close(666)
			endif
			call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
			call MPI_BCAST(tstepadj,itmax+1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
			call MPI_BCAST(time,1,MPI_REAL8,0,MPI_COMM_WORLD,IERROR)
			call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		endif
		
	!	Impose temperature boundary conditions.
		do id=1,nd
			do ii=1,(nt+1)**2
				temp(ii,id,1) = tb(1)
				if(tb(2)>0.0) temp(ii,id,nr+1) = tb(2)
			enddo
		enddo
		
		! write temp. field to file as a reference model (initial stage)
		if(nadj>0) then
			cname='c511.'//char1//'.00'
			open(61, file=trim(lpath)//cname, status='unknown')
			call vecout(temp,rshl,1,nr,nt,nr,61,1)
			close(61)
		endif
		call fldsout(idump0,iter)
		call fldsout2(idump0,iadj)
	else

		write(char1,'(I4.4)') mynum
		cname='c'//casenum//'.'//char1//'.00'
		open(61, file=trim(lpath)//cname, status='unknown')
		call vecin(temp,otitl,orshl,1,nr,61,1)
		close(61)
		call nulvec(u,   (nt+1)**2*nd*3*(nr+1))
		call nulvec(pres,(nt+1)**2*nd*(nr+1))
		
	endif
	
	call layrav(temp,tav)
	call gmtout2(temp,1,1,idump0,iadj-1,'t')
	tnrm  = smean(temp)

	end subroutine


*dk fluxmax
	subroutine fluxmax(fmax,tstep)
!	This routine finds the maximum value for the ratio of flux through
!	a radial cell face to the adjacent cell volume for the purpose of
!	time step control.

	include 'size.h'
	include 'pcom.h'

	common /velo/ upb(nt+1), u((nt+1)**2,nd,3,nr+1),  upe(nt+1)
	common /mesh/ xn((nt+1)**2,nd,3)
	common /radl/ rshl(nr+1), ird
	common /ndar/ arn((nt+1)**2),  arne((nt+1)**2),
     &              rarn((nt+1)**2), rarne((nt+1)**2)
	common /volm/ vol((nt+1)**2,nr+1,2)

	real tstep,fmax

	fmax = 0.0
	do ir=1,nr

		jr = ir + 1
		aa = 0.5*tstep*(0.5*(rshl(ir) + rshl(jr)))**2*arn(2)
     &        *min(vol(2,ir,2), vol(2,jr,2))
 
		do id=1,nd
 
			ii1 = 2
			if(mod(id,5)==1) ii1 = 1
			do ii=ii1,(nt+1)**2
				vfl  = aa*(((u(ii,id,1,ir) + u(ii,id,1,jr))*xn(ii,id,1)
     &                   + (u(ii,id,2,ir) + u(ii,id,2,jr))*xn(ii,id,2))
     &                   + (u(ii,id,3,ir) + u(ii,id,3,jr))*xn(ii,id,3))
				fmax = max(fmax, abs(vfl))
			enddo
		enddo
	enddo

	if(nproc>1) call pmax(fmax)

	end subroutine
      

*dk history2
	subroutine history2(step,tstep,time)
!	This routine writes on output file 6 data describing the time history
!	of the convection calculation.
 
	include 'size.h'
	include 'pcom.h'

	common /velo/ upb(nt+1), u((nt+1)**2,nd,3,nr+1),  upe(nt+1)
	common /heat/ htop, hbot, hrad, heat, hshr, hvol, hnet, tnrm,
     &              tav(nr+1), qc(nr)

	real step, tstep, time

	if(mynum==0) then
		write(6,30)
		write(6,40) (tav(ir),ir=1,nr+1)
 30   format(27x,'MEAN LAYER TEMPERATURES')
 40   format(1x,9f8.1)
	endif

	if(time>1000.0) then
		if(mynum==0) write(6,80) step, 3.1688e-8*tstep, time
 80      format(3x,'ADVECTION STEP =',f7.4,3x,'TIME STEP =',1pe10.3,
     &          3x,'TIME =',e10.3,' YR')

	else
		if(mynum==0) write(6,90) step, 1.1574e-5*tstep, time*365.25
 90      format(3x,'ADVECTION STEP =',f7.4,3x,'TIME STEP =',1pe10.3,
     &          3x,'TIME =',e10.3,' DAYS')

	endif

	call surfvel(usurf,u)
	if(mynum==0) write(6,990) usurf, 3.1558e9*usurf	! years -> seconds and meters -> centimeters
990   format (16x,'RMS  SURFACE VELOCITY     ',1p,e10.3,' M/SEC ',
     &       /16x,'                          ',1p,e10.3,' CM/YR ')

	end subroutine


*dk raylnuss
      subroutine raylnuss(time)
 
c     This routine computes the Rayleigh and Nusselt numbers.
 
      include 'size.h'
      real bmr(nr+1), rcp(nr+1)
      common /radl/ rshl(nr+1), ird
      common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
      common /heat/ htop, hbot, hrad, heat, hshr, hvol, hnet, tnrm,
     &              tav(nr+1), qc(nr)
      common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
      common /eos2/ rhocv(nr+1), gamma(nr+1), alpha(nr+1), prf(nr+1),
     &               cond(nr+1)
      common /vis2/ rdvsc(nr+1), tactv(nr+1), vscl(nr+1)
      common /rayl/ rayl, anuss
 
	real time

      rmax = rshl(1)
      rmin = rshl(nr+1)
      tadb = 0.
 
      if(ieos.gt.5) then
 
         rl     = 0.
         rlh    = 0.
         sheat  = 0.
         texpn  = 0.
         tcond  = 0.
 
         do ir=1,nr+1
            rcp(ir) = rhocv(ir)*(1.+alpha(ir)*gamma(ir)*tav(ir))
            rcp(ir) = rcp(ir)/rdvsc(ir)
         end do
 
         do ir=1,nr
            rt = rshl(ir)
            rb = rshl(ir+1)
            r3 = (rt**3 - rb**3)/3.
            r4 = (rt**4 - rb**4)/4.
            at = alpha(ir  )*grv(ir  )*rhorf(ir  )*rcp(ir  )/cond(ir  )
            ab = alpha(ir+1)*grv(ir+1)*rhorf(ir+1)*rcp(ir+1)/cond(ir+1)
            bt = rhorf(ir  )/cond(ir  )*at
            bb = rhorf(ir+1)/cond(ir+1)*ab
            ct = rhocv(ir  )/rhorf(ir  )
            cb = rhocv(ir+1)/rhorf(ir+1)
            et = alpha(ir  )
            eb = alpha(ir+1)
            rdr   = 1./(rt - rb)
            rl    = rl    + (r4*(at - ab) + r3*(ab*rt - at*rb))*rdr
            rlh   = rlh   + (r4*(bt - bb) + r3*(bb*rt - bt*rb))*rdr
            sheat = sheat + (r4*(ct - cb) + r3*(cb*rt - ct*rb))*rdr
            texpn = texpn + (r4*(et - eb) + r3*(eb*rt - et*rb))*rdr
            tcond = tcond + 0.5*(cond(ir) + cond(ir+1))*(rt - rb)
         end do
 
         rvol  = 3./(rmax**3 - rmin**3)
         rl    = rvol*rl
         rlh   = rvol*rlh
         sheat = rvol*sheat
         texpn = rvol*texpn
         tcond = tcond/(rmax - rmin)
 
         dx      = 0.05
         bmr(1) = (1. + alpha( 1)*gamma( 1)*tav( 1))/drhdp( 1)
 
         do ir=2,nr+1
            bmr(ir) = (1. + alpha(ir)*gamma(ir)*tav(ir))/drhdp(ir)
            dr      = (rshl(ir-1) - rshl(ir))*dx
            x       = -0.5*dx
            do i=1,20
               x       =  x + dx
               omx     = 1. - x
               tadb    = tadb + dr*(x*gamma(ir-1) + omx*gamma(ir))
     &                            *(x*  grv(ir-1) + omx*  grv(ir))
     &                            *(x*  tav(ir-1) + omx*  tav(ir))
     &                            /(x*  bmr(ir-1) + omx*  bmr(ir))
            end do
         end do
 
      else
 
         vc  = 0.
 
         do ir=1,nr
            rt = rshl(ir)
            rb = rshl(ir+1)
            r3 = (rt**3 - rb**3)/3.
            r4 = (rt**4 - rb**4)/4.
            at = rdvsc(ir  )
            ab = rdvsc(ir+1)
            vc = vc + (r4*(at - ab) + r3*(ab*rt - at*rb))/(rt - rb)
         end do
 
         rv  = (rmax**3 - rmin**3)/(3.*vc)
         rl  = texpn*grav*rho0**2*sheat*rv/tcond
         rlh = rl*rho0/tcond
 
      endif
 
      if(ieos/10 .eq. 2) tadb = 1200.
 
      if(tb(2).ne.0.) then
         rayl  = rl*(tb(2) - tb(1) - tadb)*(rmax - rmin)**3/visc
         anuss = htop*(rmax - rmin)/(12.56637*tcond*rmax*rmin
     &           *(tb(2) - tb(1)))
      else
         rayl  = rlh*hgen*(rmax - rmin)**5/visc
         eta   = rmin/rmax
         anuss = hgen*rho0*rmax**2*(1. - eta**2*(3. - 2.*eta))/
     &           (6.*tcond*(tav(nr) - tb(1)))
      endif
 
      write(6,10) anuss, rayl, tadb, time
 10   format(/6x,'NUSSELT NUMBER =',1pe10.3,6x,'RAYLEIGH NUMBER =',
     &   e10.3/16x,'TADB =',e10.3,9x,'ELAPSED TIME =',e10.3,' YR'/)
 
      write(6,20) (rdvsc(ir),ir=1,nr+1)
 20   format(22x,'RADIAL VISCOSITY FACTORS:',9(/1x,8f9.3))
 
	end subroutine


*dk wrapup
      subroutine wrapup(time)
 
      include 'size.h'
      include 'pcom.h'
      common /velo/ upb(nt+1), u((nt+1)**2,nd,3,nr+1),  upe(nt+1)
      common /temp/ tpb(nt+1), temp((nt+1)**2,nd,nr+1), tpe(nt+1)
      common /radl/ rshl(nr+1), ird
      common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
      common /heat/ htop, hbot, hrad, heat, hshr, hvol, hnet, tnrm,
     &              tav(nr+1), qc(nr)
      common /eos1/ rhorf(nr+1), tmprf(nr+1), drhdp(nr+1), grv(nr+1)
      common /eos2/ rhocv(nr+1), gamma(nr+1), alpha(nr+1), prf(nr+1),
     &               cond(nr+1)
      common /name/ titl(4,4)
      common /rayl/ rayl, anuss
      character*8 titl
      real tempav(nr+1), time
 
      if(mynum .eq. 0) write(6,10) titl(1,1), rshl(1), rshl(nr+1), rho0,
     &   visc, texpn, tcond, sheat, hgen, grav, tb
 10   format(1h1,28x,a8,' SUMMARY'//32x,'INPUT DATA'/
     &/9x,'OUTER RADIUS OF SHELL              ',1pe10.3,' M',
     &/9x,'INNER RADIUS OF SHELL              ',1pe10.3,' M',
     &/9x,'REFERENCE DENSITY                  ',1pe10.3,' KG/M**3',
     &/9x,'DYNAMIC VISCOSITY                  ',1pe10.3,' PA-SEC',
     &/9x,'COEFFICIENT OF THERMAL EXPANSION   ',1pe10.3,' 1/K',
     &/9x,'THERMAL CONDUCTIVITY               ',1pe10.3,' W/M/K',
     &/9x,'SPECIFIC HEAT                      ',1pe10.3,' J/KG/K',
     &/9x,'RADIOGENIC HEAT PRODUCTION RATE    ',1pe10.3,' W/KG',
     &/9x,'GRAVITATIONAL ACCELERATION         ',1pe10.3,' M/SEC*2',
     &/9x,'TEMPERATURE OF OUTER BOUNDARY      ',1pe10.3,' K',
     &/9x,'TEMPERATURE OF INNER BOUNDARY      ',1pe10.3,' K')
 
      if(ieos.ge.10 .and. mynum.eq.0) write(6,30)
 30   format(/14x,'SHELL MATERIAL IS TREATED AS A COMPRESSIBLE,',
     &/14x,'ISOTROPIC, HOMOGENEOUS, LINEAR VISCOUS FLUID',
     &/15x,'WITH A MORSE POTENTIAL EQUATION OF STATE.')
      if(ieos.le.5 .and. mynum.eq.0) write(6,40)
 40   format(/13x,'SHELL MATERIAL IS TREATED AS AN INCOMPRESSIBLE,',
     &/14x,'ISOTROPIC, HOMOGENEOUS, LINEAR VISCOUS FLUID.')
 
      if(mynum .eq. 0) write(6,50)
 50   format(/14x,'SHELL BOUNDARIES ARE TREATED AS UNDEFORMABLE,',
     &/21x,'ISOTHERMAL, AND TRACTION-FREE.')
 
      if(hgen.eq.0. .and. mynum.eq.0) write(6,60)
 60   format(/12x,'SHELL IS HEATED STRICTLY FROM ITS INNER BOUNDARY.')
      if(hgen.gt.0. .and. mynum.eq.0) write(6,70)
 70   format(/9x,'SHELL IS HEATED PARTIALLY FROM INTERNAL ',
     &   'RADIOACTIVITY.')
 
      if(mynum .eq. 0) write(6,80) (10*mt**2+2)*(nr+1), 20*mt**2*nr, nr
 80   format(/16x,'GRID DERIVED FROM REGULAR ICOSAHEDRON HAS',
     &/24x,i10,' NODES'/24x,i10,' ELEMENTS',
     &/21x,i3,' RADIAL LAYERS OF ELEMENTS'/)
 
      call surfvel(usurf,u)
      tngarea = 12.56637*rshl(1)**2
      if(mynum .eq. 0) write(6,90) rayl, anuss, usurf, 3.1558e9*usurf,
     &   htop/tngarea, time
 90   format(/30x,'OUTPUT DATA'/
     &/16x,'RAYLEIGH NUMBER           ',1pe10.3,
     &/16x,'NUSSELT NUMBER            ',1pe10.3,
     &/16x,'RMS  SURFACE VELOCITY     ',1pe10.3,' M/SEC ',
     &/16x,'                          ',1pe10.3,' CM/YR ',
     &/16x,'MEAN SURFACE HEAT FLUX    ',1pe10.3,' W/M2  ',
     &/16x,'TOTAL ELAPSED RUN TIME    ',1pe10.3,' YR    ')
 
      lvr = 1.45*log(real(nr))
      call layrav(temp,tempav)
      if(mynum .eq. 0) write(6,120)
 120  format(/60x,'RADIAL'/32x,'LAYER AVERAGES',12x,'CONDUCTIVE',
     &/16x,'DEPTH    DENSITY   PRESSURE  TEMPERATURE   HEAT FLOW',
     &/16x,' (KM)    (KG/M3)     (PA)        (K)          (W)')
 
      do ir=1,nr+1
      depth = rshl(1) - rshl(ir)
      if(mynum .eq. 0) write(6,140) depth,rhorf(ir),prf(ir),tempav(ir)
 140  format(5x,'---------- ',-3pf6.1,1p3e11.3,' -------------')
      if(ir.le.nr .and. mynum.eq.0) write(6,160) ir, qc(ir)
 160  format(6x,'LAYER ',i3,42x,1pe11.3)
      end do
 
	end subroutine

