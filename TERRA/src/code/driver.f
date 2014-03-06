*dk driver
	program terra
 
	include 'size.h'
	include 'pcom.h'
	include 'para.h'

	common /flds/ f(nv*8)
	common /velo/ upb(nt+1), u(3*nv), upe(nt+1)
	common /temp/ tpb(nt+1), temp(nv), tpe(nt+1)
	common /pres/ ppb(nt+1), pres(nv), ppe(nt+1)
	common /name/ titl(4,4)
	common /radl/ rshl(nr+1), ird
	common /radi/ rmin, rmax
	common /mgrd/ itlimit, convtol, itsolve
	common /conv/ step0, stepmin, stepmax
	common /prop/ ieos,  rho0, visc,  grav,  texpn, tcond,
     &              sheat, hgen, tb(2), cl410, cl660
	common /vis1/ vscmax, rvscscl, tvscscl, pwrlawn, pwrlawsr, yldstrs
	common /volc/ tmpmlt
	common /moon/ imoon, radis(5), tcore, qcgn, visc0
	common /io01/ casenum, gpath,  lpath, opath, tpath
	common /io02/ idump0, ngpath, nlpath, vtkw_output, cfile_output

	integer iadj, ieos0, itforw, i,j, itmp, istart
	real tstepadj(itmax+1), tmp, orshl(nr+1)
	character*8  titl, otitl
	character gpath*80, lpath*80, opath*80, tpath*80
	character casenum*3, ctmp, char1*4, char2*3, cname1*15
	logical vtkw_output, cfile_output, check_vtkw_usage

	mynum  = 0
	nproc  = (mt/nt)**2*10/nd
	lvproc = 1.45*log(real(mt/nt))
	mproc  = 2**(2*lvproc)

!	Initialize parallel communications.
	!if(nproc>1) call pinit
	call pinit
 
!	Read input data.
	open(5,file='interra',status='old')
 
	read (5,10) titl, casenum, ird, rmin, rmax, itlimit, convtol,
     & idump0, step0, stepmin, stepmax, ieos, rho0, visc,
     & grav, texpn, tcond, sheat, hgen, tb, cl410, cl660, vscmax,
     & rvscscl, tvscscl, pwrlawn, pwrlawsr, yldstrs, tmpmlt

 10   format(4(4a8/),7x,a3/1(i10/),2(1pe10.3/),i10/,e10.3/1(i10/),
     &   3(e10.3/),i10/17(e10.3/),e10.3)

	write(char1,'(I4.4)') nproc
	write(char2,'(I3.3)') mt
! #### SuperMUC ####
!	gpath = '/scratch/pr32mu/di34feq/output_tmp/'
!	opath = './'
!	tpath = './tomo/'
!	ngpath = 35

! #### Tethys-2g ####
	gpath = './output/'
	opath = '/SCRATCH/horbach/'//casenum//'/'
	tpath = '/SCRATCH/horbach/TERRA'//char1//'_mt'//char2//'/tomo/'
	ngpath = 9
!	gpath='/SCRATCH/horbach/'//casenum//'/c-files/'
!	ngpath = 29

! #### borgcube, buserror ####
!	gpath = './output/'
!	opath = './'//casenum//'/'
!	tpath = './tomo/'
!	tpath = '/import/netapp-m-02-terra/horbach/adjoint/
!     &Endfeld/Endfeld0032_mt032/'
!	ngpath = 9
	
!	determine desired output style
	call get_output_types(5,vtkw_output,cfile_output)
	
	close(5)
 
	lpath  =  gpath
	nlpath = ngpath
 
!	open output files
	if(mynum==0) then
		open(6,file=trim(opath)//'out'//casenum,status='unknown')
		open(8,file=trim(opath)//'outtime'//casenum,status='unknown')
		open(456,file=trim(opath)//'outheat'//casenum,status='unknown')
		open(222,file=trim(opath)//'outmult'//casenum,status='unknown')
		if(nadj>0) then
			open(7,file=trim(opath)//'outadj'//casenum,status='unknown')
			open(223,file=trim(opath)//'outnorm'//casenum,status='unknown')
		endif
     		if(ibc==6) then
     			open(311,file=trim(opath)//'outplate'//casenum,status='unknown')
			open(312,file=trim(opath)//'outrot'//casenum,status='unknown')
		endif
	endif

!	Echo input data to output file.
	if(mynum==0) then
		write(6,20) titl, nproc, nr, nt, ird, ibc,
     &   rmin, rmax, itlimit, convtol, init, itmax, idump0,
     &   step0, stepmin, stepmax, ieos, rho0, visc,
     &   grav, texpn, tcond, sheat, hgen, tb, cl410, cl660,
     &   vscmax, rvscscl, tvscscl, pwrlawn, pwrlawsr, yldstrs,
     &   tbeg, tsim, tmpmlt
	
		if(vtkw_output.and.cfile_output) then
			write(6,'(''EXPTYPE  ='',A11)') 'both'
		else if ( vtkw_output ) then
			write(6,'(''EXPTYPE  ='',A11)') 'vtkw'
		else if ( cfile_output ) then
			write(6,'(''EXPTYPE  ='',A11)') 'c-file'
		else
			write(6,'(''EXPTYPE  ='',A11)') 'none!!!'
		endif
	endif

 20   format(4(4a8/)/,      'NPROC    = ',    i10/,
     &  'NR       = ',  i10/'NT       = ',    i10/'IRD      = ',  i10/,
     &  'IBC      = ',  i10/'RMIN     = ',1pe10.3/'RMAX     = ',e10.3/,
     &  'ITLIMIT  = ',  i10/'CONVTOL  = ',e10.3/,
     &  'INIT    = ',  i10/'ITMAX    = ',    i10/'IDUMP0   = ',  i10/,
     &  'STEP0    = ',  e10.3/,
     &  'STEPMIN  = ',e10.3/'STEPMAX  = ',  e10.3/,
     &  'IEOS     = ',  i10/'RHO0     = ',  e10.3/'VISC     = ',e10.3/,
     &  'GRAV     = ',e10.3/'TEXPN    = ',  e10.3/'TCOND    = ',e10.3/,
     &  'SHEAT    = ',e10.3/'HGEN     = ',  e10.3/'TB(1)    = ',e10.3/,
     &  'TB(2)    = ',e10.3/'CL410    = ',  e10.3/'CL660    = ',e10.3/,
     &  'VSCMAX   = ',e10.3/'RVSCSCL  = ',  e10.3/'TVSCSCL  = ',e10.3/,
     &  'PWRLAWN  = ',e10.3/'PWRLAWSR = ',  E10.3/'YLDSTRS  = ',e10.3/,
     &  'TBEG     = ',e10.3/'TSIM     = ',  e10.3/'TMPMLT   = ',e10.3 )

	rshl(1) = rmax
	rshl(nr+1) = rmin

	! Read old time steps from file in case of restart in adjoint run
!	itforw=346
!	if(mynum==0) then
!		open(666,file=trim(opath)//'outtime'//casenum//'_1',status='old')
!		read(666,*) ctmp
!		do i=1,itforw
!			read(666,*) itmp, tmp, tstepadj(i)
!			tstepadj(i)=tstepadj(i)/3.1688e-8*1000000.0
!			if(i==itforw) tstepadj(i+1)=tmp*1000000.0
!		enddo
!		close(666)
!	endif
!	call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
!	call MPI_BCAST(tstepadj,itmax+1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)
!	call MPI_BCAST(time,1,MPI_REAL8,0,MPI_COMM_WORLD,IERROR)
!	call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	
	if(init<0) then
		istart=-init+1
	else
		istart=1
	endif

	do iadj=istart,nadj+1

!	Perform forward convection calculation.
!	if(iadj>1) then
	if(mynum==0) then
		write(6,*)
		write(6,'(/24x,"BEGIN FORWARD RUN ",i2/)') iadj	
		write(8,'(/24x,"BEGIN FORWARD RUN ",i2/)') iadj
		write(456,*) 1, nadj+1, iadj
	endif
	call convect(tstepadj,iadj,itforw)
!	endif

!	Perform adjoint backward calculation.
		if(iadj<=nadj) then
		
!			write(char1,'(I4.4)') mynum
!			cname1='c608.'//char1//'.01_00'
!			open (99, file=trim(gpath)//cname1, status='unknown')
!			call vecin(temp,otitl,orshl,1,nr,99,1)
!			close(99)
!			call diffout(iadj,iter)
			
			if(mynum==0) then
				write(6,*)
				write(6,'(/24x,"BEGIN ADJOINT RUN ",i2/)') iadj
				write(8,'(/24x,"BEGIN ADJOINT RUN ",i2/)') iadj
			endif
			ieos0 = ieos
			ieos = 1
			call convectadjoint(tstepadj,iadj,itforw)
			ieos = ieos0
		endif
		
	enddo
 
!	Leave PVM before exiting.
	if(mynum==0) then
		write(6,*)
		write(6,*) "Servus!"
		write(6,*)
		close(6)
		if(nadj>0) then
			close(7)
			close(223)
		endif
		close(8)
		close(456)
		close(222)
		if(ibc==6) then
			close(311)
			close(312)
		endif
	endif
 	!if(nproc>1) call pexit
 	call pexit

	end program


! =============================================================================
! Fixed data for moon stuff
! =============================================================================

	blockdata bdterra
 
	common /moon/ imoon, radis(5), tcore, qcgn, visc0

	data imoon/0/, radis/0.033e-6, 0.033e-6, 0.125e-6, 83.e-6, 0./
	data tcore/2500./, qcgn/0.e12/, visc0/2.e20/
 
	end blockdata

