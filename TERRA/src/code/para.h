! ##### TERRA simulation parameters #####

	integer, parameter:: nadj=6				! number of adjoint iterations (0: just forward simulation)
	integer, parameter:: init=0				! type of initialization (-x: restart, 0: tomography, 2: c-file, 17: cube 1-50: see prev. version, 20-49: random field)
	
	integer, parameter:: res_skp_appl=150.0	! residual will not be applied above (km) 150
	integer, parameter:: res_skp_calc=150.0	! residual will not be calculated above (km) 150
	integer, parameter:: dpth_damp=3000.0	! depth, which damping is applied up to (km)
	
	integer, parameter:: adj_skp=4			! use every 'adj_skp' output file (reduction of number of a-files)
	integer, parameter:: time_skp=4			! calculate u every 'time_skp' timestep (performance increase, ind. of adj_skp)
	integer, parameter:: adj_u=0			! use Stokes solver in adjoint iteration (0/1)

	integer, parameter:: dpth_shcut=200.0	! cut depth of shear heating calculation (km)
	
	real, parameter:: tbeg=100.00e+06		! starting point in the past (years)
	real, parameter:: tsim=100.00e+06		! simulation time (years)
	real, parameter:: velfac=1.21			! velocity scaling factor (1.21)

	integer, parameter:: ibc=6				! boundary condition (1: free slip, 6: plates)
	integer, parameter:: buff=4				! number of plate buffer zones
	integer, parameter:: plateskp=3			! how many plate stages should be skipped? (use every 'plateskp+1' stages)

	integer, parameter:: itmax=500000		! maximum number of time steps
	integer, parameter:: npres=10			! number of pressure iterations in forward iteration
	integer, parameter:: npres_adj=1		! number of pressure iterations in backward iteration
		
	integer, parameter:: nout=11			! number of output files incl. first/final step (equal time intervals)
	integer, parameter:: nout_gmt=21		! number of gmt-output files incl. first/final step (equal time intervals)
	integer, parameter:: isamp=2			! downsampling level for gmt plots (max. 2)
	integer, parameter:: nout_surf=41		! number of timesteps for heatflux and surface velocity output
	
! (velfac: 1.21 (normal))


!       411  casenum  three-digit case number used in naming i/o files.
!         1  ird      radial discretization index--values between 1 and 5
! 3.480e+06  rmin     inner radius of spherical shell
! 6.370e+06  rmax     outer radius of spherical shell
!        20  itlimit  maximum number of multigrid iterations
! 1.000e-02  convtol  convergence tolerance for multigrid solver
!        00  idump0   dump number for restart case.
! 1.000e-10  step0    initial time step fraction of advection limit
! 1.000e-14  stepmin  minimum time step fraction of advection limit
! 3.500e-01  stepmax  maximum time step fraction of advection limit
!        22  ieos     index specifying EOS type--1 for Boussinesq case
! 4.500e+03  rho0     reference density
! 1.000e+21  visc     dynamic viscosity
! 1.000e+01  grav     gravitational acceleration
! 2.500e-05  texpn    volume coefficient of thermal expansion
! 3.000e+00  tcond    thermal conductivity
! 1.000e+03  sheat    specific heat at constant volume
! 6.000e-12  hgen     specific radiogenic heat production rate
! 3.000e+02  tb(1)    temperature at outer shell boundary
! 4.200e+03  tb(2)    temperature on inner shell boundary
! 0.000e+06  cl410    Clapeyron slope (dp/dT) for 410 km transition region
!-0.000e+06  cl660    Clapeyron slope (dp/dT) for 660 km transition region
! 1.000e+03  vscmax   maximum value for viscosity variation
!-1.000e-02  rvscscl  radial scaling for viscosity activation energy
! 0.000e-00  tvscscl  tangential scaling for viscosity activation energy
! 0.000e+00  pwrlawn  power-law exponent (zero turns off this feature)
! 3.000e-15  pwrlawsr transition strain rate for power-law rheology
! 0.000e+08  yldstrs  plastic yield stress (zero turns off this feature)
! 0.000e+03  tmpmlt   reference asthenospheric melting temperature
!   both  exptype  sets type data export with fldsout (c-file or vtkw)

