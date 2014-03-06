program truepol_terra
implicit none

include 'mpif.h'

integer i,i2,j,k,l,m,n,begdeg, enddeg, grid, level_terra, n2, n_mod_tpw
integer n_mod, degree_mod, start1, start2, end1, end2, null, dump
integer  begASP, endASP, time, ip1, ip2, imax, ntime, past, disc, disc2
integer imax1, nmax, ns, ms, ipo, ownvisc, nm, nw, Tplot, grid2, islice
integer adjbeg, adjend, iheat, itopo, nheat, dat_fac, ind_B, visc_Brian
integer itrend, ntime_trend, T_dist, nrfac, k_tmp, k_end, shrplot, velplot
integer np, myid, ierr, k_id_beg, k_id_end

integer, allocatable:: current(:), ipf(:)

integer, parameter:: bdry=1				! boundary condition at the surface (free slip=1, no slip=0)
integer, parameter:: nstep=5000

real(8), parameter:: G=6.67e-11			! gravitational constant
real(8), parameter:: rho0=4500			! reference mantle density
real(8), parameter:: rho_c=11.0e3			! density of the core
real(8), parameter:: pi=3.141592			! pi :-)
real(8), parameter:: Emass=5.974d24		! mass of the Earth [kg]
real(8), parameter:: R=6370000.0			! radius of the Earth [m]
real(8), parameter:: const=1.0e25
   
real(8), allocatable:: gc(:,:), gs(:,:), anom(:,:,:), depth(:), rho_mean(:), theta(:), phi(:)
real(8), allocatable:: A(:,:,:), P(:,:,:), kernel(:,:), velokernel(:,:), velo(:,:), geoid(:,:), dat(:,:), geo(:,:), leg(:,:,:)
real(8), allocatable:: layer(:), kernel2(:,:), refgeo(:,:), ann(:,:,:), ann2(:,:,:), geoid0(:,:)
real(8), allocatable:: rho(:), visc(:), inert(:,:,:), tdelay(:), eij(:,:), ref(:,:), velo_abs(:,:,:)
real(8), allocatable:: topsurf(:,:), topcmb(:,:), topsurf_kernel(:,:), topcmb_kernel(:,:), topsurf_coef(:,:), topcmb_coef(:,:)
real(8), allocatable:: T_mean(:), T_anom(:,:,:), geo0(:,:), topsurf0(:,:), trend_time(:)
real(8), allocatable:: spec_dyn(:), spec_etopo(:), dens_min_ges(:), dens_max_ges(:), dens_min(:), dens_max(:)
real(8) P_ab(6,6), C(8,8), B(8), u(8), corr, sum1, sum2, diff, diff1, diff2, dpth, dpth_end, densdpth_beg, densdpth_end
real(8) grav, fac, grav_a, grav_c, cmb, surf, laysize, sphharm, vec(4,3), h, x, yRK(3), sqt, rhofac
real(8) T1, k_T, eval(3), evec(3,3), emax, vprinc(3), omega(nstep,3), begdens, enddens
real(8) temp, temp1, x20, x21, x22, the, ph, r1, r2, amam, anno, ampdyn
real(8) visctmp, visc0, viscDM, viscLIT, ownvisc_ref, min_dens, max_dens, sqpi
real(8) lat1, lat2, long1, long2, sum3, corr2, gridfac1, gridfac2, dens_abs
real(8) dpth_B(15), visc_B(15), depth2, tdiff, velfac, T_abs, age
! ##### love numbers #####
real(8), allocatable:: aaa(:), ap(:), am(:), an(:), bbb(:), ccc(:)
real(8), allocatable:: ahs(:,:), aks(:,:), als(:,:), amode(:)
real(8), allocatable:: freq(:), freq2(:), determ(:), determ2(:), hflux(:,:)
real(8), allocatable:: etopo_grid(:,:), facto(:), facto2(:), sqr(:), T_min_ges(:), T_max_ges(:)
real(8), allocatable:: velo_abscoef(:,:,:), velo_surf(:,:), T_min(:), T_max(:)
real(8), allocatable:: shrh_coef(:,:,:), shrh(:,:,:)
! ##### ETOPO ##########
real(8), allocatable:: rad(:,:), rad2(:,:), theta2(:), theta3(:), phi2(:), phi3(:), etopo(:,:)
real(8), allocatable:: theta_top(:), phi_top(:), topsurf_egrid(:,:), topsurf_egrid_0(:,:)
real(8), allocatable:: theta4(:), phi4(:)

real(8) sk, lat(nstep), vt(nstep-1)
real(8) aky, apm, ddx, ddx1, rap, amplmax, ampl, vmax

 character(len=8) model
 character(len=9) mname, bb, ee, jj
 character(len=3) plot, temp2, vinput, topstep
 character(len=4) timestep, depth_file
 character(len=2) T_out,dat_sh
 character(len=200) path, ownvisc_file
 
 ! MPI initialization
 call MPI_INIT(ierr) 
 call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
 call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)


!*****************INPUT************************************
! In welchem Abstand sollen Daten gelesen werden (Dateiname)
dat_fac=50

! In welchem Abstand sollen T-Plotdaten generiert werden?
T_dist=1

! Wie soll der sh-Grad der Daten sein, die eingelesen werden sollen?
write(dat_sh,'(I2.2)') 50

! Viskosität von Brian Kennett inkl. Interpolation?
visc_Brian=0

! velfac aus TERRA-Rechnungen für Geoid-Trend
velfac=1.21

open(1610,file='./input_geo',status='unknown')

	read(1610,*) temp2
	read(1610,*) path
	read(1610,*) temp2
	read(1610,*) begdeg
	read(1610,*) enddeg
	read(1610,*) dpth
	read(1610,*) dpth_end		
	read(1610,*) null
	read(1610,*) ntime
	read(1610,*) past
	read(1610,*) itrend
	read(1610,*) temp2	
	read(1610,*) begASP
	read(1610,*) endASP
	read(1610,*) viscDM
	read(1610,*) viscLIT
	read(1610,*) visc0
	read(1610,*) ownvisc
	read(1610,*) ownvisc_file
	read(1610,*) ownvisc_ref
	read(1610,*) temp2
	read(1610,*) grid	
	read(1610,*) lat1
	read(1610,*) lat2
	read(1610,*) long1
	read(1610,*) long2
	read(1610,*) temp2
	read(1610,*) Tplot
	read(1610,*) velplot
	read(1610,*) shrplot
	read(1610,*) islice
	read(1610,*) densdpth_beg
	read(1610,*) densdpth_end
	read(1610,*) temp2
	read(1610,*) iheat
	read(1610,*) nheat
	read(1610,*) temp2
	read(1610,*) itopo	
	read(1610,*) level_terra	
	
 close(1610)
 
sqt=sqrt(4.0*pi)
rhofac=1.0
gridfac1=(long2-long1)/(lat2-lat1)
grid2=(grid-1)*gridfac1+1

allocate(tdelay(ntime), inert(ntime,3,3),eij(ntime,6))

write(timestep,'(I4.4)') (ntime-1)*dat_fac
!write(*,*) timestep, ntime
	
! We take the mean densities of the model at 'today' as our reference model.
open(12,file=''//trim(path)//'/'//trim(timestep)//'_'//trim(dat_sh)//'.dat',status='unknown')

	read(12,*) degree_mod, n_mod, dump
		
	allocate(anom(0:degree_mod,-degree_mod:degree_mod,n_mod))
	allocate(T_anom(0:degree_mod,-degree_mod:degree_mod,n_mod))
	allocate(velo_abscoef(0:degree_mod,-degree_mod:degree_mod,n_mod))
	allocate(shrh_coef(0:degree_mod,-degree_mod:degree_mod,n_mod))
	allocate(depth(n_mod),rho_mean(n_mod),visc(n_mod),rho(n_mod))
	allocate(T_mean(n_mod),T_min(n_mod),T_max(n_mod))
	allocate(T_min_ges(n_mod),T_max_ges(n_mod))
	allocate(dens_min(n_mod),dens_max(n_mod))
	allocate(dens_min_ges(n_mod),dens_max_ges(n_mod))
	do j=1,n_mod
     		read(12,*) depth(j)									 ! depth of layer 'j'
     		depth(j)=anint(depth(j)/100)*100 					 ! same rounding like the depths of the layers		
     		!write(*,*) depth(j)
	enddo
	do j=1,n_mod
     		read(12,*) rho_mean(j)        						 ! mean density of layer 'j'
     		rho_mean(j)=rho_mean(j)*rhofac
		rho(j)=rho_mean(j)/rho0
     	enddo
 
 close(12)
 
! The viscosity (already visc/visc_0) profile of the model is read in
! REMARK:
! The absolute value of the viscosity doesn't have any influence on the shape
! of the kernels, only the size of the differences between the layers does matter!
begASP=depth(1)-begASP*1000
endASP=depth(1)-endASP*1000
begdens=depth(1)-densdpth_beg*1000      
enddens=depth(1)-densdpth_end*1000

diff=abs(depth(1)-begASP)
diff1=abs(depth(1)-begdens)
diff2=abs(depth(1)-enddens)

start1=1
start2=1
end1=1
end2=1

do j=1,n_mod
	if(abs(depth(j)-begASP)<=diff) then
		diff=abs(depth(j)-begASP)
		start1=j
	endif
	if(abs(depth(j)-begdens)<=diff1) then
		diff1=abs(depth(j)-begdens)
		start2=j
	endif
	if(abs(depth(j)-enddens)<=diff2) then
		diff2=abs(depth(j)-enddens)
		end2=j
	endif
enddo
diff=abs(depth(1)-endASP)
do j=1,n_mod
	if(abs(depth(j)-endASP)<diff) then
		diff=abs(depth(j)-endASP)
		end1=j
	endif
enddo
do j=1,n_mod
	if(j<start1) then
		visc(j)=viscLIT/visc0
	else if(j>=end1) then
		visc(j)=viscDM/visc0
	else
		visc(j)=1
	endif
	!write(*,*) (6370000.0-depth(j))/1000.0, visc(j)
enddo

!open(99,file='./visc_temp512.dat', status='replace')
! visc. profile from 'visc.dat'
if(ownvisc==1) then
	visc0=ownvisc_ref
	open(111,file=''//trim(ownvisc_file)//'',status='unknown')
	if(visc_Brian==1) then
		do i=1,15
			read(111,*) dpth_B(i), visc_B(i)
		enddo
		do i=1,n_mod
			ind_B=1
			depth2=(depth(1)-depth(i))/1000.0
			do j=1,15
				if(depth2>=dpth_B(j)) ind_B=min(j+1,15)
			enddo
			!if(depth2<700.0) then
			!	visc(i)=visc_B(ind_B)
			!else
				visc(i)=10.0**(log10(visc_B(ind_B-1))+(depth2-dpth_B(ind_B-1))*&
     &					(log10(visc_B(ind_B))-log10(visc_B(ind_B-1)))/(dpth_B(ind_B)-dpth_B(ind_B-1)))
			!endif
			!write(99,'(F8.3,8X,F8.3)') visc(i)/visc0
		enddo
		visc=visc/visc0
	
	else
	
		do i=1,n_mod
			read(111,*) visc(i)
			!write(99,'(F8.3)') visc(i)
		enddo
		
	endif
 	
 	close(111)
 endif

! close(99)

 ! The sph.harm.coefficients of the reference geoid EGM96 are read in.
open(97,file='eigen5c2',status='unknown')

allocate(refgeo(0:enddeg,-enddeg:enddeg))	
do l=0,enddeg
	do m=0,l
		if (m==0) then
			read(97,*) refgeo(l,m)
		else
			read(97,*) refgeo(l,-m),refgeo(l,m)
     		endif
     		if(m==0.and.l==2) refgeo(l,m)=refgeo(l,m) + 1072.618e-6/sqrt(2.0*l+1.0)		! Nakiboglu 1982
     		if(m==0.and.l==4) refgeo(l,m)=refgeo(l,m) - 2.992e-6/sqrt(2.0*l+1.0)			! Nakiboglu 1982
     		!refgeo(l,m)=refgeo(l,m)*(l-1)
		!if(m/=0) refgeo(l,-m)=refgeo(l,-m)*(l-1)
	enddo
enddo
refgeo=refgeo*Emass*G/R/sqrt(2.0*pi)		! geodetic normalization contained in EIGEN data
!refgeo=refgeo/R**3.0*G*Emass
 close(97)
 
allocate(etopo(0:enddeg,-enddeg:enddeg))
 call topocoef(etopo,enddeg)

 
!*****************OUTPUT*************************
if(myid==0) then
open(14, file=''//trim(path)//'/kernel', status='replace')
open(15, file=''//trim(path)//'/velokernel', status='replace')
open(16, file=''//trim(path)//'/geoid',status='replace')
open(17, file=''//trim(path)//'/geoid_coef',status='replace')
open(18, file=''//trim(path)//'/velo_coef',status='replace')
open(19,file=''//trim(path)//'/topo_surf',status='replace')
open(20,file=''//trim(path)//'/topo_cmb',status='replace')
open(21,file=''//trim(path)//'/topo_surf_kernel',status='replace')
open(22,file=''//trim(path)//'/results',status='replace')
open(23,file=''//trim(path)//'/topo_cmb_kernel',status='replace')
open(24, file=''//trim(path)//'/principal',status='replace')
open(25, file=''//trim(path)//'/refgeo',status='replace')
open(26, file=''//trim(path)//'/topo_cmb_coef',status='replace')
open(110, file=''//trim(path)//'/TPW',status='replace')
open(167, file=''//trim(path)//'/heatflux',status='replace')
open(168, file=''//trim(path)//'/velo_surf',status='replace')
open(166, file=''//trim(path)//'/T_cmb',status='replace')

if(Tplot==1) then
open(161, file=''//trim(path)//'/dens_new',status='replace')
open(162, file=''//trim(path)//'/dens_old',status='replace')
open(163, file=''//trim(path)//'/T_new',status='replace')
open(164, file=''//trim(path)//'/T_old',status='replace')
open(203,file=''//trim(path)//'/T_minmax',status='replace')
endif
endif !myid
!*****************Love number*********************
! ### parameters and allocation ###
n_mod_tpw=5			!n_mod
nm=30					! number of modes
nw=180000

ip1=-21
ip2=-7
imax=180
imax1=1800
nmax=(ip2-ip1+1)*imax

allocate(ipf(0:nmax),amode(nm))
allocate(freq(-nw:-1),freq2(-nw:-1))
allocate(determ(-nw:-1),determ2(-nw:-1))
allocate(aaa(0:n_mod_tpw+1),ap(0:n_mod_tpw+1),am(0:n_mod_tpw+1))
allocate(an(0:n_mod_tpw+1),bbb(0:n_mod_tpw+1),ccc(0:n_mod_tpw+1))
allocate(ahs(0:n_mod_tpw,0:nm+1), aks(0:n_mod_tpw,0:nm+1))
allocate(als(0:n_mod_tpw,0:nm+1))
! ####################################

! depth
do i=1,n_mod_tpw
	aaa(i-1)=depth(i)
enddo
aaa(n_mod_tpw-1)=3480000.0		! CMB 
aaa(n_mod_tpw)=1225500.0			! ICB 	 
aaa(n_mod_tpw+1)=0.0				! center 

! density
ap=0.366650E+04

ap(0)=0.0
!ap(1)=0.323235E+04
!ap(2)=0.366650E+04
!ap(3)=0.400367E+04
ap(n_mod_tpw-1)=0.490367E+04
ap(n_mod_tpw)=0.109010E+05
ap(n_mod_tpw+1)=0.128939E+05
!if(n_mod_tpw>5) then
!	do i=4,n_mod_tpw-2
!		fac=(ap(n_mod_tpw-1)-ap(3))/(n_mod_tpw-4)*(i-3)
!		ap(i)=ap(3)+fac
!	enddo
!endif

! rigidity
am=0.222466E+12

am(0)=0.100000E+00 
!am(1)=0.611380E+11
!am(2)=0.916891E+11
!am(3)=0.222466E+12
am(n_mod_tpw-1)=0.238530E+12
am(n_mod_tpw)=0.000000E+00
am(n_mod_tpw+1)=0.164384E+12
!if(n_mod_tpw>5) then
!	do i=4,n_mod_tpw-2
!		fac=(am(n_mod_tpw-1)-am(3))/(n_mod_tpw-4)*(i-3)
!		am(i)=am(3)+fac
!	enddo
!endif

! viscosity
do i=2,n_mod_tpw
	an(i-1)=visc0*visc(i)
enddo
an(0)=0.1
an(n_mod_tpw-1)=visc0*visc(n_mod)
an(n_mod_tpw)=0.0
an(n_mod_tpw+1)=0.100e+14

! Test
aaa(1)=begASP
aaa(2)=endASP
aaa(3)=3771000.0      
aaa(4)=3480000.0 

an(1)=visc0*100
an(2)=visc0
an(3)=an(1)
an(4)=an(1)

! input: aaa, ap
! output: bbb, ccc, apm
 call sources(n_mod_tpw+1,aaa,ap,bbb,ccc,apm)

! input: aaa, ap, am, an, ip1, ip2, imax, nw, ddx, apm, bbb, ccc
! output: determ, freq
ddx=0.05
 call determinant(2,aaa,ap,am,an,n_mod_tpw+1,ip1,ip2,imax,nw,ddx,apm,bbb,ccc,determ,freq)

amode=-1.0e-12
ns=0
ms=0
ipf=0
ddx1=ddx/10.0
aky=86400.0*365.25*1000.0
i=-nmax
do while(i<=-2.and.ns<nm)
	rap=determ(i)/determ(i+1)
	if(rap<0) then
		ms=ms+1
		ipf(ms)=dint(dlog10(dabs(freq(i)/aky))-0.01)-1
		if(ms==1.or.(ipf(ms)/=ipf(ms-1))) then
			ipo=ipf(ms)
			call determinant(2,aaa,ap,am,an,n_mod_tpw+1,ipo,ipo,imax1,nw,ddx1,apm,bbb,ccc,determ2,freq2)
			k=-imax1
			do while(k<=-2.and.ns<nm)
				if(determ2(k)/determ2(k+1)<0) then
					ns=ns+1
					amode(ns)=freq2(k)
					!write(*,*) "ns = ", ns, amode(ns)
				endif
				k=k+1
			enddo	
		endif
	endif
	i=i+1
enddo

! input: aaa,ap,am,an,bbb,ccc,apm,amode
! output: ahs, aks, als
 call love(2,n_mod_tpw+1,nm,aaa,ap,am,an,bbb,ccc,apm,amode,ahs,aks,als)

sk=0.0
do k=1,ns-1
	sk=sk-aks(0,k)/amode(k)
	!write(*,*) sk, aks(0,k)
enddo

k_T=aks(0,nm+1)
T1=abs(sk/k_T)*1000.0

if(myid==0) write(*,*) T1, k_T

!**********************************************************************************************
! The Geoid kernel is computed:
!**********************************************************************************************
allocate(kernel(0:enddeg,n_mod),A(6,6,n_mod),P(6,6,0:n_mod-1))
allocate(velokernel(0:enddeg,n_mod))
allocate(topsurf_kernel(0:enddeg,n_mod),topcmb_kernel(0:enddeg,n_mod))

if(myid==0) then
write(14,*) enddeg, n_mod
write(15,*) enddeg, n_mod
write(21,*) enddeg, n_mod
write(23,*) enddeg, n_mod

A=0
 C=0
P(:,:,0)=0
 cmb=depth(n_mod)
surf=depth(1)

 call gravity(rho*rho0,depth,n_mod,surf,grav_a)
 call gravity(rho*rho0,depth,n_mod,cmb,grav_c)

do i=1,6
	do j=1,6
		if(i==j) P(i,j,0)=1
		if(i==j) A(i,j,n_mod)=1
	enddo
enddo

! Calculation of the propagator matrices and solvation of the LGS containing the boundary conditions
! for each degree of spherical harmonics.
! For l=0, the matrix C is singular, the system of lin. eq. cannot be solved
do l=2,min(enddeg,degree_mod)

! the matrix A containing the fundamental equations is calculated for each dens./visc. layer
! be aware of: depth(1)=surface, depth(n)=CMB
! values at cmb not used (n-1 layers, but n dens./visc. values available)
	do j=1,n_mod-1
		A(1,1,j)=-2
		A(1,2,j)=l*(l+1)
		A(2,1,j)=-1
		A(2,2,j)=1
		A(2,4,j)=1/visc(j)
		A(3,1,j)=12*visc(j)
		A(3,2,j)=-6*l*(l+1)*visc(j)
		A(3,3,j)=1
		A(3,4,j)=l*(l+1)
		A(3,6,j)=-rho(j)
		A(4,1,j)=-6*visc(j)
		A(4,2,j)=2*(2*l*(l+1)-1)*visc(j)
		A(4,3,j)=-1
		A(4,4,j)=-2
		A(4,5,j)=-rho(j)
		A(5,5,j)=1
		A(5,6,j)=1
		A(6,5,j)=l*(l+1)

! the propagator matrices are determined for each layer step and the respective
! products are stored in the field 'P' for later use during the calculation of P(a,b).
		call propma(A(:,:,j),depth(j)/depth(j+1),l,P(:,:,j))
		P(:,:,j)=matmul(P(:,:,j-1),P(:,:,j))
	enddo

	fac=4*pi*G/(2*l+1)

! the system of linear equations Cu=B is set up
	C(2,2)=bdry
	C(3,1)=1
	C(4,2)=1-bdry
	C(5,3)=1
	C(6,4)=1
	C(7,1)=fac*surf/grav_a*rho0
	C(7,3)=1
	C(7,6)=-fac*cmb/grav_c*rho0*((cmb/surf)**l)
	C(7,7)=C(7,6)*rho_c/rho0
	C(8,1)=C(7,1)*((cmb/surf)**(l+1))
	C(8,6)=-fac*cmb/grav_c*rho0
	C(8,7)=1-fac*cmb/grav_c*rho_c
	do k=1,6
		do m=5,8
			n=m-3
			if(m>6) n=n+1
			C(k,m)=-P(k,n,n_mod-1)
		enddo
	enddo
		
	do i=1,n_mod
	
! Since the depths of the mass anomalies correspond to the general layer structure,
! P(a,b) has already been calculated above and been stored in the field 'P'.
		P_ab=P(:,:,i-1)
		call gravity(rho*rho0,depth,n_mod,depth(i),grav)

		do k=1,6
			B(k)=P_ab(k,3)*depth(i)*grav/visc0-P_ab(k,6)*4*pi*depth(i)*depth(i)*G*rho0/visc0
		enddo
		B(7)=rho0/visc0*fac*((depth(i)/surf)**l)*depth(i)**2
		B(8)=rho0/visc0*fac*((cmb/depth(i))**(l-1))*cmb**2

! the system is solved, the spherical harmonic coefficients of the disturbed potential can be found
! in the third component of the solution vector 'u'.		
		call solve_lgs(C,B,8,u)
		kernel(l,i)=u(3)*visc0/rho0/surf/grav_a
		velokernel(l,i)=u(2)*sqrt(l*(l+1.0))
		topsurf_kernel(l,i)=-u(1)*visc0/grav_a/surf/rho0							!????? why rho0 ?????
		!write(*,*) topsurf_kernel(l,i), kernel(l,i)
		topcmb_kernel(l,i)=(u(6)+u(7)/rho0*rho_c)*visc0/cmb/grav_c/rho_c		!????? why rho_c ????
		write(14,'(4e20.10)') real(l),depth(i)/1000,kernel(l,i), visc(i)*visc0
		write(15,*) l,depth(i)/1000,velokernel(l,i)
		write(21,*) l,depth(i)/1000,topsurf_kernel(l,i)
		write(23,*) l,depth(i)/1000,topcmb_kernel(l,i)
	enddo
enddo
 close(14)
 close(15)
 close(21)
 close(23)
endif !myid

!*****************************************************
! geoid calculation and plotting preparation
!*****************************************************
allocate(geoid(0:enddeg,-enddeg:enddeg))
allocate(geoid0(0:enddeg,-enddeg:enddeg))
allocate(topsurf_coef(0:enddeg,-enddeg:enddeg))
allocate(topcmb_coef(0:enddeg,-enddeg:enddeg))
allocate(theta(grid),phi(grid2),theta2(grid),phi2(grid2))
allocate(geo(grid,grid2),leg(0:enddeg,0:enddeg,grid))
allocate(topsurf(grid,grid2),topcmb(grid,grid2))
allocate(geo0(grid,grid2),topsurf0(grid,grid2))
allocate(etopo_grid(grid,grid2))
allocate(spec_dyn(0:enddeg),spec_etopo(0:enddeg))
allocate(velo(0:enddeg,-enddeg:enddeg))
allocate(ann(grid,grid2,n_mod),ref(grid,grid2))
allocate(ann2(grid,grid2,n_mod),velo_surf(grid,grid2))
allocate(velo_abs(grid,grid2,n_mod),shrh(grid,grid2,n_mod))
allocate(facto(0:2*enddeg),facto2(0:2*enddeg),sqr(2*enddeg+1))

laysize=(surf-cmb)/(n_mod-1)
geo=0
leg=0

! Grid definition (corr. to Matlab axes) of the resulting geoid plot
do i=1,grid
	theta(i)=(1.0-(((lat1+(grid-i)*(lat2-lat1)/(grid-1)))+90.0)/180.0)*pi
enddo
do i=1,grid2
	phi(i)=(long1+(i-1)*(long2-long1)/(grid2-1))/180.0*pi
enddo

! convert to 'degrees'
do i=1,grid
	theta2(i)=anint(((grid-i)*(lat2-lat1)/(grid-1)+lat1)*100.0)/100.0
enddo
do i=1,grid2
	phi2(i)=anint(((i-1)*(long2-long1)/(grid2-1)+long1)*100.0)/100.0
enddo

! Calculation of the associated Legendre polynomials using some recurrence formula.
! The Condon-Shortley-phase (-1)**m used by these recurrence formulae will be cancelled again
! during the calculation of the spherical harmonics later (except for the case l=1, m=1)
! Summarised, this means that the associated Leg. polynomials are defined without the Condon-Shortley-phase here.
do i=1,grid
	leg(0,0,i)=1
	leg(1,0,i)=cos(theta(i))
	leg(1,1,i)=-sqrt(1-cos(theta(i))**2)
enddo

if(enddeg>1) then
	do l=2,enddeg
		do m=0,l
			do i=1,grid
				if (m==0) then
					leg(l,m,i)=((2*l-1)*cos(theta(i))*leg(l-1,m,i)-(l-1)*leg(l-2,m,i))/l
				else
					leg(l,m,i)=(leg(l-2,m,i)-(2*l-1)*sqrt(1-cos(theta(i))**2)*leg(l-1,m-1,i))
				endif
			enddo
		enddo
	enddo
endif

! Calculation of the sph.harm. coefficients of the geoid 'geoid(l,m)' and multiplication
! of these coefficients with the resp. sph. harmonics, evaluated at the points of the
! previous defined grid. 'geo' contains the real geoid heights at the grid points.
if(myid==0) then
write(16,*) lat1, lat2, grid, ntime
write(16,*) long1, long2, grid2, past
write(19,*) lat1, lat2, grid, ntime
write(19,*) long1, long2, grid2, past
write(20,*) lat1, lat2, grid, ntime
write(20,*) long1, long2, grid2, past
write(166,*) lat1, lat2, grid, ntime
write(166,*) long1, long2, grid2, past
write(25,*) grid, ntime, past
write(24,*) ntime
write(18,*) level_terra, begdeg, enddeg
write(26,*) enddeg

! temperature output
if(Tplot==1) then
write(161,*) lat1, lat2, grid, end2-start2+1
write(161,*) long1, long2, grid2
write(162,*) lat1, lat2, grid, end2-start2+1
write(162,*) long1, long2, grid2
write(163,*) lat1, lat2, grid, end2-start2+1
write(163,*) long1, long2, grid2
write(164,*) lat1, lat2, grid, end2-start2+1
write(164,*) long1, long2, grid2
endif

write(168,*) grid, ntime, past
write(22,*) trim(path)

do l=0,begdeg-1
	do m=-l,l
		write(18,*) l, m, 0
	enddo
enddo
endif !myid

! all timesteps where data is given
do i=1,ntime
	if(ntime==1) then
		tdelay(i)=0
	else
		tdelay(i)=(i-1)*(past)/(ntime-1.0)*1.0d6
	endif
enddo

! all square root factorials and square roots
facto=1.0
facto2=1.0
do l=1,2*enddeg
	if(l>170) then
		facto2(l)=l*facto2(l-1)
		facto(l)=facto(l-1)
	else
		facto(l)=l*facto(l-1)
	endif
enddo
facto=sqrt(facto)
facto2=sqrt(facto2)

do i=1,2*enddeg+1
	sqr(i)=sqrt(dble(i))
enddo
sqpi=1.0/sqrt(pi)

! ##################################
! read heat flux from output file and write it to disk
! ##################################
adjend=1
if(iheat==0) allocate(hflux(nheat,adjend))
if(iheat==1) then
	write(temp2,'(I3.3)') dump
	open(51,file=''//trim(path)//'/outheat'//temp2//'',status='unknown')
		read(51,*) adjbeg, adjend
		adjend=adjend-adjbeg+1
		write(167,*) adjend, nheat, past
		allocate(hflux(nheat,adjend))
		do k=1,adjend
			if(k>1) read(51,*) temp
			do i=1,nheat
				read(51,*) temp, hflux(i,k), temp
				write(167,*) hflux(i,k), temp
			enddo
		enddo
	close(51)
	close(167)
endif

! #####################################
! read geoid trend time steps from file 'trend_timesteps'
! #####################################
ntime_trend=dat_fac*(ntime-1)+1
if(itrend==1) then
	allocate(trend_time(ntime_trend))
	open(52,file=''//trim(path)//'/trend_timesteps',status='unknown')
	do i=1,ntime_trend
		read(52,*) temp, trend_time(i)
		trend_time(i)=trend_time(i)/velfac
	enddo
	close(52)
endif

! ###############################
! main loop over all time steps
! ###############################
time=(ntime-1)*dat_fac
amam=0

! at most 65 slices as T output and parallelisation
nrfac=(n_mod-1)/64

! k_ids laufen von 1 bis k_end
k_end=floor(real(end2-start2)/real(nrfac))+1
k_id_beg=myid*k_end/np+1
k_id_end=(myid+1)*k_end/np
		
do n=1,ntime
	n2=ntime-n+1
	! age: Ma in the past
	if(ntime==1) then
		age=0.0
	else
		age=past-(n2-1)*past/(ntime-1)
	endif
	write(T_out,'(I2.2)') int(age)

	if(itopo==1.and.myid==0) then
		write(topstep,'(I3.3)') n-1
		open(117,file=''//trim(path)//'/topo'//T_out//'.dat',status='unknown')
	endif
	write(timestep, '(i4.4)') time
	open(12,file=''//trim(path)//'/'//timestep//'_'//trim(dat_sh)//'.dat',status='unknown')

	time=time-dat_fac

	do j=1,n_mod+1
     		read(12,*) temp
	enddo
	do j=1,n_mod
     		read(12,*) rho_mean(j)
     		rho_mean(j)=rho_mean(j)*rhofac
	enddo

	do j=1,n_mod
     		do l=0,degree_mod
			read(12,*) anom(l,0,j)
			! The density anomalies are given in percentual deviation,
			! they have to be multiplied with the average density of the layer 'j'
			! and with the TERRA sph.harm normalisation sqrt(4*pi) (cf. Stacey)
			if((depth(1)-depth(j))/1000.0<dpth.or.(depth(1)-depth(j))/1000.0>dpth_end) then
				anom(l,0,j)=0
			else
				anom(l,0,j)=anom(l,0,j)*rho_mean(j)*sqt
			endif
			if (l>0) then			
     				do m=1,l
     					read(12,*) anom(l,-m,j),anom(l,m,j)
     	     				if((depth(1)-depth(j))/1000.0<dpth.or.(depth(1)-depth(j))/1000.0>dpth_end) then
     						anom(l,m,j)=0
     						anom(l,-m,j)=0
     					else
						anom(l,m,j)=anom(l,m,j)*rho_mean(j)*sqt
						anom(l,-m,j)=anom(l,-m,j)*rho_mean(j)*sqt
					endif
    				enddo
			endif
     		enddo
	enddo
	! ### read in temperature anomalies (in %) #####
	do j=1,n_mod
		read(12,*) T_mean(j)
	enddo
	do j=1,n_mod
		do l=0,degree_mod
			read(12,*) T_anom(l,0,j)
			if (l>0) then			
				do m=1,l
					read(12,*) T_anom(l,-m,j), T_anom(l,m,j)
				enddo
			endif
		enddo
	enddo
	T_anom=T_anom*sqt!*100.0
	if(velplot==1) then
	do j=1,n_mod
     		do l=0,degree_mod
			read(12,*) velo_abscoef(l,0,j)
			if (l>0) then			
     				do m=1,l
     					read(12,*) velo_abscoef(l,-m,j), velo_abscoef(l,m,j)
    				enddo
			endif
     		enddo
	enddo
	velo_abscoef=velo_abscoef*sqt
	endif
	if(shrplot==1) then
		do j=1,n_mod
     		do l=0,degree_mod
			read(12,*) shrh_coef(l,0,j)
			if (l>0) then			
     				do m=1,l
     					read(12,*) shrh_coef(l,-m,j), shrh_coef(l,m,j)
    				enddo
			endif
     		enddo
		enddo
		shrh_coef=shrh_coef*sqt
	endif
	! #################################
	close(12)

	topsurf_coef=0.0
	topsurf=0.0
	topcmb_coef=0.0
	topcmb=0.0
 	geoid=0.0
 	geo=0.0
 	ref=0.0
 	velo=0.0
 	ann=0.0
 	ann2=0.0
 	velo_surf=0.0
 	velo_abs=0.0
 	shrh=0.0
 	
	do l=0,min(degree_mod,enddeg)
		do m=-l,l
			if(l>=begdeg) then
			if(myid==0) then
			do j=1,n_mod
      	 		geoid(l,m)=geoid(l,m)+anom(l,m,j)*kernel(l,j)*laysize
      	 		!geoid(l,m)=geoid(l,m)+surf/grav_a*4.0*pi*G/(2.0*l+1.0)*(depth(j)/surf)**(l+2.0)*anom(l,m,j)*laysize
      	 		velo(l,m)=velo(l,m)+anom(l,m,j)*velokernel(l,j)*laysize
             		topsurf_coef(l,m)=topsurf_coef(l,m)+anom(l,m,j)*topsurf_kernel(l,j)*laysize
       			topcmb_coef(l,m)=topcmb_coef(l,m)+anom(l,m,j)*topcmb_kernel(l,j)*laysize
			enddo
			write(18,*) l, m, velo(l,m)
			endif
			! #### topcmb coefficients in TERRA format ######
			if(n==ntime.and.m>=0) then
				 if(m==0) then
					 write(26,'(2E20.10)') topcmb_coef(l,0)/sqt, 0.0
				else
					 write(26,'(2E20.10)') topcmb_coef(l,-m)/sqt, topcmb_coef(l,m)/sqt
				endif
			endif
			endif !myid
			do i=1,grid
				do j=1,grid2
					if(m==0) then
						sphharm=sqr(2*l+1)*0.5*sqpi*leg(l,m,i)
					elseif(m<0.and.m>-l-1) then
	 					sphharm=sqr(2*l+1)*sqpi/sqr(2)*facto(l+m)/facto(l-m)*(-1)**m*leg(l,-m,i)*&
	 								& facto2(l+m)/facto2(l-m)*cos(m*phi(j));
					elseif(m>0.and.m<l+1) then
	 					sphharm=sqr(2*l+1)*sqpi/sqr(2)*facto(l-m)/facto(l+m)*(-1)**m*leg(l,m,i)*&
	 								& facto2(l-m)/facto2(l+m)*sin(m*phi(j));
					endif
					if(l>=begdeg.and.myid==0) then
					geo(i,j)=geo(i,j)+geoid(l,m)*sphharm
					ref(i,j)=ref(i,j)+refgeo(l,m)*sphharm
					topsurf(i,j)=topsurf(i,j)+topsurf_coef(l,m)*sphharm
					topcmb(i,j)=topcmb(i,j)+topcmb_coef(l,m)*sphharm
					!velo_surf(i,j)=velo_surf(i,j)+velo_abscoef(l,m,1)*sphharm
					
					! n_mod-1 (Schicht über CMB) nur berechnen, wenn es nicht schon unten gemacht wird
					if(Tplot==0) ann2(i,j,n_mod-1)=ann2(i,j,n_mod-1)+T_anom(l,m,n_mod-1)*sphharm
					endif
					
					if(Tplot==1.and.l>=begdeg) then
					if(mod(n-1,T_dist)==0.or.islice==1) then
						do k_tmp=k_id_beg,k_id_end
							k=start2+(k_tmp-1)*nrfac
							ann(i,j,k)=ann(i,j,k)+anom(l,m,k)*sphharm 						! Dichte
							ann2(i,j,k)=ann2(i,j,k)+T_anom(l,m,k)*sphharm					! Temperatur
							velo_abs(i,j,k)=velo_abs(i,j,k)+velo_abscoef(l,m,k)*sphharm		! Geschwindigkeit
							if(shrplot==1) shrh(i,j,k)=shrh(i,j,k)+shrh_coef(l,m,k)*sphharm	! shear heating
						enddo
					endif
					endif
					
					if(l==min(degree_mod,enddeg).and.m==l.and.myid==0) then
						
						! Berechnung Geoid-Trend
						if(itrend==1) then
						if(n==1) then
							geo0(i,j)=geo(i,j)
							topsurf0(i,j)=topsurf(i,j)
						else
							tdiff=trend_time(ntime_trend)-trend_time(ntime_trend-dat_fac*(n-1))
							geo(i,j)=(geo0(i,j)-geo(i,j))/tdiff
							topsurf(i,j)=(topsurf0(i,j)-topsurf(i,j))/tdiff
						endif
						endif
						!if(geo(i,j)>1.0) geo(i,j)=1.0
						!if(geo(i,j)<-3.0) geo(i,j)=-3.0
						
						if(itopo==1) write(117,*) theta2(i), phi2(j), (topsurf(i,j)-geo(i,j))
						
						write(16,"(E30.20)",advance="no") geo(i,j)
						write(19,"(E30.20)",advance="no") topsurf(i,j)
						write(20,"(E30.20)",advance="no") topcmb(i,j)
						write(25,"(E30.20)",advance="no") ref(i,j)
						write(166,"(E30.20)",advance="no") ann2(i,j,n_mod-1)!+1.0)*T_mean(n_mod-1)
						!write(168,"(E30.20)",advance="no") velo_surf(i,j)*10.0e10*sqt
						
						if(amam<abs(geo(i,j)).and.n==1) amam=abs(geo(i,j))
						if(ampdyn<abs(topsurf(i,j)).and.n==1) ampdyn=abs(topsurf(i,j))
						if(j<grid2) then
							write(16,"(A)",advance="no") " "
							write(19,"(A)",advance="no") " "
							write(20,"(A)",advance="no") " "
							write(25,"(A)",advance="no") " "
							write(166,"(A)",advance="no") " "
							!write(168,"(A)",advance="no") " "
						endif
					endif
				enddo
				if(l==min(degree_mod,enddeg).and.m==l.and.myid==0) then
					write(16,*)
					write(19,*)
					write(20,*)
					write(25,*)
					write(166,*)
					!write(168,*)
				endif
			enddo
		enddo
	enddo
	
	! write temperature to file at every time step for slices
	if(islice==1) then
	open(169, file=''//trim(path)//'/T_'//T_out,status='replace')
	write(169,*) grid, end2-start2+1
	do k=start2,end2
		write(169,*) (depth(1)-depth(k))/1000.0
		do i=1,grid
			do j=1,grid2
				write(169,"(E30.20)",advance="no") ann2(i,j,k)	! (ann2(i,j,k)/100.0+1.0)*T_mean(k)
				if(j<grid2) write(169,"(A)",advance="no") " "
			enddo
			write(169,*)
		enddo
		enddo
	close(169)
	endif
	
! density/temperature/velocity/shear heating output
	if(mod(n-1,T_dist)==0.and.Tplot==1) then
		if(n==ntime) then
			disc=162
			disc2=164
		else
			disc=161
			disc2=163
		endif

		do k_tmp=k_id_beg,k_id_end
		k=start2+(k_tmp-1)*nrfac

		write(depth_file,"(I4.4)") (nint((depth(1)-depth(k))/10000.0))*10
		if(myid==0) then
			write(disc,*) depth_file
			write(disc2,*) depth_file
		endif
		
		open(199, file=''//trim(path)//'/plotdata/T.'//depth_file//'.'//T_out//'',status='replace')
		open(202, file=''//trim(path)//'/plotdata/rho.'//depth_file//'.'//T_out//'',status='replace')
		if(velplot==1) then
			open(200, file=''//trim(path)//'/plotdata/vel.'//depth_file//'.'//T_out//'',status='replace')
		endif
		if(shrplot==1) then
			open(201, file=''//trim(path)//'/plotdata/shrh.'//depth_file//'.'//T_out//'',status='replace')
		endif
!		do i=1,grid
!			do j=1,grid2
!				anno=ann(i,j,k)/rho_mean(k)*100.0
!				write(disc,"(E30.20)",advance="no") anno
!				write(disc2,"(E30.20)",advance="no") ann2(i,j,k)		! (ann2(i,j,k)/100.0+1.0)*T_mean(k)	
!				if(j<grid2) write(disc,"(A)",advance="no") " "
!			enddo
!			write(disc,*)
!			write(disc2,*)
!		enddo

		do j=1,grid2
			do i=1,grid
				i2=grid-i+1
				write(199,*) phi2(j), theta2(i2), ann2(i2,j,k)*100.0
				write(202,*) phi2(j), theta2(i2), ann(i2,j,k)/rho_mean(k)*100.0
				if(velplot==1) then
					write(200,*) phi2(j), theta2(i2), velo_abs(i2,j,k)*3.15576e9
				endif
				if(shrplot==1) then
					write(201,*) phi2(j), theta2(i2), shrh(i2,j,k)*3.15576e9
				endif
				T_abs=ann2(i,j,k)*100			!+1.0)*T_mean(k)
				dens_abs=ann(i,j,k)*100			!+1.0)*T_mean(k)
				if(T_abs<T_min(k).or.(j==1.and.i==1)) T_min(k)=T_abs
				if(T_abs>T_max(k).or.(j==1.and.i==1)) T_max(k)=T_abs
				if(dens_abs<dens_min(k).or.(j==1.and.i==1)) dens_min(k)=dens_abs/rho_mean(k)
				if(dens_abs>dens_max(k).or.(j==1.and.i==1)) dens_max(k)=dens_abs/rho_mean(k)
			enddo
		enddo
		
		close(199)
		if(velplot==1) close(200)
		close(202)
		if(shrplot==1) close(201)
		enddo
		
		! temperature min/max per layer
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		call MPI_REDUCE(T_min,T_min_ges,n_mod,MPI_REAL8,MPI_MIN,0,MPI_COMM_WORLD,ierr)
		call MPI_REDUCE(T_max,T_max_ges,n_mod,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		call MPI_REDUCE(dens_min,dens_min_ges,n_mod,MPI_REAL8,MPI_MIN,0,MPI_COMM_WORLD,ierr)
		call MPI_REDUCE(dens_max,dens_max_ges,n_mod,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		
		if(myid==0) then
		write(203,*) T_out//" Ma"
		write(203,*)
		do k_tmp=1,k_end
			k=start2+(k_tmp-1)*nrfac
			write(203,*) (depth(1)-depth(k))/1000.0, T_min_ges(k), T_max_ges(k)
			write(203,*) (depth(1)-depth(k))/1000.0, dens_min_ges(k), dens_max_ges(k)
		enddo
		write(203,*)
		endif
	endif
	
! Calculation of the correlation of the computed and the reference geoid
	if(myid==0) then
	corr=0.0
	corr2=0.0
	sum1=0.0
	sum2=0.0
	
	call spectral(topsurf_coef,2,enddeg,spec_dyn)
	call spectral(etopo,2,enddeg,spec_etopo)
		
	do l=begdeg,enddeg
		do m=-l,l
			corr=corr+geoid(l,m)*refgeo(l,m)
			if(n==1) geoid0(l,m)=geoid(l,m)
			corr2=corr2+geoid(l,m)*geoid0(l,m)
			sum1=sum1+geoid(l,m)**2
			sum2=sum2+refgeo(l,m)**2
		enddo
		!write(*,*) l, spec_dyn(l), spec_etopo(l)
	enddo
	
	if(n==1) sum3=sum1
	corr=corr/sqrt(sum1*sum2)
	corr2=corr2/sqrt(sum1*sum3)
	write(16,*) corr, corr2, hflux(n2,adjend)
	write(25,*) corr, corr2, hflux(n2,adjend)
	write(22,*)
	if(ntime==1) then
		write(*,*) 0, corr, corr2
		write(22,*) 0, corr, corr2
	else
		write(*,*) T_out, corr, corr2
		write(22,*) T_out, corr, corr2
	endif
	
	!if(n==ntime) then
	!	if (enddeg<degree_mod) then
	!		do l=enddeg+1,degree_mod
	!			do m=-l,l
	!				write(18,*) 0
	!			enddo
	!		enddo
	!	endif
	!endif

      ! inertia tensor due to mantle convection
      ! REMARK: the factor 'grav_a' is contained in the equation below because we need here 
      !		    the coef. of the grav.potential and not those of the geoid      
      fac=surf**3.0/G*sqrt(5.0/12.0/pi)*grav_a	
      inert(n,1,1)=fac*(sqrt(1.0/3.0)*geoid(2,0)-geoid(2,-2))
      inert(n,2,2)=fac*(sqrt(1.0/3.0)*geoid(2,0)+geoid(2,-2))
      inert(n,3,3)=-2.0*fac*sqrt(1.0/3.0)*geoid(2,0)
      inert(n,1,2)=-fac*geoid(2,2)
      inert(n,1,3)=-fac*geoid(2,-1)
      inert(n,2,3)=-fac*geoid(2,1)
      inert(n,2,1)=inert(n,1,2)
      inert(n,3,1)=inert(n,1,3)
      inert(n,3,2)=inert(n,2,3)

	! computation of the eigenvectors and values of the inertia tensor
	! the largest eigenvalues corresponds to the maximal principal axis of inertia
	call jacobi(inert(n,:,:),3,eval,evec) 
	
	emax=max(eval(1), eval(2), eval(3))  
	do j=1,3
		if(eval(j)==emax) i=j
	enddo
	do j=1,3
		vprinc(j)=evec(j,i)
	enddo
	if(n==1) write(22,*) vprinc

	if(n==1.and.null==1) then
		omega(1,1)=0
		omega(1,2)=0
		omega(1,3)=1
	endif
	if(n==1.and.null==0) omega(1,:)=vprinc
	
	if(n==1) then
		write(24,*) tdelay(n)
		write(24,*) omega(1,1),omega(1,2), omega(1,3)
	else
		write(24,*) tdelay(n)
		write(24,*) vprinc(1), vprinc(2), vprinc(3)
	endif
	if(itopo==1) close(117)
	endif !myid
enddo

! solution of the equation of motion with Runge Kutta, 'nstep' steps used 
if(myid==0) then

x=tdelay(1)
h=(tdelay(ntime)-tdelay(1))/nstep
vmax=0
amplmax=0
lat(1)=90.0-(pi/2.0-atan(omega(1,3)/sqrt(omega(1,1)**2+omega(1,2)**2)))/pi*180.0
write(110,*) nstep
write(110,*) omega(1,1), omega(1,2), omega(1,3)

do k=1,nstep-1
	call rightHS(x,tdelay,ntime,inert,omega(k,:),k_T,T1,surf,vec(1,:))
	do j=1,3
		if(j==3) then
			fac=h
		else
			fac=h/2.0
		endif
		do i=1,3
			yRK(i)=omega(k,i)+fac*vec(j,i)
		enddo
	     	call rightHS(x+fac,tdelay,ntime,inert,yRK,k_T,T1,surf,vec(j+1,:))
	enddo
	
	do i=1,3
		omega(k+1,i)=omega(k,i)+h/6.0*(vec(1,i)+2.0*(vec(2,i)+vec(3,i))+vec(4,i))
	enddo
	x=x+h
	write(110,*) omega(k+1,1), omega(k+1,2), omega(k+1,3)

	lat(k+1)=90.0-(pi/2.0-atan(omega(k+1,3)/sqrt(omega(k+1,1)**2+omega(k+1,2)**2)))/pi*180.0

	vt(k)=abs(lat(k+1)-lat(k))/h	
	ampl=abs(lat(1)-lat(k+1))
	
	if(vt(k)>vmax) vmax=vt(k)
	if(ampl>amplmax) amplmax=ampl	
enddo

 close(16)
 close(19)
 close(20)
 close(25)
 close(26)

write(*,*)
write(*,*) 'max amplitude geoid= ', amam
write(*,*) 'max amplitude dyntop= ', ampdyn
write(*,*) 'max amplitude TPW= ', amplmax
write(*,*) 'max speed TPW= ', vmax*1000000.0
write(22,*) 'max amplitude geoid= ', amam
write(22,*) 'max amplitude dyntop= ', ampdyn
write(22,*) 'max amplitude TPW= ', amplmax
write(22,*) 'max speed TPW= ', vmax*1000000.0
 
 close(22)
 close(166)
 
if(Tplot==1) then
 close(161)
 close(162)
 close(163)
 close(164)
 close(203)
endif

 endif !myid==0
 
 call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
if(myid==0) then
	write(*,*)
	write(*,*) "Servus!"
endif 
 call MPI_FINALIZE(ierr)

end program truepol_terra


!**********************************************************************************************
!    'propma' calculates the propagator matrix P
!**********************************************************************************************    
subroutine propma(A,rad,l,P)  
implicit none

integer i, j, k, l, id(3)
real(8) P(6,6), A(6,6), uni(6,6), e(4), rad
		
e(1)=l+1
e(2)=-l
e(3)=l-1
e(4)=-l-2
uni=0
P=0

! unit matrix
do i=1,6
	do j=1,6
		if(i==j) uni(i,j)=1
	enddo
enddo

! matrix functional
do i=1,4
	k=1
	do j=1,4
		if(j/=i) then
			id(k)=j
			k=k+1
		endif
	enddo
	P=P+matmul(matmul(A-e(id(1))*uni,A-e(id(2))*uni),A-e(id(3))*uni)&
	&	/(e(i)-e(id(1)))/(e(i)-e(id(2)))/(e(i)-e(id(3)))*rad**e(i)
enddo

end subroutine propma


!**********************************************************************************************
!	'gravity(rho,depth,n,rad,grav)' calculates gravitational acceleration at radius 'rad'
!	within a model consisting of 'n' density layers 'rho' at the respective depth 'depth'
!**********************************************************************************************
subroutine gravity(rho,depth,n,rad,grav)
implicit none

integer n, i
real(8) grav, rad, rho(n), depth(n)

real(8), parameter:: rho_c=11.0e3					! density of the core
real(8), parameter:: fac=4*3.14159/3*6.67e-11		! 4/3*pi*G

i=n-1

! gravitational acc. generated by the core
grav=fac*rho_c*(depth(n)**3)/(rad**2)

! gravitational acc. generated by the density layers below the current position 
do while(depth(i)<rad)
	grav=grav+fac*rho(i)*((depth(i)**3)-(depth(i+1)**3))/(rad**2)	
	i=i-1
enddo

! gravitational acc. generated by the area from the nearest layer boundary to the current position
grav=grav+fac*rho(i)*((rad**3)-depth(i+1)**3)/(rad**2)

end subroutine gravity


!**********************************************************************************************
!	'solve_lgs(A,b,n,x)' solves a system of linear equations Ax=b (with A an nxn matrix)
!	(including partial pivoting)
!**********************************************************************************************
subroutine solve_lgs(A,b,n,x)
implicit none
	
integer n, i, j, k, row
real(8) A_temp(n,n+1), A(n,n), b(n), x(n), s, maxx
	
A_temp(:,1:n)=A
A_temp(:,n+1)=b

! Gauß algorithm including partial pivoting
do i=1,n-1
	maxx=0
	do j=i,n
		s=0
		do k=i,n
			s=s+abs(A_temp(j,k))
		enddo
		if(abs(A_temp(j,i))/s>maxx) then
			maxx=abs(A_temp(j,i))/s
			row=j
		endif
	enddo
	if(maxx==0) then
		write(*,*) "!!!!!! system of linear equations is not solvable  !!!!!!"
		write(*,*) "!!!!!! program aborted !!!!!"
		stop
	endif

	if(row/=i) then
		do j=1,n+1
			maxx=A_temp(i,j)
			A_temp(i,j)=A_temp(row,j)
			A_temp(row,j)=maxx
		enddo
	endif
		
	do j=i+1,n
		A_temp(j,i)=A_temp(j,i)/A_temp(i,i)
		do k=i+1,n+1
			A_temp(j,k)=A_temp(j,k)-A_temp(j,i)*A_temp(i,k)
		enddo
	enddo
enddo

! backsubstitution
do i=1,n
	j=n-i+1
	maxx=0
	if(i>1) then
		do k=j+1,n
			maxx=maxx+x(k)*A_temp(j,k)
		enddo
	endif

	x(j)=(A_temp(j,n+1)-maxx)/A_temp(j,j)

enddo

end subroutine solve_lgs


!*****************************************************************************************
!	eigenvalues and eigenvectors
!	a is the matrix, that you put it (inertiatensor),np i its dimension,  n=3
!	eval(np) are the three eigenvalues, evec(np,np) are the egenvectors for the
!	three eigenvalues ( v(coord, eigenvalues)).
!*****************************************************************************************
subroutine jacobi(Atemp,n,eval,evec)
implicit none

integer, parameter:: nmax=500
integer i, j, k, l, n, nrot
real(8) Atemp(n,n), a(n,n), eval(n), evec(n,n)
real(8) b(nmax), z(nmax), sm, tresh, g, h, t, theta
real(8) c,s,tau

do i=1,n
	do j=1,n
		A(i,j)=Atemp(i,j)
	enddo
enddo

do i=1,n
	do j=1,n
		evec(i,j)=0.0
	enddo
	evec(i,i)=1.0
enddo

do i=1,n
	b(i)=a(i,i)
	eval(i)=b(i)
	z(i)=0.0
enddo

nrot=0
do i=1,50
	
	sm=0.0
	do j=1,n-1
		do k=j+1,n
			sm=sm+abs(a(j,k))
		enddo
	enddo

	if(sm==0.0) return

	if(i<4) then
		tresh=0.2*sm/n**2.0
	else
		tresh=0.0
	endif

	do j=1,n-1
		do k=j+1,n
			g=100.0*abs(a(j,k))
			if((i>4).and.(abs(eval(j))+g==abs(eval(j)))&
				&	.and.(abs(eval(k))+g==abs(eval(k)))) then
				a(j,k)=0.0
			else if(abs(a(j,k))>tresh) then
				h=eval(k)-eval(j)
				if(abs(h)+g==abs(h)) then
					t=a(j,k)/h
				else
					theta=0.5*h/a(j,k)
					t=1.0/(abs(theta)+sqrt(1.0+theta**2.0))
					if(theta<0.0) t=-t
				endif
              		
				c=1.0/sqrt(1.0+t**2.0)
				s=t*c
				tau=s/(1.0+c)
				h=t*a(j,k)
				z(j)=z(j)-h
				z(k)=z(k)+h
				eval(j)=eval(j)-h
				eval(k)=eval(k)+h
				a(j,k)=0.0

				do l=1,j-1
					g=a(l,j)
					h=a(l,k)
					a(l,j)=g-s*(h+g*tau)
					a(l,k)=h+s*(g-h*tau)
				enddo

				do l=j+1,k-1
					g=a(j,l)
					h=a(l,k)
					a(j,l)=g-s*(h+g*tau)
					a(l,k)=h+s*(g-h*tau)
				enddo
				
				do l=k+1,n
					g=a(j,l)
					h=a(k,l)
					a(j,l)=g-s*(h+g*tau)
					a(k,l)=h+s*(g-h*tau)
				enddo
				
				do l=1,n
					g=evec(l,j)
					h=evec(l,k)
					evec(l,j)=g-s*(h+g*tau)
					evec(l,k)=h+s*(g-h*tau)
				enddo
				
				nrot=nrot+1
			endif
		enddo
	enddo

	do j=1,n
          b(j)=b(j)+z(j)
          eval(j)=b(j)
          z(j)=0.0
	enddo
enddo

stop 'too many iterations in jacobi'
return

end subroutine jacobi


!**************************************************************************************************
!	rightHS calculates the right hand side of the diff.equation derived
!	from the Liouville equation. The time x is needed for the values of the time-dependent
!	inertia tensor 'inert' due to mantle convection, which is available for the times 'tdelay' (dim.'ntime'),
!	'omega' is the given rotational axis, 'sol' is the solution vector
!	more parameters needed: k_T: the tidal Love number, T1: the relaxation time, surf: the Earth's radius
!**************************************************************************************************
subroutine rightHS(x,tdelay,ntime,inert,omega,k_T,T1,surf,sol)
implicit none

real(8), parameter:: G=6.67e-11    		    	    	  ! gravitational constant
real(8), parameter:: Emass=5.974d24  			  ! mass of the Earth (kg)
real(8), parameter:: pi=3.141592					  ! just pi ;-)

integer ntime, i, j, m
real(8) x, omega(3), sol(3), w2, fac, MatA(3,3), MatB(3,3), gamm(3), k_T, T1, MatAinv(3,3)
real(8) I_0, dt, tdelay(ntime), E(3,3), erg(3,3), inert(ntime,3,3), surf, fac1, fac2, fac3, xk

w2=omega(1)**2+omega(2)**2+omega(3)**2
fac1=2.0*pi
fac3=2.0*pi/24.0/3600.0

fac=k_T*T1*surf**5.0/3.0/G*fac1*fac3**2.0
I_0=0.33*Emass*surf**2
 
MatA(1,1)=I_0/fac
MatA(2,2)=I_0/fac
MatA(3,3)=I_0/fac

MatA(1,2)=w2*omega(3)
MatA(3,1)=w2*omega(2)
MatA(2,3)=w2*omega(1)
MatA(2,1)=-MatA(1,2)
MatA(1,3)=-MatA(3,1)
MatA(3,2)=-MatA(2,3)

 call invers3(MatA,MatAinv)

! linear interpolation of the inertia tensor to find the matrix B at time t.  
if(ntime>1) then
	do i=1,ntime-1
		if(x>=tdelay(i).and.x<=tdelay(i+1)) j=i 
	enddo
	dt=(x-tdelay(j))/(tdelay(j+1)-tdelay(j))
	do i=1,3
		do m=1,3
			E(i,m)=inert(j,i,m)+(inert(j+1,i,m)-inert(j,i,m))*dt
		enddo
	enddo
	do i=1,3
		MatB(i,i)=(inert(j+1,i,i)-inert(j,i,i))/(tdelay(j+1)-tdelay(j))
	enddo
else
	do i=1,3
		do m=1,3
			E(i,m)=inert(1,i,m)
		enddo
		MatB(i,i)=0
	enddo
endif 

 call matvec(E,omega,3,3,gamm)
MatB(1,2)=gamm(3)
MatB(3,1)=gamm(2)
MatB(2,3)=gamm(1)
MatB(2,1)=-MatB(1,2)
MatB(1,3)=-MatB(3,1)
MatB(3,2)=-MatB(2,3)

do i=1,3
	MatB(i,i)=MatB(i,i)/fac1
enddo
      
erg=-1.0*matmul(MatAinv,2.0*pi*MatB)/fac
 call matvec(erg,omega,3,3,sol)

end subroutine rightHS


!******************************************************************************************
!	invers3 calculates the invers 'InvA' of a 3x3 matrix 'A'
!******************************************************************************************
subroutine invers3(A,InvA)
implicit none

	real(8) A(3,3), InvA(3,3), det
	
	det=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)*A(3,1)-&
		&A(1,2)*A(2,1)*A(3,3)-A(1,1)*A(2,3)*A(3,2)
	! cofactor matrix, adjunct theorem	
	InvA(1,1)=A(2,2)*A(3,3)-A(2,3)*A(3,2)
	InvA(1,2)=A(1,3)*A(3,2)-A(1,2)*A(3,3)
	InvA(1,3)=A(1,2)*A(2,3)-A(1,3)*A(2,2)
	InvA(2,1)=A(2,3)*A(3,1)-A(2,1)*A(3,3)
	InvA(2,2)=A(1,1)*A(3,3)-A(1,3)*A(3,1)
	InvA(2,3)=A(1,3)*A(2,1)-A(1,1)*A(2,3)
	InvA(3,1)=A(2,1)*A(3,2)-A(2,2)*A(3,1)
	InvA(3,2)=A(1,2)*A(3,1)-A(1,1)*A(3,2)
	InvA(3,3)=A(1,1)*A(2,2)-A(1,2)*A(2,1)
	
	InvA=InvA/det
 
end subroutine invers3


!**********************************************************************************************
!	'matvec(A,b,n,m,c)' calculates the product of an matrix A(dim.nxm) and a vector b(dim.m)
!	A*b=c(dim.n)
!**********************************************************************************************
subroutine matvec(A,b,n,m,c)
implicit none

integer:: n,m,i,j
real(8):: A(n,m), b(m), c(n)

do i=1,n
	c(i)=0
	do j=1,m
		c(i)=c(i)+A(i,j)*b(j)
	enddo
enddo

end subroutine matvec


!**********************************************************************************************
!    sources
!**********************************************************************************************  
subroutine sources(n,a,ap,b,c,apm)
implicit none

integer i,k,n
real(8) a(0:n),ap(0:n),b(0:n), c(0:n)
real(8) apm, c1, bm, s

 c1=2.79e-10

!	calcul des contrates de densite
do i=0,n-1
	b(i)=(ap(i+1)-ap(i))/ap(1)
enddo

!	calcul de bm
bm=0.0
do i=0,n-1
	bm=bm+b(i)*((a(i)/a(0))**3)
enddo

apm=bm*ap(1)

!	calcul du champ de gravite en a(i)
do i=0,n-1
	s=0.0
	do k=i+1,n-1
		s=s+(ap(k+1)-ap(k))*(a(k)**3)/(a(i)**2)
	enddo
	c(i)=c1*(ap(i+1)*a(i)+s)
enddo

end subroutine


!**********************************************************************************************
!   determinant
!**********************************************************************************************
subroutine determinant(l,a,ap,am,an,n,ip1,ip2,imax,nw,ddx,apm,b,c,&
     & determ,freq)
implicit none

integer i,j,n,l,nn,ip,ip1,ip2,nw
integer iw, imax, indx(6*n-3)

real(8) G1(6*n-3,6*n-3), G(6*n-3,6*n-3)
real(8) a(0:n),ap(0:n),am(0:n),an(0:n)
real(8) b(0:n),c(0:n)
real(8) freq(-nw:-1),determ(-nw:-1)
real(8) det,dx,ddx,apm,aky,xd,w,d,apo

nn=6*n-3
aky=86400.0*365.25*1000.0
!	introduction de la boucle sur les frequences
do ip=ip1,ip2
	
	dx=ddx*(10d0)**(dfloat(ip))
	xd=(10.)**(dfloat(ip))
	do iw=-imax,-1
		w=(-xd+dx*dfloat(iw))*aky
		call matriceG(l,a,ap,am,an,n,w,b,c,apm,G)
!	introduction de G1 pour sauver les 2n-2 premieres lignes de G
		G1=G
		call ludcmp(G1,nn,indx,d)
		do i=1,nn
			d=d*G1(i,i)
		enddo
		det=d
!	on va multiplier le determinant par le polynome dont les racines sont
!	les temps de Maxwell associes aux couches viscoelastiques du manteau

!	quelque que soit le modele choisi,
!	la lithosphere est la couche 1
!	le manteau est stratifie en couches de 2 a n-2
!	le noyau fluide parfait est la couche n-1
!	la graine est la couche n
		apo=1.0
		do i=1,n-2
			apo=apo*((w+am(i)/an(i)*aky)**4)/w
			if((ip==ip1+1.or.ip==ip2-1).and.(iw==-1)) then
				!write(*,*) apo
			endif
		enddo
				
		if((ip==ip1+1.or.ip==ip2-1).and.(iw==-1)) then
			!write(*,*) apo
		endif
		apo=apo*((w+am(n)/an(n)*aky)**2)/w

		
		if((ip==ip1+1.or.ip==ip2-1).and.(iw==-1)) then
			!write(*,*) apo, det,w
			!write(*,*) "apo: ", apo/(w**2.0)
			!write(*,*) "det: ", apo/(w**2.0)*det
		endif		
		det=(apo)*det/(w**2.0)
!	rentree du vecteur frequence dans un vecteur de double precision 
!	(ip2-ip1)/ddx
		freq((ip1-ip)*imax+iw)=w
		determ((ip1-ip)*imax+iw)=det
	enddo
enddo

end subroutine


!**********************************************************************************************
!    calculation of love numbers
!**********************************************************************************************  
subroutine love(l,n,nm,a,ap,am,an,b,c,apm,amode,ahS,akS,alS)
implicit none

integer i,j,k,l,n,nm,nn,iw
integer indx(6*n-3), indxw(nm)

real(8) a(0:n),ap(0:n),am(0:n),an(0:n)
real(8) c(0:n),b(0:n), G(6*n-3,6*n-3)
real(8) XX(6*n-1), XX1(nm+2), amode(nm)
real(8) ahS(0:n-1,0:nm+1), akS(0:n-1,0:nm+1)
real(8) alS(0:n-1,0:nm+1), amat(nm,nm)
real(8) d,dw,w,apm

nn=6*n-3

!	rentree de la matrice G pour w=amode(i)
!	boucle sur les frequences
do iw=1,nm
	XX=0.0
	XX(3)=real(2.0*l+1.0)
	w=-amode(iw)
	call matriceG(l,a,ap,am,an,n,w,b,c,apm,G)

	indx=0
	call ludcmp(G,nn,indx,d)
	do i=1,nn
		d=d*G(i,i)
	enddo
	call lubksb(G,nn,indx,XX)

!	calcul des deplacements a l interface i pour une surcharge 
	do i=0,n-3
		ahS(i,iw)=a(i)/a(0)*(XX(6*i+1)+XX(6*i+2)+XX(6*i+3)+XX(6*i+4))
		alS(i,iw)=a(i)/a(0)*(-dfloat(l-2)/dfloat(l*(l+1))*XX(6*i+1)-&
     &  1./dfloat(l+1)*XX(6*i+2)+dfloat(l+3)/dfloat(l*(l+1))*XX(6*i+3)+&
     &  1./dfloat(l)*XX(6*i+4))
	enddo
	ahS(n-2,iw)=a(n-2)/a(0)*( ((a(n-3)/a(n-2))**(l+1))*&
	&  XX(6*(n-2)-5)+((a(n-3)/a(n-2))**(l+3))*XX(6*(n-2)-4)+&
	&  ((a(n-2)/a(n-3))**l)*XX(6*(n-2)-3)+((a(n-2)/a(n-3))&
	&  **(l-2))*XX(6*(n-2)-2))
	alS(n-2,iw)=a(n-2)/a(0)*(-((a(n-3)/a(n-2))**(l+1))*&
	&  dfloat(l-2)/dfloat(l*(l+1))*XX(6*(n-2)-5)-1./dfloat(l+1)*&
	&  ((a(n-3)/a(n-2))**(l+3))*XX(6*(n-2)-4)+&
	&  dfloat(l+3)/dfloat(l*(l+1))*((a(n-2)/a(n-3))**l)*XX(6*(n-2)-3)+&
	&  1./dfloat(l)*((a(n-2)/a(n-3))**(l-2))*XX(6*(n-2)-2))
	ahS(n-1,iw)=a(n-1)/a(0)*(XX(6*n-5)+XX(6*n-4))
	alS(n-1,iw)=a(n-1)/a(0)*(dfloat(l+3)/dfloat(l*(l+1))*XX(6*n-5)+&
     &  1./dfloat(l)*XX(6*n-4))

!	calcul du potentiel de redistribution des masses a l interface i
!	pour une surcharge 
	do i=0,n-3
		akS(i,iw)=XX(6*i+5)+XX(6*i+6)-(a(i)/a(0))**l
	enddo
	akS(n-2,iw)=((a(n-2)/a(n-3))**l)*XX(6*(n-2)-1)+&
	&  ((a(n-3)/a(n-2))**(l+1))*XX(6*(n-2))-(a(n-2)/a(0))**l
	akS(n-1,iw)=XX(6*n-3)-(a(n-1)/a(0))**(l)
enddo

do i=1,nm
	do j=1,nm
	amat(i,j)=amode(j)/(amode(i)+amode(j))
	enddo
enddo

 call ludcmp(amat,nm,indxw,dw)
do k=1,nm
	dw=dw*amat(k,k)
enddo

do i=0,n-1

	do k=1,nm
		XX1(k)=ahS(i,k)-ahS(i,0)
	enddo
	call lubksb(amat,nm,indxw,XX1)
	do k=1,nm
		ahS(i,k)=XX1(k)
	enddo

	do k=1,nm
		XX1(k)=alS(i,k)-alS(i,0)
	enddo
	call lubksb(amat,nm,indxw,XX1)
	do k=1,nm
		alS(i,k)=XX1(k)
	enddo

	do k=1,nm
		XX1(k)=akS(i,k)-akS(i,0)
	enddo
	call lubksb(amat,nm,indxw,XX1)
	do k=1,nm
		akS(i,k)=XX1(k)
	enddo
		
!	calcul des limites fluides
	ahS(i,nm+1)=0.0
	alS(i,nm+1)=0.0
	akS(i,nm+1)=0.0
	
	do k=0,nm
		ahS(i,nm+1)=ahS(i,nm+1)+ahS(i,k)
		alS(i,nm+1)=alS(i,nm+1)+alS(i,k)
		akS(i,nm+1)=akS(i,nm+1)+akS(i,k)
		!if(i==0) write(*,*) akS(i,nm+1), akS(i,k)
	enddo
enddo

end subroutine


!**********************************************************************************************
!   matriceG
!**********************************************************************************************  
subroutine matriceG(l,a,ap,am,an,n,w,b,c,apm,G)
implicit none

integer i,k,n,nn,l
integer ik, ik1

real(8) gC1(0:n-1,n),gC2(0:n-1,n),gC3(0:n,n)
real(8) gC4(0:n-1,n),gC5(0:n-1,n),gC6(0:n,n)
real(8) G(6*n-3,6*n-3)
real(8) a(0:n),b(0:n),c(0:n),ap(0:n),am(0:n),an(0:n)
real(8) amu(0:n),delta(0:n)
real(8) apm,w,c1,aky,apg

 c1=2.79e-10	! 4/3*pi*G  
aky=86400.*365.25*1000.
apg=a(0)*c(0)*apm
nn=6*n-3

do i=0,n
	if(i==0.or.i==n-1) then
		amu(i)=am(i)
	else
		if(w>1.0e+05) then
			amu(i)=am(i) 
		else
			amu(i) = am(i)*w/(w + am(i)/an(i)*aky)
		endif
	endif
enddo

delta=0.0
delta(0)=1.0
do i=1,n-3
	do k=1,n
		ik = iabs(k-i)
		ik1 = iabs(k-i-1)
		gC1(i,k) =((a(k-1)/a(k))**(l+1))*delta(ik) - delta(ik1)
		gC2(i,k) =((a(k-1)/a(k))**(l+3))*delta(ik) - delta(ik1)
		gC3(i,k) =((a(k)/a(k-1))**(l))*delta(ik) - delta(ik1)
		gC4(i,k) =((a(k)/a(k-1))**(l-2))*delta(ik) - delta(ik1)
		gC5(i,k)=0.0
		gC6(i,k)=0.0
	enddo
enddo

!	C.L. a la CMB : i=n-2
do k=1,n
	ik=iabs(k-n+2)
	ik1=iabs(k-n+1)
	gC1(n-2,k)=((a(n-3)/a(n-2))**(l+1))*delta(ik)-(c(0)/c(n-2))*&
	&  (a(0)/a(n-2))*delta(ik1)
	gC2(n-2,k)=((a(n-3)/a(n-2))**(l+3))*delta(ik)-(c(0)/c(n-2))*&
	&  (a(0)/a(n-2))*delta(ik1)
	gC3(n-2,k)=((a(n-2)/a(n-3))**(l))*delta(ik)-delta(ik1)
	gC4(n-2,k)=((a(n-2)/a(n-3))**(l-2))*delta(ik)
	gC5(n-2,k)=0.0
	gC6(n-2,k)=0.0
enddo

!	C.L. a l ICB : i=n-1
do k=1,n
	ik=iabs(k-n+1)
	ik1=iabs(k-n)
	gC1(n-1,k)=(c(0)/c(n-1))*(a(0)/a(n-1))*((a(n-1)/a(n-2))**(l))*delta(ik)
	gC2(n-1,k)=(c(0)/c(n-1))*(a(0)/a(n-1))*((a(n-2)/a(n-1))**(l+1))*&
	&  delta(ik)
	gC3(n-1,k)=-delta(ik1)
	gC4(n-1,k)=-delta(ik1)-delta(ik)
	gC5(n-1,k)=0.0
	gC6(n-1,k)=0.0
enddo

!	stockage dans G
do i=1,n-1
	do k=1,n-1
		G(4+6*(i-1),6*(k-1)+1)=gC1(i,k)
		G(4+6*(i-1),6*(k-1)+2)=gC2(i,k)
		G(4+6*(i-1),6*(k-1)+3)=gC3(i,k)
		G(4+6*(i-1),6*(k-1)+4)=gC4(i,k)
		G(4+6*(i-1),6*(k-1)+5)=gC5(i,k)
		G(4+6*(i-1),6*(k-1)+6)=gC6(i,k)
	enddo
	G(4+6*(i-1),6*n-5)=gC3(i,n)
	G(4+6*(i-1),6*n-4)=gC4(i,n)
	G(4+6*(i-1),6*n-3)=gC5(i,n)
enddo

!	Conditions aux limites en deplacement tangentiel
do i=1,n-3
	do k=1,n
		ik = iabs(k-i)
		ik1 = iabs(k-i-1)
		gC1(i,k)=-dfloat(l-2)/dfloat(l*(l+1))*(((a(k-1)/a(k))**(l+1))*delta(ik)&
		&   - delta(ik1) )
		gC2(i,k) = -1./dfloat(l+1)*( ((a(k-1)/a(k))**(l+3))*delta(ik) - &
		&  delta(ik1) )
		gC3(i,k) =dfloat(l+3)/dfloat(l*(l+1))*(((a(k)/a(k-1))**(l))*delta(ik) -&
		&  delta(ik1) )
		gC4(i,k) = 1./dfloat(l)*( ((a(k)/a(k-1))**(l-2))*delta(ik) -&
		&  delta(ik1) )
		gC5(i,k)=0.0
		gC6(i,k)=0.0
	enddo
enddo

!	C.L. a la CMB : i=n-2
do k=1,n
	ik=iabs(k-n+2)
	ik1=iabs(k-n+1)
	gC1(n-2,k)=-dfloat(l-2)/dfloat(l*(l+1))*((a(n-3)/a(n-2))**(l+1))*&
	&  delta(ik)
	gC2(n-2,k)=-1./dfloat(l+1)*((a(n-3)/a(n-2))**(l+3))*delta(ik)
	gC3(n-2,k)=dfloat(l+3)/dfloat(l*(l+1))*((a(n-2)/a(n-3))**(l))*delta(ik)
	gC4(n-2,k)=1./dfloat(l)*((a(n-2)/a(n-3))**(l-2))*delta(ik)
	gC5(n-2,k)=-delta(ik1)
	gC6(n-2,k)=0.0
enddo

!	C.L. a l ICB : i=n-1
do k=1,n
	ik=iabs(k-n+1)
	ik1=iabs(k-n)
	gC1(n-1,k)=0.0
	gC2(n-1,k)=0.0
	gC3(n-1,k)=-dfloat(l+3)/dfloat(l*(l+1))*delta(ik1)
	gC4(n-1,k)=-1./dfloat(l)*delta(ik1)
	gC5(n-1,k)=0.0
	gC6(n-1,k)=delta(ik)
enddo

!	stockage dans G
do i=1,n-1
	do k=1,n-1
		G(5+6*(i-1),6*(k-1)+1)=gC1(i,k)
		G(5+6*(i-1),6*(k-1)+2)=gC2(i,k)
		G(5+6*(i-1),6*(k-1)+3)=gC3(i,k)
		G(5+6*(i-1),6*(k-1)+4)=gC4(i,k)
		G(5+6*(i-1),6*(k-1)+5)=gC5(i,k)
		G(5+6*(i-1),6*(k-1)+6)=gC6(i,k)
	enddo
	G(5+6*(i-1),6*n-5)=gC3(i,n)
	G(5+6*(i-1),6*n-4)=gC4(i,n)
	G(5+6*(i-1),6*n-3)=gC5(i,n)
enddo

!	Conditions aux limites en y4 

!	C.L. en surface : i=0
do	k=1,n
	ik1 = iabs(k-1)
	gC1(0,k) =(2.*amu(k)/apg)*dfloat(l-1)/dfloat(l)*(- delta(ik1) )
	gC2(0,k) = (2.*amu(k)/apg)*dfloat(l+2)/dfloat(l+1)*(-delta(ik1))
	gC3(0,k) = (2.*amu(k)/apg)*dfloat(l+2)/dfloat(l+1)*(-delta(ik1))
	gC4(0,k) = (2.*amu(k)/apg)*dfloat(l-1)/dfloat(l)*(-delta(ik1))
	gC5(0,k)=0.0
	gC6(0,k)=0.0
enddo

!	C.L. dans le manteau : i=1 .. n-3
do i=1,n-3
	do k=1,n
		ik = iabs(k-i)
		ik1 = iabs(k-i-1)
		gC1(i,k)=(2.*amu(k)/apg)*dfloat(l-1)/dfloat(l)*(((a(k-1)/a(k))**(l+1))&
		&  *delta(ik) - delta(ik1) )
		gC2(i,k) = (2.*amu(k)/apg)*dfloat(l+2)/dfloat(l+1)*( ((a(k-1)/a(k))**&
		&  (l+3))*delta(ik) - delta(ik1) )
		gC3(i,k) = (2.*amu(k)/apg)*dfloat(l+2)/dfloat(l+1)*( ((a(k)/a(k-1))**&
		&  (l))*delta(ik) - delta(ik1) )
		gC4(i,k) = (2.*amu(k)/apg)*dfloat(l-1)/dfloat(l)*( ((a(k)/a(k-1))**&
		&  (l-2))*delta(ik) - delta(ik1) )
		gC5(i,k)=0.0
		gC6(i,k)=0.0
	enddo
enddo

!	C.L. a la CMB : i=n-2
do k=1,n
	ik=iabs(k-n+2)
	ik1=iabs(k-n+1)
	gC1(n-2,k)=(2.*amu(k)/apg)*dfloat(l-1)/dfloat(l)*((a(n-3)/a(n-2))**&
	&  (l+1))*delta(ik)
	gC2(n-2,k)=(2.*amu(k)/apg)*dfloat(l+2)/dfloat(l+1)*((a(n-3)/a(n-2))&
	&  **(l+3))*delta(ik)
	gC3(n-2,k)=(2.*amu(k)/apg)*dfloat(l+2)/dfloat(l+1)*((a(n-2)/a(n-3))&
	&  **(l))*delta(ik)
	gC4(n-2,k)=(2.*amu(k)/apg)*dfloat(l-1)/dfloat(l)*((a(n-2)/a(n-3))&
	&  **(l-2))*delta(ik)
	gC5(n-2,k)=0.
	gC6(n-2,k)=0.
enddo

!	C.L. a l ICB : i=n-1
do k=1,n
	ik=iabs(k-n+1)
	ik1=iabs(k-n)
	gC1(n-1,k)=0.
	gC2(n-1,k)=0.
	gC3(n-1,k)=-(2.*amu(k)/apg)*dfloat(l+2)/dfloat(l+1)*delta(ik1)
	gC4(n-1,k)=-(2.*amu(k)/apg)*dfloat(l-1)/dfloat(l)*delta(ik1)
	gC5(n-1,k)=0.
	gC6(n-1,k)=0.
enddo

!	stockage dans G
do i=0,n-1
	do k=1,n-1
		G(1+6*i,6*(k-1)+1)=gC1(i,k)
		G(1+6*i,6*(k-1)+2)=gC2(i,k)
		G(1+6*i,6*(k-1)+3)=gC3(i,k)
		G(1+6*i,6*(k-1)+4)=gC4(i,k)
		G(1+6*i,6*(k-1)+5)=gC5(i,k)
		G(1+6*i,6*(k-1)+6)=gC6(i,k)
	enddo
	G(1+6*i,6*n-5)=gC3(i,n)
	G(1+6*i,6*n-4)=gC4(i,n)
	G(1+6*i,6*n-3)=gC5(i,n)
enddo

!	Conditions aux limites en  y2

!	C.L. en surface : i=0
do k=1,n
	ik1 = iabs(k-1)
	gC1(0,k) =( -(2.*amu(k)/apg)*dfloat(l*l+3*l-1)/dfloat(l+1) +&
	&  ap(1)/apm )*(-delta(ik1) ) 
	gC2(0,k) = ( -(2.*amu(k)/apg)*dfloat(l+2)+ ap(1)/apm )*(-delta(ik1) )
	gC3(0,k) =  ( (2.*amu(k)/apg)*dfloat(l*l-l-3)/dfloat(l)+ap(1)/apm )*&
	&  (-delta(ik1) )
	gC4(0,k) = ( (2.*amu(k)/apg)*dfloat(l-1) +ap(1)/apm )*(-delta(ik1) )
	gC5(0,k)=ap(1)/apm*delta(ik1)
	gC6(0,k)=ap(1)/apm*delta(ik1)
enddo

!	C.L. dans le manteau : i=1 .. n-3
do i=1,n-3
	do k=1,n
		ik = iabs(k-i)
		ik1 = iabs(k-i-1)
		gC1(i,k) =( (2.*amu(k)/apg)*dfloat(-l*l-3*l+1)/dfloat(l+1) +& 
		& (ap(k)*c(i)*a(i)/apg))*(((a(k-1)/a(k))**(l+1))*delta(ik)& 
		&  - delta(ik1) )
		gC2(i,k) =( (2.*amu(k)/apg)*dfloat(-l-2) +(ap(k)*c(i)*a(i)/apg))*&
		& (((a(k-1)/a(k))**(l+3))*delta(ik) - delta(ik1) )
		gC3(i,k) = ( (2.*amu(k)/apg)*dfloat(l*l-l-3)/dfloat(l) +(ap(k)*c(i)*&
		& a(i)/apg) )*( ((a(k)/a(k-1))**(l))*delta(ik) - delta(ik1) )
		gC4(i,k) = ( (2.*amu(k)/apg)*dfloat(l-1) +(ap(k)*c(i)*a(i)/apg) )*&
		& (((a(k)/a(k-1))**(l-2))*delta(ik) - delta(ik1) )
		gC5(i,k) =-ap(k)/apm*( ((a(k)/a(k-1))**(l))*delta(ik) - delta(ik1) )
		gC6(i,k) =-ap(k)/apm*( ((a(k-1)/a(k))**(l+1))*delta(ik) - delta(ik1))
	enddo
enddo

!	C.L. a la CMB : i=n-2
do k=1,n
	ik=iabs(k-n+2)
	ik1=iabs(k-n+1)
	gC1(n-2,k) =( (2.*amu(n-2)/apg)*dfloat(-l*l-3*l+1)/dfloat(l+1) + &
	&  (ap(k)*c(n-2)*a(n-2)/apg) )*((a(k-1)/a(k))**(l+1))*delta(ik)
	gC2(n-2,k) =( (2.*amu(n-2)/apg)*dfloat(-l-2)+(ap(k)*c(n-2)*a(n-2)/apg))&
	&  *((a(k-1)/a(k))**(l+3))*delta(ik)
	gC3(n-2,k) = ( (2.*amu(n-2)/apg)*dfloat(l*l-l-3)/dfloat(l) +&
	&  (ap(k)*c(n-2)*a(n-2)/apg) )*((a(k)/a(k-1))**(l))*delta(ik) - &
	&  (c(n-2)*ap(k)*a(n-2)/apg)*delta(ik1)
	gC4(n-2,k) = ( (2.*amu(n-2)/apg)*dfloat(l-1)+(ap(k)*c(n-2)*a(n-2)/apg))&
	&  *((a(k)/a(k-1))**(l-2))*delta(ik) 
	gC5(n-2,k) =-ap(k)/apm*((a(k)/a(k-1))**(l))*delta(ik)
	gC6(n-2,k) =-ap(k)/apm*((a(k-1)/a(k))**(l+1))*delta(ik)
enddo

!	C.L. a l ICB : i=n-1
do k=1,n
	ik=iabs(k-n+1)
	ik1=iabs(k-n)
	gC1(n-1,k)=0.
	gC2(n-1,k)=0.
	gC3(n-1,k) =-( (2.*amu(n)/apg)*dfloat(l*l-l-3)/dfloat(l) +&
	&  (ap(k)*c(n-1)*a(n-1)/apg) )*delta(ik1)
	gC4(n-1,k) =-(ap(k)*c(n-1)*a(n-1)/apg)*delta(ik) - ( (2.*amu(n)/apg)*&
	&  dfloat(l-1) +(ap(k)*c(n-1)*a(n-1)/apg) )*delta(ik1) 
	gC5(n-1,k) =ap(k)/apm*delta(ik1)
	gC6(n-1,k) =0.
enddo

!	stockage dans G
do i=0,n-1
	do k=1,n-1
		G(2+6*i,6*(k-1)+1)=gC1(i,k)
		G(2+6*i,6*(k-1)+2)=gC2(i,k)
		G(2+6*i,6*(k-1)+3)=gC3(i,k)
		G(2+6*i,6*(k-1)+4)=gC4(i,k)
		G(2+6*i,6*(k-1)+5)=gC5(i,k)
		G(2+6*i,6*(k-1)+6)=gC6(i,k)
	enddo
	G(2+6*i,6*n-5)=gC3(i,n)
	G(2+6*i,6*n-4)=gC4(i,n)
	G(2+6*i,6*n-3)=gC5(i,n)
enddo

!	Conditions aux limites en y5
do i=1,n-3
	do k=1,n
		ik = iabs(k-i)
		ik1 = iabs(k-i-1)
		gC1(i,k)=0.
		gC2(i,k)=0.
		gC3(i,k)=0.
		gC4(i,k)=0.
		gC5(i,k) =((a(k)/a(k-1))**(l))*delta(ik) - delta(ik1)
		gC6(i,k) =((a(k-1)/a(k))**(l+1))*delta(ik) - delta(ik1)
	enddo
enddo

!	C.L. a la CMB : i=n-2
do k=1,n
	ik=iabs(k-n+2)
	ik1=iabs(k-n+1)
	gC1(n-2,k) =- delta(ik1)
	gC2(n-2,k) =- delta(ik1)
	gC3(n-2,k)=0.
	gC4(n-2,k)=0.
	gC5(n-2,k) =((a(k)/a(k-1))**(l))*delta(ik)
	gC6(n-2,k) =((a(k-1)/a(k))**(l+1))*delta(ik)
enddo

!	C.L. a l ICB : i=n-1
do k=1,n
	ik=iabs(k-n+1)
	ik1=iabs(k-n)
	gC1(n-1,k) =((a(k)/a(k-1))**(l))*delta(ik)
	gC2(n-1,k) =((a(k-1)/a(k))**(l+1))*delta(ik)
	gC3(n-1,k)=0.
	gC4(n-1,k)=0.
	gC5(n-1,k) = -delta(ik1)
	gC6(n-1,k) = 0.
enddo

!	stockage dans G
do i=1,n-1
	do k=1,n-1
		G(6+6*(i-1),6*(k-1)+1)=gC1(i,k)
		G(6+6*(i-1),6*(k-1)+2)=gC2(i,k)
		G(6+6*(i-1),6*(k-1)+3)=gC3(i,k)
		G(6+6*(i-1),6*(k-1)+4)=gC4(i,k)
		G(6+6*(i-1),6*(k-1)+5)=gC5(i,k)
		G(6+6*(i-1),6*(k-1)+6)=gC6(i,k)
	enddo
	G(6+6*(i-1),6*n-5)=gC3(i,n)
	G(6+6*(i-1),6*n-4)=gC4(i,n)
	G(6+6*(i-1),6*n-3)=gC5(i,n)
enddo

!	Conditions aux limites en y6 

!	C.L. en surface : i=0  
!	Combinaison de y5 et de y6 : y6(a) +(l+1)/a y5(a)=(2l+1)/a*Vl
do k=1,n
	ik1 = iabs(k-1)
	gC1(0,k) = -3.*ap(1)/apm*delta(ik1)
	gC2(0,k) =-3.*ap(1)/apm*delta(ik1)
	gC3(0,k) =-3.*ap(1)/apm*delta(ik1)
	gC4(0,k) =-3.*ap(1)/apm*delta(ik1)
	gC5(0,k)=dfloat(2*l+1)*delta(ik1)
	gC6(0,k)=0.
enddo

!	C.L. dans le manteau : i=1 .. n-3
do i=1,n-3
	do k=1,n
		ik = iabs(k-i)
		ik1 = iabs(k-i-1)
		gC1(i,k) =-3.*ap(k)/apm*a(k-1)/a(0)* ( ((a(k-1)/a(k))**(l))*delta(ik)&
		&  - delta(ik1) )
		gC2(i,k) =-3.*ap(k)/apm*a(k-1)/a(0)*(((a(k-1)/a(k))**(l+2))*delta(ik)&
		&  - delta(ik1) )
		gC3(i,k) =-3.*ap(k)/apm*a(k-1)/a(0)*( ((a(k)/a(k-1))**(l+1))*delta(ik)&
		&  - delta(ik1) )
		gC4(i,k) =-3.*ap(k)/apm*a(k-1)/a(0)*( ((a(k)/a(k-1))**(l-1))*delta(ik)&
		&  - delta(ik1) )
		gC5(i,k)=a(0)/a(k-1)*dfloat(l)*( ((a(k)/a(k-1))**(l-1))*delta(ik) &
		&  - delta(ik1) )
		gC6(i,k)=-a(0)/a(k-1)*dfloat(l+1)*( ((a(k-1)/a(k))**(l+2))*delta(ik) &
		&  - delta(ik1) )
	enddo
enddo

!	C.L. a la CMB : i=n-2
do k=1,n
	ik=iabs(k-n+2)
	ik1=iabs(k-n+1)
	gC1(n-2,k) =-3.*ap(n-2)/apm*a(n-3)/a(0)*((a(n-3)/a(n-2))**(l))*&
	&  delta(ik)-delta(ik1)*(dfloat(l)*a(0)/a(n-2)-3.*c(0)/c(n-2)*&
	&  ap(n-1)/apm)
	gC2(n-2,k) =-3.*ap(n-2)/apm*a(n-3)/a(0)*((a(n-3)/a(n-2))**(l+2))*&
	&  delta(ik)+delta(ik1)*(dfloat(l+1)*a(0)/a(n-2)+3.*&
	&  c(0)/c(n-2)*ap(n-1)/apm)
	gC3(n-2,k) =-3.*ap(n-2)/apm*a(n-3)/a(0)*((a(n-2)/a(n-3))**(l+1))*&
	&  delta(ik) + 3.*ap(n-1)/apm*a(n-2)/a(0)*delta(ik1)
	gC4(n-2,k) =-3.*ap(n-2)/apm*a(n-3)/a(0)*((a(n-2)/a(n-3))**(l-1))*&
	&  delta(ik)
	gC5(n-2,k)=a(0)/a(n-3)*dfloat(l)*((a(n-2)/a(n-3))**(l-1))*delta(ik) 
	gC6(n-2,k)=-a(0)/a(n-3)*dfloat(l+1)*((a(n-3)/a(n-2))**(l+2))*delta(ik)
enddo

!	C.L. a l ICB : i=n-1
do k=1,n
	ik=iabs(k-n+1)
	ik1=iabs(k-n)
	gC1(n-1,k) = delta(ik)*((a(k)/a(k-1))**l)*(dfloat(l)*a(0)/a(k)-&
	&  3.*c(0)/c(k)*ap(k)/apm)
	gC2(n-1,k) = -delta(ik)*((a(k-1)/a(k))**(l+1))*(dfloat(l+1)*a(0)/a(k)&
	&  +3.*c(0)/c(k)*ap(k)/apm)
	gC3(n-1,k) =3.*ap(n)/apm*a(n-1)/a(0)*delta(ik1) 
	gC4(n-1,k) =3.*ap(n)/apm*a(n-1)/a(0)*delta(ik1) + 3.*ap(n-1)/apm*&
	&  a(n-1)/a(0)*delta(ik)
	gC5(n-1,k)=-a(0)/a(n-1)*dfloat(l)*delta(ik1) 
	gC6(n-1,k)=0.
enddo

!	stockage dans G
do i=0,n-1
	do k=1,n-1
		G(3+6*i,6*(k-1)+1)=gC1(i,k)
		G(3+6*i,6*(k-1)+2)=gC2(i,k)
		G(3+6*i,6*(k-1)+3)=gC3(i,k)
		G(3+6*i,6*(k-1)+4)=gC4(i,k)
		G(3+6*i,6*(k-1)+5)=gC5(i,k)
		G(3+6*i,6*(k-1)+6)=gC6(i,k)
	enddo
	G(3+6*i,6*n-5)=gC3(i,n)
	G(3+6*i,6*n-4)=gC4(i,n)
	G(3+6*i,6*n-3)=gC5(i,n)
enddo

end subroutine


!**********************************************************************************************
!   ludcmp
!**********************************************************************************************
subroutine ludcmp(a,n,indx,d)
implicit none

integer i,j,k,n,imax
integer indx(n)

real, parameter:: tiny=1.0E-20
real(8) a(n,n), vv(n)
real(8) aamax, sum, dum, d

d=1.0

do i=1,n
	aamax=0.0
	do j=1,n
		if(abs(a(i,j))>aamax) aamax=abs(a(i,j))
	enddo
	if(aamax==0.0) stop
	vv(i)=1.0/aamax
enddo

do j=1,n
	do i=1,j-1
		sum=a(i,j)
		do k=1,i-1
			sum=sum-a(i,k)*a(k,j)
		enddo
		a(i,j)=sum
	enddo
	aamax=0.0
	do i=j,n
		sum=a(i,j)
		do k=1,j-1
			sum=sum-a(i,k)*a(k,j)
		enddo
		a(i,j)=sum
		dum=vv(i)*abs(sum)
		if(dum>=aamax) then
			imax=i
			aamax=dum
		endif
	enddo
	if(j/=imax) then
		do	k=1,n
			dum=a(imax,k)
			a(imax,k)=a(j,k)
			a(j,k)=dum
		enddo
		d=-d
		vv(imax)=vv(j)
	endif
	indx(j)=imax
	if(a(j,j)==0.0) a(j,j)=tiny
	if(j/=n) then
		dum=1.0/a(j,j)
		do i=j+1,n
			a(i,j)=a(i,j)*dum
		enddo
	endif
enddo

end subroutine
         
 
!**********************************************************************************************
!   lubksb
!**********************************************************************************************
subroutine lubksb(a,n,indx,b)
implicit none

integer i,j,k,ii,n
integer indx(n)

real(8) a(n,n), b(n+2), sum1

ii=0
do i=1,n
	k=indx(i)
	sum1=b(k)
	b(k)=b(i)
	if(ii/=0) then
		do	j=ii,i-1
			sum1=sum1-a(i,j)*b(j)
		enddo
	else if(sum1/=0) then
		ii=i
	endif
	b(i)=sum1
enddo

do	i=n,1,-1
	sum1=b(i)
	if(i<n) then
		do j=i+1,n
			sum1=sum1-a(i,j)*b(j)
		enddo
	endif
	b(i)=sum1/a(i,i)
enddo

end subroutine

!#####################################
subroutine topocoef(vec,n)
implicit none

integer n,m,l,nmax
real(8) vec(0:n,-n:n), temp

open(97,file='../../Crust/ETOPO/topocoef',status='unknown')

read(97,*) nmax
do l=0,min(n,nmax)
	do m=-l,l
		read(97,*) temp, temp, vec(l,m)	
	enddo
enddo

 close(97) 
 
end subroutine

!#####################################
subroutine spectral(v1,n1,n2, spec)
implicit none

integer l, m, n1, n2
double precision spec(0:n2), v1(0:n2,-n2:n2)

spec=0.0

do l=n1,n2
	do m=-l,l
		spec(l)=spec(l)+v1(l,m)**2
	enddo
	spec(l)=sqrt(spec(l))/(2.0*l+1.0)
enddo

end subroutine

