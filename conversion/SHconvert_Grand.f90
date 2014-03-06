program SHconvert_Grand
implicit none

include 'size.h'
include 'pcom.h'

integer xxx, psp, byebye, idump, pro, lmax_tomo
integer i,j,k,ir, id, ii,locsize, maxls, tzt, nproc2, kt
integer ind(nr+1), ind_mod(nr+1)

integer, parameter:: nn=(nt+1)**2*nd      
integer, parameter:: lmax=50
integer, parameter:: nr_tomo=128

real, parameter:: PI = 3.141592653589793
real, parameter:: RHOAV = 5514.3d0		 ! average density in the full Earth to normalize equation
real, parameter:: GRAV = 6.6723d-11		 ! gravitational constant

real rmean_T(nr+1), rshl(nr+1), rmean_rho(nr+1), rshl_tomo(nr_tomo+1)
real temp, x, rad, fac, depth_crit, depth_cut, rmax_loop(nr+1), rmax(nr+1)
real rmean_T_loop(nr+1), border_pos, border_neg, T_mean, rmean_rho_loop(nr+1)
real, allocatable:: seis(:,:), shar_tomo(:,:,:), shar(:,:,:), T(:,:)
real T_pt(1,nn), T_bottom, var_scale(nr+1), var_scale_max, rmax_loop_proc(nr+1)
real depth(nr+1), meantemp(129), meantemp_bs(129), var_scale2(nr+1)
real max_model(129), max_model_orig(129)
real R_max, R_min, max_mod(nr+1), meantmp, quot

 character MP_TYP, comp*12, gpath*80, lpath*50, cname*3, cstage*4, cstage2*2, fname*40, char1*4
 character titl(4,4)*8, cstage3*3, char2*3

! ################################################
! ################# MAIN PART #######################
! ################################################

! ########################################
! Kommentar zur Konvertierung:
! In der Tomografie werden die Strukturen zu grob aufgelöst dargestellt, d.h. das Volumen überschätzt.
! In Folge dessen könnte die gesamte Dichteanomalie zu groß ausfallen. Dies könnte man korrigieren, in dem
! man die Dichteanomalien runterskaliert (ca. bis Faktor 1/2).
! Auf der anderen Seite sind Tomografien auf Grund der Inversionsprozesse gedämpft, was die Anomalie-Amplituden
! wiederum abschwächt.
! Hier muss man nun einen Kompromiss finden. Ein Vergleich zwischen Dichte-Anomalien aus der Tomografie und den
! Anomalien aus einer Konvektionsrechnung sollte einen ersten Anhaltspunkt geben.
! ########################################

common /mtemp/ meantemp, meantemp_bs, max_model, max_model_orig
 
! ######################    
! Initialize parallel communications.   

 CALL MPI_INIT(ierr)
 CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid2, ierr)
 CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc2, ierr)
nproc=(mt/nt)**2*10/nd
lvproc = 1.45*log(real(mt/nt))
mproc  = 2**(2*lvproc)

! ######################    
 
R_max=6.370e+06 
R_min=3.480e+06
T_bottom=4200.0 
allocate(T(nn*(nr+1),nproc/nproc2),seis(nn*(nr+1),nproc/nproc2))
         
! read in sh coefficients of GRAND s-wave
open(120,file='./Grand129.dat',status='unknown')
	read(120,*) lmax_tomo
	allocate(shar_tomo(2,(lmax_tomo+1)*(lmax_tomo+2)/2,nr_tomo+1))
	allocate(shar(2,(lmax_tomo+1)*(lmax_tomo+2)/2,nr+1))
	do i=1,2*129
		read(120,*) temp
	enddo
	read(120,*) shar_tomo
 close(120)

! skip layers of tomography model if nr<128
! interpolate linearly between tomography layers if nr>128
quot=real(nr_tomo)/real(nr)
do ir=1,nr+1
	ind(ir)=int(quot*(ir-1))+1
	if(quot>1) then
		ind_mod(ir)=0.0
	else
		ind_mod(ir)=mod(ir-1,int(1.0/quot))
	endif
	if(ir<nr+1) then
		shar(:,:,ir)=shar_tomo(:,:,ind(ir))+(shar_tomo(:,:,ind(ir)+1)-shar_tomo(:,:,ind(ir)))*ind_mod(ir)*quot
	else
		shar(:,:,ir)=shar_tomo(:,:,ind(ir))
	endif
	rshl(ir)=R_min+(R_max-R_min)*real(nr+1-ir)/real(nr)

	! interpolate in case of nr>128 (mt>256)
	! do not interpolate between last layer and CMB (val=0.0) to prevent unappropriate values
	if(ind(ir)<128) then
		max_mod(ir)=max_model(ind(ir))+(max_model(ind(ir)+1)-max_model(ind(ir)))*ind_mod(ir)*quot
	else
		max_mod(ir)=max_model(ind(ir))
	endif
enddo

rmean_T=0.0
rmean_rho=0.0
rmax=0.0

border_pos=1000000000.0
border_neg=-1000000000.0
depth_crit=-1.0
depth_cut=-1.0

do pro=1,nproc/nproc2
	
	myid=myid2*nproc/nproc2+pro-1
	
! generate s-wave velocity anomalies on TERRA grid from sh coefficients
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	call gridinit
	call shtofe(seis(:,pro),shar,lmax_tomo)
	
! convert anomalies to absolute values using PREM
k=1
do ir=1,nr+1
	! One additional kilometer to fit the range of the min. models
	depth(ir)=rshl(1)+1000.0-rshl(ir)
	rad=rshl(ir)/1000.0
	x=(rad+1.0)/6371.0
	do id=1,nd
		do ii=1,(nt+1)**2
			if(depth(ir)/1000.0<=200.0.and.seis(k,pro)>0.0) seis(k,pro)=0.0
			if((id==1.or.id==2.or.id==5.or.id==6.or.id==10)&
     &		   .and.depth(ir)/1000.0<=350.0.and.seis(k,pro)>0.0) seis(k,pro)=0.0
			!if(depth(ir)/1000.0<=200.0.and.seis(k,pro)<0.0) seis(k,pro)=seis(k,pro)/2.0
			
			fac=1.0+seis(k,pro)/100.0
			if(rad>=6151.0) then
				seis(k,pro)=(2.1519+2.3481*x)*fac
			else if(rad>=5971.0) then
				seis(k,pro)=(8.9496-4.4597*x)*fac
			else if(rad>=5771.0) then
				seis(k,pro)=(22.3512-18.5856*x)*fac
			else if(rad>=5701.0) then
				seis(k,pro)=(9.9839-4.9324*x)*fac
			else if(rad>=5600.0) then
				seis(k,pro)=(22.3459-17.2473*x-2.0834*x**2.0+0.9783*x**3.0)*fac
			else if(rad>=3630.0) then
				seis(k,pro)=(11.1671-13.7818*x+17.4575*x**2.0-9.2777*x**3.0)*fac
			else
				seis(k,pro)=(6.9254+1.4672*x-2.0834*x**2.0+0.9783*x**3.0)*fac
			endif
			k=k+1
		enddo
	enddo
enddo

! convert absolute velocities to absolute temperatures
 call read_MP_DATA

do ir=1,nr+1
	i = (ir-1)*nn+1
	j = ir*nn
	
	if(ir==1) then
		T(i:j,pro)=300.0
	else if(ir==nr+1) then
		T(i:j,pro)=T_bottom
	else
		if(ir==2) T_mean=850.0
		call convert_vpvsrho(depth(ir),1,seis(i:j,pro),nn,'s',T_pt,T_mean,ir)
		
		! At that point, temperatures could theoretically exceed boundary values (> 4200 K)
		T(i:j,pro)=T_pt(1,:)
	endif
enddo

! calculate the mean radial temperatures
 call rad_mean(T(:,pro),rmean_T_loop)
rmean_T=rmean_T+rmean_T_loop

enddo !pro


do pro=1,nproc/nproc2

! convert absolute temperatures into temperature anomalies

	!border_pos=0.0
	!border_neg=-0.0
	!depth_crit=300.0
	!depth_cut=150.0
	
	rmax_loop=0.0
	rmax_loop_proc=0.0
	
	k=1
	do ir=1,nr+1
		do id=1,nd
			do ii=1,(nt+1)**2
				T(k,pro)=(T(k,pro)-rmean_T(ir))/rmean_T(ir)
				if(depth(ir)<=depth_cut*1000.0) then
					T(k,pro)=0.0
				else if(depth(ir)<=depth_crit*1000.0) then
					if(T(k,pro)*100<border_neg) then
						T(k,pro)=border_neg/100.0
					elseif(T(k,pro)*100>border_pos) then
						T(k,pro)=border_pos/100.0
					endif
				endif
				
				! determine maximum positive variations in each layer
				if(T(k,pro)>rmax_loop_proc(ir)) rmax_loop_proc(ir)=T(k,pro)

				k=k+1
			enddo
		enddo
	enddo
	
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(rmax_loop_proc,rmax_loop,nr+1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(rmax_loop,(nr+1),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! determine maximal variations in each layer
! and respective scaling factor to reduce variations in case of a too high temperature
	do ir=1,nr+1
		! maximum over all loops
		if(rmax_loop(ir)>rmax(ir)) rmax(ir)=rmax_loop(ir)
		
		! calculate scaling factor 'var_scale' according to convection model
		if(max_mod(ir)==0.0) then
			var_scale(ir)=rmax(ir)
		else
			var_scale(ir)=rmax(ir)/max_mod(ir)*100.0
		endif
		! determine maximal scaling factor (currently not used during further proceedings)
		if(var_scale(ir)>var_scale_max) var_scale_max=var_scale(ir)
	enddo
	
enddo !pro

if(myid2==0) write(*,*) var_scale

rmean_T=0.0
rmax=0.0

do pro=1,nproc/nproc2

	myid=myid2*nproc/nproc2+pro-1
	
	k=1
	do ir=1,nr+1
		do id=1,nd
			do ii=1,(nt+1)**2
			
				! restrict the anomalies to values of the convection model
				if(var_scale(ir)/=0.0.and.var_scale(ir)>1.0) T(k,pro)=T(k,pro)/var_scale(ir)
				
				! use the (scaled) variations from the tomography model on a convection geotherm
				! interpolate geotherm in case of nr>128
				if(ir<nr+1) then
					meantmp=meantemp(ind(ir))+(meantemp(ind(ir)+1)-meantemp(ind(ir)))*ind_mod(ir)*quot
				else
					meantmp=meantemp(ind(ir))
				endif
				T(k,pro)=meantmp*(1.0+T(k,pro))		! absolute temperatures

				if(T(k,pro)>4200.0.and.myid2==0.and.ir<nr+1) write(*,*) T(k,pro)
				k=k+1
			enddo
		enddo
	enddo
	
	if(myid==0) write(*,*) "begin with output"
	write(char1,'(I4.4)') myid
	write(char2,'(I3.3)') mt
	open(unit=314,file='../TOMO1024_mt512/tomo.'//trim(char1)//'.00',status='replace')

	call vecout(T(:,pro),rshl,314,1)
	close(314)

	! calculate the mean radial temperatures
 	call rad_mean(T(:,pro),rmean_T_loop)
	rmean_T=rmean_T+rmean_T_loop
enddo !pro

if(myid2==0) write(*,*) rmean_T

! leave parallel communication before exiting
if(myid2==0) write(*,*) "Servus!"
 call MPI_FINALIZE(ierr)

end program SHconvert_Grand


! ########################################################
! ############# SUBROUTINES #################################
! ########################################################
  
subroutine rad_mean(s,rm)
implicit none

include 'size.h'
include 'pcom.h'

integer k, ir, id, ii
real s((nt+1)**2*nd*(nr+1))
real rmean(nr+1), rm(nr+1)

k=1
rmean=0.0
rm=0.0

do ir=1,nr+1	
	do id=1,nd
		do ii=1,(nt+1)**2
			rmean(ir)=rmean(ir)+s(k)
			k=k+1
		enddo
	enddo
enddo
	
rmean=rmean/nd/(nt+1)**2/nproc
	
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 call MPI_REDUCE(rmean,rm,nr+1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(rm,(nr+1),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)

end subroutine rad_mean


! #####################################################
! ############# SUBROUTINE vecout  #########################
! #####################################################
subroutine vecout(u,rshl,nf,sz)
implicit none
 
!...  This routine writes the nodal field u to logical unit nf
!...  using 1pe10.3 format when ifmt = 0 and f10.3 when ifmt = 1.
 
      include 'size.h'
      
      integer i, nf, sz
      real rshl(nr+1), propr(20), u((nt+1)**2*nd*(nr+1))
      character*8 titl(4,4)
      
      titl=""
 
      write(nf,10) nr, nt
 10   format(2i5)
 
      write(nf,20) titl
 20   format(4a8)
 
      write(nf,30) (rshl(i),i=1,nr+1)
      write(nf,30) propr
 30   format(1p10e15.8)

	write(nf,50) (u(i),i=1,(nt+1)**2*nd*(nr+1)*sz)
 50      format(15f10.3)
 
end subroutine vecout


! #####################################################
! ############# SUBROUTINE read_MP_DATA ###### ##############
! #####################################################
subroutine read_MP_DATA
	implicit none

	include 'min.h'
	
      if(MP_TYP=='A') then
! Antonio`s MP dataset is used
      	call read_MP_DATA_Antonio
      elseif (MP_TYP=='S') then
! Lars Stixrude`s MP dataset is used
      	call read_MP_DATA_Stixrude
      else
      	write(*,*) "Desired mineralogical model not available!"
      	return
      endif

end subroutine read_MP_DATA


! #####################################################
! ############# SUBROUTINE read_MP_DATA_Stixrude ##############
! #####################################################
subroutine read_MP_DATA_Stixrude
implicit none
	
	include 'min.h'
	
	integer i

	real ttt(nT_S), depth_S(nD_S), temp_S(nT_S), vs_S(nD_S, nT_S)
	real vsvprho_S(3,nD_S,nT_S), rho_S(nD_S, nT_S)
	real vp_S(nD_S, nT_S)

	character MP_dir*150
	
	common /mine_S/ depth_S, temp_S, vsvprho_S
		 
	if(comp=='PYROLITE') then
		MP_dir='../MP_DATA/MP_Stixrude/Pyrolite'
	elseif(comp=='PICLOGITE') then
		MP_dir='../MP_DATA/MP_Stixrude/Piclogite'
	else
		write(*,*) "Desired composition not available!"
		return
	endif

	open(803,file=trim(MP_dir)//'/MP_Stixrude_Pyrolite_depth.dat',action='read',status='old')
	open(804,file=trim(MP_dir)//'/MP_Stixrude_Pyrolite_T.dat',action='read',status='old')
	open(805,file=trim(MP_dir)//'/MP_Stixrude_Pyrolite_vs.dat',action='read',status='old')
	open(806,file=trim(MP_dir)//'/MP_Stixrude_Pyrolite_vp.dat',action='read',status='old')
	open(807,file=trim(MP_dir)//'/MP_Stixrude_Pyrolite_rho.dat',action='read',status='old')

	read(804,*) temp_S

	do i=1,nD_S
		read(803,*) depth_S(i)
		read(805,*) ttt
	  	vs_S(i,:)=ttt
		read(806,*) ttt
	 	vp_S(i,:)=ttt
		read(807,*) ttt
	  	rho_S(i,:)=ttt 
	enddo

	 close(803)
	 close(804)
	 close(805)
	 close(806)
	 close(807)

	! convert km in m and km/s in m/s
	depth_S=depth_S*1e3
	vsvprho_S(1,:,:)=vp_S*1e3
	vsvprho_S(2,:,:)=vs_S*1e3
	vsvprho_S(3,:,:)=rho_S*1e3

end subroutine read_MP_DATA_Stixrude


! #####################################################
! ############# SUBROUTINE read_MP_DATA_Antonio ###############
! #####################################################
subroutine read_MP_DATA_Antonio
implicit none

	include 'min.h'

      integer i,j

	real depth_A(nD_A), temp_A(nT_A), K_A(nD_A, nT_A)
	real vsvprho_A(3,nD_A,nT_A), rho_A(nD_A,nT_A)
	real xxx, G_A(nD_A, nT_A)
     
	common /mine_A/ depth_A, temp_A, vsvprho_A

	 character MP_dir*150

	if(comp=='PYROLITE') then
		MP_dir='../MP_DATA/MP_Antonio/Pyrolite'
	elseif(comp=='PICLOGITE') then
		MP_dir='../MP_DATA/MP_Antonio/Piclogite'
	else
		write(*,*) "Desired composition not available!"
		return
	endif

	open(804,file=trim(MP_dir)//'/MP_Antonio_Pyrolite_Depth.dat',status='old')
	do i=1,nD_A
	read(804,*) depth_A(i)
	enddo
	close(804)

	open(805,file=''//trim(MP_dir)//'/MP_Antonio_Pyrolite_T.dat',status='old')
	read(805,*) (temp_A(i),i=1,nT_A)
	close(805)

	open(806,file=trim(MP_dir)//'/MP_Antonio_Pyrolite_G.dat',status='old')
	do j=1,nD_A
	read(806,*) (G_A(j,i),i=1,nT_A)
	enddo
	close(806)

	open(807,file=trim(MP_dir)//'/MP_Antonio_Pyrolite_K.dat',status='old')
	do j=1,nD_A
	read(807,*) (K_A(j,i),i=1,nT_A)
	enddo
	close(807)

	open(808,file=trim(MP_dir)//'/MP_Antonio_Pyrolite_rho.dat',status='old')
	do j=1,nD_A
	read(808,*) (rho_A(j,i),i=1,nT_A)
	enddo
	 close(808)

	! convert km in m, GPa in Pa, 1e3 kg/m^3 in kg/m^3 and km/s in m/s
	depth_A=depth_A*1e3
	G_A=G_A*1e9
	K_A=K_A*1e9
	rho_A=rho_A/1e3

	do i=1,nD_A
	do j=1,nT_A

	xxx=G_A(i,j)/rho_A(i,j)
	if(xxx<0) xxx=0
	vsvprho_A(1,i,j)=sqrt((K_A(i,j) + 4./3.*G_A(i,j))/rho_A(i,j))
	vsvprho_A(2,i,j)=sqrt(xxx)
	vsvprho_A(3,i,j)=rho_A(i,j)

	enddo
	enddo

end subroutine read_MP_DATA_Antonio


! #####################################################
! ################## SUBROUTINE convert_T ###################
! #####################################################
subroutine convert_T(dpth,size_D,temp_cv,size_T,vp_new,vs_new,rho_new)
! This routine converts temperatures stored in array temp_cv at given depths dpth
! into vp, vs and density.
	
	include 'min.h'
	
	common /mine_A/ depth_A(nD_A), temp_A(nT_A), vsvprho_A(3,nD_A,nT_A)
	common /mine_S/ depth_S(nD_S), temp_S(nT_S), vsvprho_S(3,nD_S,nT_S)
	
	integer i, size_D, size_T
	
	real dpth(size_D), temp_cv(size_T), rho_new(size_D,size_T)
	real vp_new(size_D,size_T), vs_new(size_D,size_T)
	
	real, allocatable:: mat_tmp(:,:)
	real mat_interpol(size_D, size_T)
	
	do i=1,3
		if(MP_TYP=='A') then
			allocate(mat_tmp(nD_A,nT_A))
			mat_tmp=vsvprho_A(i,:,:)
			call interp_2D(depth_A, nD_A, temp_A, nT_A, mat_tmp,dpth, size_D, temp_cv, size_T, mat_interpol)
    		else if(MP_TYP=='S') then
    			allocate(mat_tmp(nD_S,nT_S))
    			mat_tmp=vsvprho_S(i,:,:)
			call interp_2D(depth_S, nD_S, temp_S, nT_S, mat_tmp, dpth, size_D, temp_cv, size_T, mat_interpol)
     		endif
     		
     		deallocate(mat_tmp)
		if(i==1) vp_new=mat_interpol/1000.0
		if(i==2) vs_new=mat_interpol/1000.0
		if(i==3) rho_new=mat_interpol/1000.0
	enddo

end subroutine convert_T


! #####################################################
! ################## SUBROUTINE convert_vpvsrho ###############
! #####################################################
subroutine convert_vpvsrho(dpth,size_D,vel_cv,size_vel,v_typ,T_new,T_mean,ir)
implicit none
! This routine converts vp, vs or rho stored in array vel_cv at given depths dpth
! into temperature.

	include 'min.h'
	include 'size.h'
	include 'pcom.h'
	
	common /mine_A/ depth_A, temp_A, vsvprho_A
	common /mine_S/ depth_S, temp_S, vsvprho_S
	
	integer i, j, k, l, size_D, size_vel, sz1, sz2, val, min_val, ir
	
	real depth_A(nD_A), temp_A(nT_A), vsvprho_A(3,nD_A,nT_A)
	real depth_S(nD_S), temp_S(nT_S), vsvprho_S(3,nD_S,nT_S)
	real dpth(size_D), vel_cv(size_vel), T_new(size_D,size_vel), frac, velo, min_dist
	real T_mean, min_dist_mean, border
	real, allocatable:: mat_tmp(:,:), dpth_tmp(:), T_tmp(:), vel_compl(:), T_pos(:), dist_mean(:)
	real, allocatable:: T_pos_dist(:)

	character v_typ, bla
	
	if(MP_TYP=='A') then
		sz1=nD_A
		sz2=nT_A
		allocate(mat_tmp(sz1,sz2),dpth_tmp(sz1),T_tmp(sz2))
		dpth_tmp=depth_A
		T_tmp=temp_A
		if(v_typ=='p') then
			mat_tmp=vsvprho_A(1,:,:)
		else if(v_typ=='s') then
			mat_tmp=vsvprho_A(2,:,:)
		else if(v_typ=='r') then
			mat_tmp=vsvprho_A(3,:,:)
		endif
    	else if(MP_TYP=='S') then
    		sz1=nD_S
		sz2=nT_S
    		allocate(mat_tmp(sz1,sz2),dpth_tmp(sz1),T_tmp(sz2))
    		dpth_tmp=depth_S
    		T_tmp=temp_S
    		if(v_typ=='p') then
			mat_tmp=vsvprho_S(1,:,:)
		else if(v_typ=='s') then
			mat_tmp=vsvprho_S(2,:,:)
		else if(v_typ=='r') then
			mat_tmp=vsvprho_S(3,:,:)
		endif
     	endif

	border=7.46
	T_new=0.0
     	allocate(vel_compl(sz2),T_pos(sz2),dist_mean(sz2),T_pos_dist(sz2))
     	do i=1,size_D
     		call hunt(dpth_tmp,sz1,dpth(i),val)
     		frac=(dpth(i)-dpth_tmp(val))/(dpth_tmp(val+1)-dpth_tmp(val))
     		do j=1,sz2
     			! 1D interpolation
     			vel_compl(j)=(mat_tmp(val,j)+(mat_tmp(val+1,j)-mat_tmp(val,j))*frac)/1000.0
     		enddo
     		do j=1,size_vel
     			velo=vel_cv(j)
     			l=1
     			min_dist_mean=5000.0
     			if(velo>border) velo=border
     			do k=1,sz2-1
     				if((velo>=vel_compl(k).and.velo<=vel_compl(k+1)).or.&
     				  &(velo<=vel_compl(k).and.velo>=vel_compl(k+1))) then
     					! 1D interpolation
     					! All possible T-values are collected in T_pos since solution is not unique.
     					T_pos(l)=T_tmp(k)+(T_tmp(k+1)-T_tmp(k))*(velo-vel_compl(k))/(vel_compl(k+1)-vel_compl(k))
     					dist_mean(l)=abs(T_pos(l)-T_mean)
     					T_pos_dist(l)=abs(vel_compl(k)-vel_compl(k+1))
     					
     					if(dist_mean(l)<min_dist_mean) min_dist_mean=dist_mean(l)
     					!if(ir==nr+1) write(*,*) velo, vel_compl(k), vel_compl(k+1)
     						
     					l=l+1
     				endif 
     			enddo
     			min_dist=10000.0
     			min_val=1
     			if(l>1) then
     				do k=1,l-1
					!if(T_pos_dist(k)<min_dist.and.abs(dist_mean(k)-min_dist_mean)<300.0) then
					if(dist_mean(k)<min_dist) then
     						!min_dist=T_pos_dist(k)
     						min_dist=dist_mean(k)
     						min_val=k
     					endif
    					!if(ir==nr.and.T_pos(k)>4200.0) then
     						!write(*,*) T_pos(k), T_pos_dist(k), T_mean
     						!write(*,*) dist_mean(k), min_dist_mean, abs(dist_mean(k)-min_dist_mean)
     						!read(*,*) bla
     					!endif     				
     				enddo
   				T_new(i,j)=T_pos(min_val)
   				!if(T_new(i,j)>4200.0) then
   				!	T_new(i,j)=4200.0
   				!else if(T_new(i,j)<300.0) then
   				!	T_new(i,j)=300.0
   				!endif
   			else if(velo<5.0.and.v_typ=='s') then
   				T_new(i,j)=300.0
   			else if(v_typ=='s') then
   				T_new(i,j)=4200.0
   			endif
	
     			!if(T_new(i,j)<940.0.and.ir>100) then
     			!	write(*,*) "solutions: ", l-1, min_val
     			!	write(*,*) velo, T_new(i,j)
     			!	write(*,*) T_new(i,j), myid
     			!	read(*,*) bla
     			!endif

     			T_mean=(T_mean*(j-1)+T_new(i,j))/j

     		enddo
     	enddo
	deallocate(mat_tmp,dpth_tmp,T_tmp,vel_compl,T_pos)
	     		
end subroutine convert_vpvsrho


! #####################################################
! ################### SUBROUTINE interp_2D ##################
! #####################################################
subroutine interp_2D(x_cur,nx,y_cur,ny,mat_cur, x_vec,nx_v,y_vec,ny_v,mat_new)
implicit none	
	
	integer i,j
	integer nx, ny, nx_v, ny_v
	integer x1, x2, y1, y2
	
	real x_cur(nx), y_cur(ny), mat_cur(nx,ny)
	real x_vec(nx_v), y_vec(ny_v), mat_new(nx_v,ny_v)
	real f_x, f_y, frac
	
	do i=1,nx_v
		call hunt(x_cur,nx,x_vec(i),x1)
		if(x1==0.or.x1==nx) write(*,*) x_vec(i), x1, "X-PROBLEM!!!!!!!!"
		x2=x1+1
		frac=(x_vec(i)-x_cur(x1))/(x_cur(x2)-x_cur(x1))
		do j=1,ny_v
			call hunt(y_cur,ny,y_vec(j),y1)
			if(y1==0.or.y1==ny) write(*,*) y_vec(j), y1, "Y-PROBLEM!!!!!!!!"
			y2=y1+1
			
			f_x=mat_cur(x1,y1)+(mat_cur(x2,y1)-mat_cur(x1,y1))*frac
    			f_y=mat_cur(x1,y2)+(mat_cur(x2,y2)-mat_cur(x1,y2))*frac
  
     			mat_new(i,j)=f_x+(f_y-f_x)*(y_vec(j)-y_cur(y1))/(y_cur(y2)-y_cur(y1))
		enddo
	enddo

end subroutine interp_2D


! #####################################################
! ################### SUBROUTINE hunt #####################
! #####################################################
subroutine hunt(x_vec,n,x,jlo)
implicit none
	
	integer n, jlo, jhi, jnew
	real x_vec(n), x
	
	jlo=1
	jhi=n
	if(x>x_vec(n)) then
		jlo=n
	else if(x<x_vec(1)) then
		jlo=0
	else
		do while(jhi-jlo>1)
			jnew=floor(real(jlo+jhi)/2.0)
			if(x<=x_vec(jnew)) then
				jhi=jnew
			else
				jlo=jnew
			endif
		enddo
	endif

end subroutine hunt


subroutine fetosh(s_tmp,shar,lmax)
!implicit none

! This sub-routine computes the coefficients of the spherical harmonics
! expansion of a given 3D grid function. The expansion is performed for each
! layer of the radial discretisation.
! \param lun  Logical Unit Number for the output file to which the computed
!             coefficients will be appended. IO is only performed by the MPI
!             process with rank 0.
! \param s    Array containing the values of the scalar function for which the
!             expansion is to be computed
! \param shar used for storing the computed coefficients (work array, or
!             INTENTOUT?)
 
!...  This routine generates spherical harmonic coefficients
		
      include 'size.h'
      include 'pcom.h'
      integer lmax, ir, id, ii, k, l, m
      real s_tmp((nt+1)**2*nd*(nr+1)), smax, r4pi, a0, phi
      real s((nt+1)**2,nd,nr+1), shar(2,(lmax+1)*(lmax+2)/2,nr+1)
      real plm((lmax+1)*(lmax+2)), csm(0:128,2)
      real shar_tmp((lmax+1)*(lmax+2)*(nr+1))
      real shar_tmp2((lmax+1)*(lmax+2)*(nr+1))
      
      common /mesh/ xn((nt+1)**2,nd,3)
      common /ndar/  arn((nt+1)**2),  arne((nt+1)**2),&
     &              rarn((nt+1)**2), rarne((nt+1)**2)
 
      smax = 0.
      r4pi = 0.125/asin(1.)
 
      shar=0.0
      
      k=0
	do ir=1,nr+1
		do id=1,nd
			do ii=1,(nt+1)**2
				k=k+1
				s(ii,id,ir)=s_tmp(k)
			enddo
		enddo
	enddo

!     This 'hack' fixes a problem in the use of arn in the loop below.
!     The loop runs over all nodes in the local sub-domain. However, Terra
!     employs an element-oriented domain decomposition. Thus, nodes along the
!     subdomain edges belong to several subdomains. arn stores zeros along the
!     upper right and lower right edge. In this fashion we can avoid multiple
!     contributions from different subdomains from the loop below, when doing
!     psumlong in the end of this routine. However, those process(es) to whose
!     subdomain(s) the north and south poles belong stores the area associated
!     with those nodes. Thus, we get five contributions from these nodes.
!     We correct this by scaling the corresponding area value in arn.
!     The hack must be undone again before leaving the routine. Otherwise
!     inconsistencies might occur in other parts of the code.
if(myid==0.OR.(nd .EQ. 5 .AND. myid .EQ. nproc/2)) THEN
         arn(1) = arn(1) / 5.
endif
     
do ii=1,(nt+1)**2
         a0 = r4pi*arn(ii)
         do id=1,nd
            phi = atan2(xn(ii,id,2) + 1.e-30, xn(ii,id,1))
            do m=0,lmax
               csm(m,1) = cos(m*phi)
               csm(m,2) = sin(m*phi)
            enddo
 
            if(mod(id, 5) .eq. 1) call plmbar(plm,plm,lmax,xn(ii,id,3),0)
		k=0
            do l=0,lmax
               do m=0,l
                  k=k+1
                  do ir=1,nr+1
                     aa = a0*plm(k)*s(ii,id,ir)
                     shar(1,k,ir) = shar(1,k,ir) + aa*csm(m,1)
                     shar(2,k,ir) = shar(2,k,ir) + aa*csm(m,2)
                  enddo
                enddo
            enddo
	enddo
enddo

! Undo 'hack' by re-setting arn to original value
if(myid .EQ. 0 .OR. (nd .EQ. 5 .AND. myid .EQ. nproc/2)) THEN
         arn(1) = arn(1) * 5.
endif

k=0
do ir=1,nr+1
	do i=1,2
		kk=0
		do l=0,lmax
			do m=0,l
				k=k+1
				kk=kk+1
				shar_tmp(k)=shar(i,kk,ir)
			enddo
		enddo
	enddo
enddo

 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 call MPI_REDUCE(shar_tmp,shar_tmp2,(lmax+1)*(lmax+2)*(nr+1),MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 
if(myid2==0) then
   k=0
	do ir=1,nr+1
		do i=1,2
			kk=0
			do l=0,lmax
				do m=0,l
					k=k+1
					kk=kk+1
					shar(i,kk,ir)=shar_tmp2(k)
				enddo
			enddo
		enddo
	enddo
endif

end subroutine


subroutine pbar(c,l,m,p)
 
!...  This routine calculates the value of the normalized associated
!...  Legendre function of the first kind of degree l and order m
!...  for the real argument c, for 0 .le. m .le. l.
 
      sqrt2 = 1.414213562373092
 
      if(m .ne. 0) then
         p = sqrt2
         s = sqrt(1. - c**2)
         do i=1,m
            p = sqrt(real(2*i + 1)/real(2*i))*s*p
         end do
      else
         p = 1.
      endif
 
      if(l .eq. m) return
 
      p1 = sqrt2
      do j=m+1,l
         p2 = p1
         p1 = p
         p  = 2.*sqrt((real(j**2) - 0.25)/real(j**2 - m**2))*c*p1&
     &         - sqrt (real((2*j + 1)*(j - m - 1)*(j + m - 1))&
     &                /real((2*j - 3)*(j - m)*(j + m)))*p2
      end do
 end


subroutine plmbar(p,dp,lmax,z,ideriv)
!    Evaluates normalized associated Legendre function P(l,m) as
!    function of z=cos(colatitude); also derivative dP/d(colatitude).
!    Uses recurrence relation starting with P(l,l) and then increasing
!    l keeping m fixed.  Normalization is:
!                  Integral(Y(l,m)*Y(l,m)) = 4.*pi,
!                  where Y(l,m) = P(l,m)*exp(i*m*longitude),
!    which is incorporated into the recurrence relation. p(k) contains
!    p(l,m) with k=(l+1)*l/2+m+1; i.e. m increments through range 0 to
!    l before incrementing l. Routine is stable in single and double
!    precision to l,m = 511 at least; timing is proportional to lmax**2
!    R.J.O'Connell 7 Sept. 1989; added dp(z) 10 Jan. 1990
!
!    Added precalculation and storage of square roots srl(k) 31 Dec 1992

	integer, parameter:: lmaxx=100
      dimension p(*),dp(*)
!     --dimensions must be p((lmax+1)*(lmax+2)/2) in calling program
 
      common /plm0/   f1((lmaxx+1)*(lmaxx+2)/2),&
     &                f2((lmaxx+1)*(lmaxx+2)/2),&
     &              fac1((lmaxx+1)*(lmaxx+2)/2),&
     &              fac2((lmaxx+1)*(lmaxx+2)/2), srt(2*lmaxx+2)
      data ifirst /1/
      save ifirst
 
      if (lmax.lt.0.or.abs(z).gt.1.) stop 'bad arguments'
!     --set up sqrt and factors on first pass
      if(ifirst.eq.1) then
          ifirst = 0
          do k=1,2*lmax+2
            srt(k) = sqrt(real(k))
         end do
 
         if (lmax .eq. 0) then
            p(1) = 1.0
            if(ideriv .ne. 0) dp(1) = 0.
            return
         end if

!        --case for m > 0
          kstart = 1
         do m=1,lmax
!        --case for P(m,m)
             kstart = kstart + m + 1
            if(m .ne. lmax) then

!              --case for P(m+1,m) 
               k = kstart + m + 1
 
!              --case for P(l,m) with l > m+1
               if(m .lt. lmax-1) then
                  do l=m+2,lmax
                     k = k + l
                     f1(k) =  srt(2*l+1)*srt(2*l-1)/(srt(l+m)*srt(l-m))
                     f2(k) = (srt(2*l+1)*srt(l-m-1)*srt(l+m-1))&
     &                      /(srt(2*l-3)*srt(l+m)*srt(l-m))
                  end do
               end if
            end if
         end do
 
         k = 3
 
         do l=2,lmax
            k = k + 1
            do m=1,l-1
               k = k + 1
               fac1(k) = srt(l-m)*srt(l+m+1)
               fac2(k) = srt(l+m)*srt(l-m+1)
               if(m .eq. 1) fac2(k) = fac2(k)*srt(2)
            end do
            k = k + 1
         end do
 
      end if
 
!     --start calculation of Plm, etc.
 
!     --case for P(l,0)
 
      pm2   = 1.
      p(1)  = 1.
      if(ideriv .ne. 0) dp(1) = 0.
 
      if(lmax .eq. 0) return
 
      pm1   = z
      p(2)  = srt(3)*pm1
      k     = 2
 
      do l=2,lmax
         k = k + l
         plm  = (real(2*l-1)*z*pm1 - real(l-1)*pm2)/real(l)
         p(k) =   srt(2*l+1)*plm
         pm2  =   pm1
         pm1  =   plm
      end do
 
!    --case for m > 0
 
      pmm    =  1.
      sintsq = (1.-z)*(1.+z)
      fnum   = -1.
      fden   =  0.
      kstart =  1
 
      do m=1,lmax
 
!        --case for P(m,m)
 
         kstart = kstart + m + 1
         fnum   = fnum + 2.
         fden   = fden + 2.0
         pmm    = pmm*sintsq*fnum/fden
         pm2    = sqrt(real(4*m+2)*pmm)
         p(kstart) = pm2
 
         if(m .ne. lmax) then
 
!           --case for P(m+1,m)
 
            pm1  = z*srt(2*m+3)*pm2
            k    = kstart + m + 1
            p(k) = pm1
 
!           --case for P(l,m) with l > m+1
 
            if(m .lt. lmax-1) then
 
               do l=m+2,lmax
                  k    = k + l
!                 f1   =  srt(2*l+1)*srt(2*l-1)/(srt(l+m)*srt(l-m))
!                 f2   = (srt(2*l+1)*srt(l-m-1)*srt(l+m-1))
!     &                 /(srt(2*l-3)*srt(l+m)*srt(l-m))
                  plm  = z*f1(k)*pm1 - f2(k)*pm2
                  p(k) = plm
                  pm2  = pm1
                  pm1  = plm
               end do
            endif
         endif
      end do
 
      if(ideriv .eq. 0) return
 
!     ---derivatives of P(z) wrt theta, where z=cos(theta)
      dp(2) = -p(3)
      dp(3) =  p(2)
      k     =  3
 
      do l=2,lmax
         k = k + 1

!        --treat m=0 and m=l separately
         dp(k)   = -srt(l)*srt(l+1)/srt(2)*p(k+1)
         dp(k+l) =  srt(l)/srt(2)*p(k+l-1)
 
            do m=1,l-1
               k     = k + 1
               dp(k) = 0.5*(fac2(k)*p(k-1) - fac1(k)*p(k+1))
            enddo
         k = k + 1
      end do
end
 
subroutine gridinit
!...  This routine initializes all arrays in common blocks /radl/,
!...  /mesh/, /ndar/, /volm/, and /xnex/.
 
      include 'size.h'
      include 'pcom.h'
      parameter (nxm=4000+(nt+1)**2*41)
      parameter (nopr=(nt/2+1)**2*ndo*189*(nr/2+1)*7/5+8000)
      common /fopr/ a(nopr+nv*ndo/nd*9+nv*ndo/nd*18*5/4),&
     &              mopr(0:10), mb(0:10)
      common /grid/ mxm(0:10), xm(nxm)
      common /mesh/ xn((nt+1)**2,nd,3)
      common /ndar/  arn((nt+1)**2),  arne((nt+1)**2),&
     &              rarn((nt+1)**2), rarne((nt+1)**2)
      common /radl/ rshl(nr+1), ird, ibc
      common /volm/ vol((nt+1)**2*(nr+1),2)
      common /xnex/ xne((nt+3)**2*nd*3)
 
!...  Generate nodal coordinate array xn for finest level.
 
      mm = (mt + 1)**2
      m0 = mm*30 + 1
 
      call grdgen(a,mt)
 
      if(nproc==1) then
         call scopy((mt+1)**2*30, a, 1, xn, 1)
      else
         call subarray(a,xn,0,10,nd,mt,nt,3)
      endif
 
!...  Generate the nodal area arrays.
 
      call ndarea(a(m0),a(m0+mm),a(m0+2*mm),a(m0+3*mm),a(m0+4*mm),a,mt)
 
      if(nproc==1) then
         call scopy((mt+1)**2, a(m0), 1, arn, 1)
         call scopy((mt+1)**2, a(m0+mm), 1, arne, 1)
         call scopy((mt+1)**2, a(m0+2*mm), 1, rarn, 1)
         call scopy((mt+1)**2, a(m0+3*mm), 1, rarne, 1)
      else
         call subarray(a(m0),     arn,  1,1,1,mt,nt,1)
         call subarray(a(m0+mm),  arne, 0,1,1,mt,nt,1)
         call subarray(a(m0+2*mm),rarn, 1,1,1,mt,nt,1)
         call subarray(a(m0+3*mm),rarne,0,1,1,mt,nt,1)
      endif
 
 end

subroutine subarray(a,ap,ibctype,kd,kdp,kt,ktp,nc)
!...  This routine reads the appropriate portion of the global array
!...  a into the sub-array ap needed for process myid.  Array a has
!...  dimensions (kt+1,kt+1,kd,nc) and ibctype specifies the type of
!...  boundary conditions to be applied along subdomain edges.  For
!...  ibctype = 0, none of the edge elements are set to zero; for
!...  ibctype = 1, all of the edge elements are set to zero; and for
!...  ibctype = 2, array ap is treated as a 7-point stencil and edge
!...  node components outside the subdomain are set to zero.
 
      include 'size.h'
      include 'pcom.h'
      real a(kt+1,kt+1,kd,nc), ap(ktp+1,ktp+1,kdp,nc)
 
!...  Determine subdomain limits.
 
      npedg = 2**lvproc
      iproc = mod(myid, mproc)
      rnm   = 1./real(npedg)
      j0    = iproc*rnm
      i0    = (real(iproc)*rnm - j0)/rnm
 
      ibeg  = i0*ktp + 1
      jbeg  = j0*ktp + 1
      iend  = ibeg + ktp
      jend  = jbeg + ktp
 
!...  Load sub-array.
 
      do ic=1,nc
         do id=1,kdp
            jd = id
            if(kdp.eq.5 .and. myid.ge.mproc) jd = id + 5
            j = 0
            do jj=jbeg,jend
               j = j + 1
               i = 0
               do ii=ibeg,iend
                  i = i + 1
                  ap(i,j,id,ic) = a(ii,jj,jd,ic)
               enddo
            enddo
         enddo
      enddo
 
!...  Apply boundary conditions.
 
      if(ibctype .eq. 0) return
 
      if(ibctype .eq. 1) then
         k  = ktp + 1
         i1 = 1
         if(iproc .eq. 0) i1 = 2
 
         do ic=1,nc
            do id=1,kdp
               do i=i1,ktp+1
                  ap(1,i,id,ic) = 0.
                  ap(i,k,id,ic) = 0.
               end do
            end do
         end do
      elseif(ibctype .eq. 2) then
         call subarraybc(ap,nc/7,ktp)
      endif
end subroutine


subroutine subarraybc(ap,nc,kt)
!...  This routine sets to zero the appropriate stencil components
!...  along the sub-domain boundaries for the operator ap.  The flag
!...  idiamond, when set to one, causes redundant components on the
!...  diamond edges to be set to zero.
 
      include 'size.h'
      include 'pcom.h'
      real ap(kt+1,kt+1,7,nc)
 
      k     = kt + 1
      npedg = 2**lvproc
      iproc = mod(myid, mproc)
      i0    = 1
      if(iproc .eq. 0) i0 = 2
 
      do j=1,nc
         do i=i0,kt+1
 
!...        Treat upper right edge.
            ap(1,i,4,j) = 0.
            ap(1,i,5,j) = 0.
 
            if(mod(iproc, npedg) .ne. 0) then
               ap(1,i,1,j) = 0.
               ap(1,i,3,j) = 0.
               ap(1,i,6,j) = 0.
            endif
 
!...        Treat upper left edge.
            ap(i,1,6,j) = 0.
            ap(i,1,7,j) = 0.
         enddo
 
         do i=1,kt+1
 
!...        Treat lower left edge.
            ap(k,i,2,j) = 0.
            ap(k,i,7,j) = 0.
 
!...        Treat lower right edge.
            ap(i,k,3,j) = 0.
            ap(i,k,4,j) = 0.
 
            if(iproc .le. mproc-1-npedg) then
               ap(i,k,1,j) = 0.
               ap(i,k,2,j) = 0.
               ap(i,k,5,j) = 0.
            endif
 
         enddo
      enddo
 end
 
 
subroutine scopy(nn,vin,iin,vout,iout)
 
!     This routine copies array vin into array vout.
      real vout(*), vin(*)
 
      do 10 ii=1,nn-2,3
      vout(ii)   = vin(ii)
      vout(ii+1) = vin(ii+1)
      vout(ii+2) = vin(ii+2)
 10   continue
 
      do 20 ii=3*(nn/3)+1,nn
      vout(ii)   = vin(ii)
 20   continue
end


subroutine grdgen(xn,nt)
!...  This routine generates the nodal coordinates xn for an
!...  icosahedral grid on the unit sphere.  The grid resolution
!...  corresponds to a subdivision of the edges of the original
!...  icosahedral triangles into nt equal parts.
 	
      real xn(nt+1,nt+1,10,3)
 
      fifthpi = 0.4*asin(1.)
      w       = 2.0*acos(1./(2.*sin(fifthpi)))
      cosw    = cos(w)
      sinw    = sin(w)
      lvt     = 1.45*log(real(nt))
      nn      = (nt+1)**2*10
 
      do id=1,10
 
         sgn = 1.
         if(id .ge. 6) sgn = -1.
         phi = (2*mod(id, 5) - 3 + (id - 1)/5)*fifthpi
 
         xn(   1,   1,id,1) =  0.
         xn(   1,   1,id,2) =  0.
         xn(   1,   1,id,3) =  sgn
         xn(nt+1,   1,id,1) =  sinw*cos(phi)
         xn(nt+1,   1,id,2) =  sinw*sin(phi)
         xn(nt+1,   1,id,3) =  cosw*sgn
         xn(   1,nt+1,id,1) =  sinw*cos(phi + fifthpi + fifthpi)
         xn(   1,nt+1,id,2) =  sinw*sin(phi + fifthpi + fifthpi)
         xn(   1,nt+1,id,3) =  cosw*sgn
         xn(nt+1,nt+1,id,1) =  sinw*cos(phi + fifthpi)
         xn(nt+1,nt+1,id,2) =  sinw*sin(phi + fifthpi)
         xn(nt+1,nt+1,id,3) = -cosw*sgn
 
         do k=0,lvt-1
 
            m  = 2**k
            l  = nt/m
            l2 = l/2
 
!...        rows of diamond--
 
            do j1=1,m+1
               do j2=1,m
                     i1 = (j1-1)*l + 1
                     i2 = (j2-1)*l + l2 + 1
                     call midpt(xn(i1,i2,id,1),xn(i1,i2-l2,id,1),&
     &                          xn(i1,i2+l2,id,1),nn)
               end do
            end do
 
!...        columns of diamond--
 
            do j1=1,m+1
               do j2=1,m
                     i1 = (j2-1)*l + l2 + 1
                     i2 = (j1-1)*l + 1
                     call midpt(xn(i1,i2,id,1),xn(i1-l2,i2,id,1),&
     &                          xn(i1+l2,i2,id,1),nn)
               end do
            end do
 
!...        diagonals of diamond--
 
            do j1=1,m
               do j2=1,m
                     i1 = (j1-1)*l + l2 + 1
                     i2 = (j2-1)*l + l2 + 1
                     call midpt(xn(i1,i2,id,1),xn(i1-l2,i2+l2,id,1),&
     &                          xn(i1+l2,i2-l2,id,1),nn)
               end do
            enddo
         enddo
      end do
 end


subroutine midpt(x,x1,x2,nn)
!...  This routine finds the midpoint x along the shorter great circle
!...  arc between points x1 and x2 on the unit sphere.
 
      real x(nn,3), x1(nn,3), x2(nn,3)
 
      do j=1,3
         x(1,j) = x1(1,j) + x2(1,j)
      enddo
 
      xnorm = 1./sqrt(x(1,1)**2 + x(1,2)**2 + x(1,3)**2)
 
      do j=1,3
         x(1,j) = xnorm*x(1,j)
      enddo
end


subroutine ndarea(arn,arne,rarn,rarne,area,xn,nt)
!...  This routine computes the areas, arn, associated with the nodes
!...  on the unit sphere as well as the reciprocal areas, rarn.  These
!...  arrays have zero values along the upper right and lower right
!...  diamond edges.  Arrays arne and rarne are identical except they
!...  contain the actual areas and reciprocal areas, respectively, in
!...  these edge locations.
 
      real  arn(nt+1,nt+1),  rarn(nt+1,nt+1), area(nt+1,nt+1,2)
      real arne(nt+1,nt+1), rarne(nt+1,nt+1), xn(*)
 
      call areacalc(area,xn,nt)
      arn=0.0
      rarn=0.0

!...  Treat interior nodes.
      do i2=2,nt
         do i1=2,nt
            arn(  i1,i2) = (area(i1  ,i2  ,1) + area(i1  ,i2  ,2)&
     &                   +  area(i1+1,i2  ,1) + area(i1  ,i2-1,2)&
     &                   +  area(i1+1,i2-1,1) + area(i1+1,i2-1,2))/3.
            rarn( i1,i2) =  1./arn(i1,i2)
            arne( i1,i2) =     arn(i1,i2)
            rarne(i1,i2) =    rarn(i1,i2)
         end do
      end do
 
!...  Treat edge nodes.
      do i=2,nt
         arn(     i,1) = (area(   i,1  ,1) + area(i   ,1,2)&
     &                 +  area( i+1,1  ,1))/1.5
         arn(  nt+1,i) = (area(nt+1,i  ,1) + area(nt+1,i,2)&
     &                 +  area(nt+1,i-1,2))/1.5
         rarn(    i,1) =  1./arn(   i,1)
         rarn( nt+1,i) =  1./arn(nt+1,i)
         arne(    i,1) =     arn(   i,1)
         arne( nt+1,i) =     arn(nt+1,i)
         rarne(   i,1) =    rarn(   i,1)
         rarne(nt+1,i) =    rarn(nt+1,i)
         arne(    1,i) =     arn(   i,1)
         arne( i,nt+1) =     arn(nt+1,i)
         rarne(   1,i) =    rarn(   i,1)
         rarne(i,nt+1) =    rarn(nt+1,i)
      end do
 
!...  Treat pentagonal nodes.
      arn(     1,   1) = 5.*area(2,1,1)/3.
      arn(  nt+1,   1) =    arn(1,1)
      rarn(    1,   1) = 1./arn(1,1)
      rarn( nt+1,   1) =   rarn(1,1)
      arne(    1,   1) =    arn(1,1)
      arne( nt+1,   1) =    arn(1,1)
      arne(    1,nt+1) =    arn(1,1)
      arne( nt+1,nt+1) =    arn(1,1)
      rarne(   1,   1) =   rarn(1,1)
      rarne(nt+1,   1) =   rarn(1,1)
      rarne(   1,nt+1) =   rarn(1,1)
      rarne(nt+1,nt+1) =   rarn(1,1)
 end subroutine
 
subroutine areacalc(area,xn,nt)
!...  This routine computes the areas of the spherical triangles in
!...  a diamond on the unit sphere.
 
      real area(nt+1,nt+1,2), xn(nt+1,nt+1,10,3), xv(3,3)

      area=0.0
       
      do i2=1,nt
         do i1=2,nt+1
            xv(1,1) = xn(i1  ,i2  ,1,1)
            xv(1,2) = xn(i1  ,i2  ,1,2)
            xv(1,3) = xn(i1  ,i2  ,1,3)
            xv(2,1) = xn(i1-1,i2+1,1,1)
            xv(2,2) = xn(i1-1,i2+1,1,2)
            xv(2,3) = xn(i1-1,i2+1,1,3)
            xv(3,1) = xn(i1-1,i2  ,1,1)
            xv(3,2) = xn(i1-1,i2  ,1,2)
            xv(3,3) = xn(i1-1,i2  ,1,3)
            t1 = xv(1,1)*xv(2,1) + xv(1,2)*xv(2,2) + xv(1,3)*xv(2,3)
            t2 = xv(2,1)*xv(3,1) + xv(2,2)*xv(3,2) + xv(2,3)*xv(3,3)
            t3 = xv(3,1)*xv(1,1) + xv(3,2)*xv(1,2) + xv(3,3)*xv(1,3)
            t1 = 0.5*acos(t1)
            t2 = 0.5*acos(t2)
            t3 = 0.5*acos(t3)
            s  = 0.5*(t1 + t2 + t3)
            a  = tan(s)*tan(s-t1)*tan(s-t2)*tan(s-t3)
            area(i1,i2,1) = 4.*atan(sqrt(a))

            xv(1,1) = xn(i1  ,i2  ,1,1)
            xv(1,2) = xn(i1  ,i2  ,1,2)
            xv(1,3) = xn(i1  ,i2  ,1,3)
            xv(2,1) = xn(i1  ,i2+1,1,1)
            xv(2,2) = xn(i1  ,i2+1,1,2)
            xv(2,3) = xn(i1  ,i2+1,1,3)
            xv(3,1) = xn(i1-1,i2+1,1,1)
            xv(3,2) = xn(i1-1,i2+1,1,2)
            xv(3,3) = xn(i1-1,i2+1,1,3)
            t1 = xv(1,1)*xv(2,1) + xv(1,2)*xv(2,2) + xv(1,3)*xv(2,3)
            t2 = xv(2,1)*xv(3,1) + xv(2,2)*xv(3,2) + xv(2,3)*xv(3,3)
            t3 = xv(3,1)*xv(1,1) + xv(3,2)*xv(1,2) + xv(3,3)*xv(1,3)
            t1 = 0.5*acos(t1)
            t2 = 0.5*acos(t2)
            t3 = 0.5*acos(t3)
            s  = 0.5*(t1 + t2 + t3)
            a  = tan(s)*tan(s-t1)*tan(s-t2)*tan(s-t3)
            area(i1,i2,2) = 4.*atan(sqrt(a))

         enddo
      enddo
end subroutine

!> Generate a grid function from spherical harmonics coefficients
!>
!> The sub-routine generates a 3D grid function from the given coefficients of
!> the expansion of that function in terms of spherical harmonics on each
!> radial layer.
subroutine shtofe(s,shar,lmax)
 
      include 'size.h'
      include 'pcom.h'
      integer lmax, lmaximum
      real s((nt+1)**2,nd,nr+1), shar(2,(lmax+1)*(lmax+2)/2,nr+1)
      real plm((lmax+1)*(lmax+2)), csm(0:128,2)
      real, parameter:: PI = 3.141592653589793
      
      common /mesh/ xn((nt+1)**2,nd,3)
	
	lmaximum=lmax 
      s=0.0
 
      do ii=1,(nt+1)**2
         do id=1,nd
 
            phi = atan2(xn(ii,id,2) + 1.e-30, xn(ii,id,1))
 
 		! ROTATION ABOUT 90 DEGREES
 		if(phi<3.0*PI/2.0) then
 			phi=phi+PI/2.0
 		else
 			phi=phi-3.0*PI/2.0
 		endif
 		
            do m=0,lmaximum
               csm(m,1) = cos(m*phi)
               csm(m,2) = sin(m*phi)
            end do
 
            if(mod(id, 5) .eq. 1)&
     &         call plmbar(plm,plm,lmaximum,xn(ii,id,3),0)
 
            k = 0
 
            do l=0,lmaximum
               do m=0,l
                  k = k + 1
 
                  do ir=1,nr+1
                     s(ii,id,ir) = ((s(ii,id,ir)&
     &                           +  (plm(k)*csm(m,1))*shar(1,k,ir))&
     &                           +  (plm(k)*csm(m,2))*shar(2,k,ir))
                  end do
 
               end do
            end do
 
         end do
	enddo
end
      

blockdata bdconvect

	real meantemp(129), meantemp_bs(129), max_model(129)
	real max_model_orig(129)
      common /mtemp/ meantemp, meantemp_bs, max_model, max_model_orig

!	Profile for nr = 128 resolution:
!     Mean temperature profile for strongly bottom heated model
!     Schuberth et al. gcube 2009

      data meantemp_bs/&
     & 3.0000000e+02, 7.2118800e+02, 1.0209510e+03, 1.2499570e+03,&
     & 1.4206910e+03, 1.5506250e+03, 1.6490020e+03, 1.7226340e+03,&
     & 1.7784700e+03, 1.8211660e+03, 1.8540640e+03, 1.8809900e+03,&
     & 1.9041580e+03, 1.9215510e+03, 1.9335260e+03, 1.9424690e+03,&
     & 1.9497720e+03, 1.9567120e+03, 1.9634020e+03, 1.9700150e+03,&
     & 1.9766080e+03, 1.9833160e+03, 1.9902840e+03, 1.9975780e+03,&
     & 2.0050360e+03, 2.0127050e+03, 2.0207260e+03, 2.0287030e+03,&
     & 2.0371460e+03, 2.0452620e+03, 2.0534720e+03, 2.0617110e+03,&
     & 2.0696770e+03, 2.0775040e+03, 2.0854490e+03, 2.0933810e+03,&
     & 2.1012680e+03, 2.1092260e+03, 2.1172470e+03, 2.1253440e+03,&
     & 2.1334370e+03, 2.1415970e+03, 2.1499480e+03, 2.1585630e+03,&
     & 2.1672930e+03, 2.1758500e+03, 2.1840060e+03, 2.1917910e+03,&
     & 2.1994580e+03, 2.2073450e+03, 2.2157030e+03, 2.2246140e+03,&
     & 2.2339600e+03, 2.2434830e+03, 2.2528890e+03, 2.2619630e+03,&
     & 2.2706260e+03, 2.2789390e+03, 2.2870550e+03, 2.2951410e+03,&
     & 2.3032670e+03, 2.3113570e+03, 2.3192350e+03, 2.3267530e+03,&
     & 2.3338810e+03, 2.3407070e+03, 2.3473300e+03, 2.3537710e+03,&
     & 2.3599610e+03, 2.3658010e+03, 2.3712310e+03, 2.3762530e+03,&
     & 2.3809140e+03, 2.3852620e+03, 2.3893430e+03, 2.3931860e+03,&
     & 2.3968060e+03, 2.4002120e+03, 2.4034520e+03, 2.4066350e+03,&
     & 2.4098850e+03, 2.4132350e+03, 2.4165720e+03, 2.4196820e+03,&
     & 2.4224000e+03, 2.4247260e+03, 2.4268280e+03, 2.4289490e+03,&
     & 2.4312780e+03, 2.4338770e+03, 2.4366870e+03, 2.4395810e+03,&
     & 2.4424480e+03, 2.4452080e+03, 2.4477990e+03, 2.4501290e+03,&
     & 2.4520790e+03, 2.4535320e+03, 2.4544260e+03, 2.4547630e+03,&
     & 2.4546100e+03, 2.4540690e+03, 2.4532190e+03, 2.4520610e+03,&
     & 2.4505210e+03, 2.4484910e+03, 2.4458500e+03, 2.4424470e+03,&
     & 2.4381390e+03, 2.4329080e+03, 2.4269790e+03, 2.4207860e+03,&
     & 2.4148000e+03, 2.4094040e+03, 2.4049250e+03, 2.4017200e+03,&
     & 2.4001760e+03, 2.4006700e+03, 2.4037220e+03, 2.4101750e+03,&
     & 2.4215510e+03, 2.4408650e+03, 2.4738140e+03, 2.5303970e+03,&
     & 2.6282240e+03, 2.7968510e+03, 3.0753860e+03, 3.4877340e+03,&
     & 4.2000000e+03/
     
!     new 1D profile - 200Ma forward simulation 08/11

      data meantemp/&
     & 300.0, 744.2, 1049.6, 1290.3, 1472.9, 1623.1, 1740.5, 1832.6, 1896.8,&
     & 1942.9, 1970.1, 1986.8, 1995.9, 2002.9, 2008.4, 2014.0, 2018.6,&
     & 2022.6, 2024.8, 2026.1, 2026.7, 2027.1, 2026.6, 2024.9, 2022.1,&
     & 2019.4, 2019.3, 2022.0, 2026.8, 2033.0, 2040.1, 2048.0, 2056.3,&
     & 2065.1, 2074.0, 2082.8, 2091.6, 2100.3, 2109.1, 2118.1, 2127.3,&
     & 2136.7, 2146.4, 2156.1, 2165.9, 2175.5, 2185.1, 2194.6, 2204.1,&
     & 2213.3, 2222.2, 2230.8, 2239.0, 2247.0, 2254.7, 2262.4, 2269.9,&
     & 2277.3, 2284.6, 2291.5, 2298.2, 2304.5, 2310.4, 2316.1, 2321.6,&
     & 2326.9, 2332.0, 2337.0, 2341.9, 2346.9, 2351.9, 2356.8, 2361.4,&
     & 2365.7, 2369.7, 2373.5, 2377.2, 2380.8, 2384.5, 2388.4, 2392.5,&
     & 2396.8, 2401.2, 2405.6, 2409.7, 2413.6, 2417.1, 2420.3, 2423.4,&
     & 2426.4, 2429.2, 2431.8, 2434.2, 2436.3, 2438.4, 2440.4, 2442.3,&
     & 2443.7, 2444.6, 2445.1, 2445.2, 2444.9, 2444.0, 2442.4, 2440.0,&
     & 2436.8, 2432.8, 2428.0, 2422.3, 2415.8, 2408.5, 2400.7, 2392.7,&
     & 2384.6, 2376.8, 2369.6, 2363.4, 2359.1, 2357.9, 2362.0, 2375.2,&
     & 2403.1, 2454.8, 2542.8, 2684.9, 2904.5, 3228.5, 3671.1, 4200.0/
     
     data max_model_orig/&
     & 0.000000000000000, 4.40322899935208, 2.18070020204347,&   
     & 1.49387138770569, 1.31523445506976, 1.25548226010607,&
     & 1.35073906726345, 1.29091639063942, 1.21334844617840,&     
     & 1.16088491964316, 1.13103543137466, 1.11319027986422,&     
     & 1.03959742449051, 0.911191818583535, 0.952412685287734,&     
     & 0.970814886956845, 0.886796097481161, 0.809886614309679,&    
     & 0.782421771288664, 0.736904198946140, 0.725422185555120,&     
     & 0.697428561093099, 0.707463410380533, 0.554896619856250,&     
     & 0.538317480235575, 0.576208520837667, 0.586269433936937,&     
     & 0.562415427557548, 0.518092859856621, 0.538676761232676,&     
     & 0.539978078023776, 0.546231405001737, 0.550002862254062,&     
     & 0.552689835305904, 0.554932410483149, 0.554749283615074,&     
     & 0.555299818184290, 0.557275816868919, 0.557685575063461,&     
     & 0.557338305518007, 0.558903063273774, 0.557943542195629,&     
     & 0.558104119935667, 0.557750269703688, 0.557217160574710,&     
     & 0.556842026452436, 0.556405618204897, 0.555925694678814,&    
     & 0.555512520154437, 0.555204907401728, 0.555026758581976,&     
     & 0.555121908430498, 0.555357021002182, 0.555810717037923,&     
     & 0.556303421521832, 0.556861781774679, 0.557380624732110,&     
     & 0.557946278796568, 0.558640802746153, 0.559624983563125,&     
     & 0.560807737128859, 0.562197158022396, 0.563748929712077,&     
     & 0.565447810308908, 0.567245809965939, 0.569132040106256,&     
     & 0.571068686944112, 0.573024033709890, 0.574944200514859,&     
     & 0.576796683776682, 0.578587204705020, 0.580387595892616,&     
     & 0.582271362934621, 0.584297896042565, 0.586462117395453,&     
     & 0.588725895519433, 0.591018006317821, 0.593270330068960,&     
     & 0.595394410441432, 0.597332379781836, 0.599061885101021,&     
     & 0.600610823978882, 0.602027443314970, 0.603390510741926,&     
     & 0.605017766945782, 0.606904389337933, 0.608946965434143,&     
     & 0.611105292841730, 0.613325135127315, 0.615799024135729,&     
     & 0.618333680182797, 0.620984520707661, 0.623765071725701,&     
     & 0.626617859493677, 0.629468759410066, 0.632320244610123,&     
     & 0.635351517262956, 0.638743452015870, 0.642424664465021,&     
     & 0.646360227841748, 0.650510274118509, 0.654907132095194,&     
     & 0.659661123643794, 0.664899094987077, 0.670710285692425,&     
     & 0.677130644825011, 0.684170011477882, 0.691850034545140,&     
     & 0.700211589971165, 0.709247995038903, 0.718883454709004,&     
     & 0.729302446033525, 0.740062795136938, 0.750885672684951,&     
     & 0.761834829801230, 0.770800383587118, 0.775385057313020,&     
     & 0.778662221620479, 0.779577882334569, 0.776490269101800,&     
     & 0.766682513695426, 0.746145142469305, 0.709404452787919,&     
     & 0.650277011830621, 0.563019587278157, 0.445035298784681,&     
     & 0.300324644967106, 0.143845123746678, 0.000000000000000/
  
     data max_model/&  
     & 0.0000000000E+00, 0.1405402631E+03, 0.7212428231E+02,&
     & 0.5210598249E+02, 0.4249548243E+02, 0.3565785990E+02,&
     & 0.3351367204E+02, 0.3530581509E+02, 0.3627522610E+02,&
     & 0.3937912352E+02, 0.4168502215E+02, 0.4268859015E+02,&
     & 0.4223104157E+02, 0.3958886725E+02, 0.3454559683E+02,&
     & 0.3045620179E+02, 0.3161903752E+02, 0.3201915032E+02,&
     & 0.2921931121E+02, 0.2462457260E+02, 0.2188891574E+02,&
     & 0.1984996793E+02, 0.1984144264E+02, 0.2026092709E+02,&
     & 0.1987296064E+02, 0.1889210995E+02, 0.1863663411E+02,&
     & 0.1909573751E+02, 0.1824021778E+02, 0.1755918716E+02,&
     & 0.1822739631E+02, 0.1883993336E+02, 0.1925143471E+02,&
     & 0.1916378304E+02, 0.1865769206E+02, 0.1974818400E+02,&
     & 0.2004180223E+02, 0.1962899456E+02, 0.1888453324E+02,&
     & 0.1798586465E+02, 0.1717375878E+02, 0.1613442271E+02,&
     & 0.1504476070E+02, 0.1429111449E+02, 0.1373888064E+02,&
     & 0.1337414145E+02, 0.1428319798E+02, 0.1523068300E+02,&
     & 0.1565171111E+02, 0.1576392033E+02, 0.1552121834E+02,&
     & 0.1568875898E+02, 0.1597435024E+02, 0.1611766228E+02,&
     & 0.1617934025E+02, 0.1651039141E+02, 0.1678379562E+02,&
     & 0.1698695065E+02, 0.1746979853E+02, 0.1791182056E+02,&
     & 0.1827088482E+02, 0.1856002894E+02, 0.1899220926E+02,&
     & 0.1933268174E+02, 0.1959325572E+02, 0.1979267746E+02,&
     & 0.2014714354E+02, 0.2063696923E+02, 0.2095780505E+02,&
     & 0.2108893509E+02, 0.2104858753E+02, 0.2129632309E+02,&
     & 0.2146606337E+02, 0.2152336873E+02, 0.2146194283E+02,&
     & 0.2126496724E+02, 0.2091901927E+02, 0.2048148257E+02,&
     & 0.2028541171E+02, 0.2013233958E+02, 0.1988468868E+02,&
     & 0.1953228426E+02, 0.1923494157E+02, 0.1915068890E+02,&
     & 0.1890306109E+02, 0.1847950667E+02, 0.1786944999E+02,&
     & 0.1718348695E+02, 0.1772430360E+02, 0.1801655690E+02,&
     & 0.1867579592E+02, 0.1950244078E+02, 0.2106814343E+02,&
     & 0.2199249942E+02, 0.2247503000E+02, 0.2242977836E+02,&
     & 0.2300585778E+02, 0.2356024086E+02, 0.2417692997E+02,&
     & 0.2436453170E+02, 0.2423462239E+02, 0.2604776142E+02,&
     & 0.2813714309E+02, 0.2970270338E+02, 0.3061979659E+02,&
     & 0.3095866571E+02, 0.3071525385E+02, 0.2986099729E+02,&
     & 0.2914199477E+02, 0.2819333046E+02, 0.2765281057E+02,&
     & 0.2889373747E+02, 0.2995938812E+02, 0.3131601788E+02,&
     & 0.3210715927E+02, 0.3365924746E+02, 0.3548225610E+02,&
     & 0.3741141719E+02, 0.3939893890E+02, 0.4129554516E+02,&
     & 0.4292100130E+02, 0.4514407381E+02, 0.4736895387E+02,&
     & 0.4834306098E+02, 0.4706403546E+02, 0.4158905797E+02,&
     & 0.3003246449E+02, 0.1438451237E+02, 0.0000000000E+00/
     
end blockdata
    