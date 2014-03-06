program SHconvert_tet
implicit none

include 'size.h'
include 'pcom.h'

integer xxx, psp, byebye, idump, pro
integer i,j,k,ir, id, ii,locsize, maxls, tzt, nproc2
integer x,y,z, nbin, ibin, indx, ind(nr+1)

integer, parameter:: nn=(nt+1)**2*nd      
integer, parameter:: mm=(mt+1)**2*nd
integer, parameter:: lmax=50

real, parameter:: PI = 3.141592653589793
real, parameter:: RHOAV = 5514.3d0		 ! average density in the full Earth to normalize equation
real, parameter:: GRAV = 6.6723d-11		 ! gravitational constant

! all vectors are "real" because of the very large size
real rmean_T(nr+1), rmean_rho(nr+1), rshl(nr+1)
real rmean_rho_loop(nr+1), border_pos, border_neg, rmean_T_loop(nr+1)
real f(2,(lmax+1)*(lmax+2)/2,nr+1), f2(2,(lmax+1)*(lmax+2)/2,nr+1)
real f3(2,(lmax+1)*(lmax+2)/2,nr+1), f4(2,(lmax+1)*(lmax+2)/2,nr+1)
real f_tmp(2,(lmax+1)*(lmax+2)/2,nr+1), f2_tmp(2,(lmax+1)*(lmax+2)/2,nr+1)
real f3_tmp(2,(lmax+1)*(lmax+2)/2,nr+1), f4_tmp(2,(lmax+1)*(lmax+2)/2,nr+1)
real, allocatable:: s2(:,:), s2_rho(:,:), u(:,:), u_ampl(:,:), s2_T(:,:)
real, allocatable:: shrh(:,:)
real vp(1,nn), vs(1,nn), rho(1,nn), depth_crit, depth_cut
real depth(nr+1), meantemp(129), meantemp_bs(129), mnb
real rmax(nr+1), rmax_loop(nr+1), rmax_loop_proc(nr+1)
real, allocatable:: bin(:), bin_cnt(:), bin_cnt_tot(:)
real val_pos, val_neg, val, R_min, R_max, quot

 character MP_TYP, comp*12, gpath*100, lpath*100, cname*3, cstage*4, cstage2*2, fname*40, char1*4
 character titl(4,4)*8, cstage3*3, cstage_out*4, gpath2*100, iter*2, iter_out*2

! ################################################
! ################# MAIN PART #######################
! ################################################

common /mtemp/ meantemp, meantemp_bs
      
! Initialize parallel communications.
      
 CALL MPI_INIT(ierr)
 CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid2, ierr)
 CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc2, ierr)
nproc=(mt/nt)**2*10/nd
lvproc = 1.45*log(real(mt/nt))
mproc = 2**(2*lvproc)
 
allocate(s2(nn*(nr+1),nproc/nproc2),s2_rho(nn*(nr+1),nproc/nproc2))
allocate(u(nn*3*(nr+1),nproc/nproc2),u_ampl(nn*(nr+1),nproc/nproc2))
allocate(s2_T(nn*(nr+1),nproc/nproc2),shrh(nn*(nr+1),nproc/nproc2))

border_pos=1000000000.0
border_neg=-1000000000.0
depth_crit=-1.0
depth_cut=-1.0
R_max=6.370e+06 
R_min=3.480e+06
quot=128.0/real(nr) 	! scaling to 'meantemp' array
do ir=1,nr+1
	ind(ir)=int(quot*(ir-1))+1
	rshl(ir)=R_min+(R_max-R_min)*real(nr+1-ir)/real(nr)
enddo

! ######### Histograms ###########
nbin=100	! number of bins
allocate(bin(0:nbin),bin_cnt(nbin*(nr+1)),bin_cnt_tot(nbin*(nr+1)))
val_pos=50.0
val_neg=-50.0
do ii=0,nbin
	bin(ii)=val_neg+ii*(val_pos-val_neg)/real(nbin)
enddo
! ############################


! multiple conversion
do psp=1,3

!gpath='~/'
gpath='../'
gpath='/SCRATCH/horbach/381/'
!gpath='/import/netapp-terra/horbach/adjoint/runs/'
!gpath='/import/netapp-terra/horbach/Bernhard/'
!gpath='/import/netapp-terra/horbach/adjoint/Endfeld/'
!gpath='/import/netapp-terra/bernhard/NEW_MCM/'
!gpath='/import/netapp-terra/horbach/Adjungierte/'

write(iter,'(I2.2)') psp-1
write(iter_out,'(I2.2)') psp-1

 cname="381"
byebye=11
!lpath='c-files_'//iter
lpath='c-files' 

do xxx=1,byebye

rmean_rho=0.0
rmean_T=0.0
	
do pro=1,nproc/nproc2
	
	myid=myid2*nproc/nproc2+pro-1

 	write(char1,'(I4.4)') myid
	
	idump=xxx-1
		
	write(cstage2,'(i2.2)') idump
	write(cstage,'(i4.4)') idump
	write(cstage3,'(i3.3)') idump
	write(cstage_out,'(i4.4)') 50*(xxx-1)
	
      if(xxx>410) then
      	fname = 'a'//cname//'.'//char1//'.'//cstage
      	open(188, file=trim(gpath)//trim(lpath)//'/'//trim(fname), form='unformatted', status='unknown')

     		call vecinunform44(s2(:,pro),188)
        else
         	fname='c'//cname//'.'//char1//'.'//cstage2
         	!fname='c'//cname//'.'//char1//'.'//cstage3
         	!fname='s'//cname//'.'//char1//'.'//cstage
         	!fname='tomo.'//char1//'.'//cstage2
         	!fname='g-files/g210.'//char1//'.'//cstage3
        	!open(188, file=trim(gpath)//trim(lpath)//'/output/'//trim(fname), form='formatted', status='unknown')
   	  	open(188, file=trim(gpath)//trim(lpath)//'/'//trim(fname)//'_'//iter, form='formatted', status='unknown')
    		call vecin11(s2(:,pro),188,1)
    		call vecin11(u(:,pro),188,3)
		!call vecin11(shrh(:,pro),188,1)
	endif
	close(188)
 
! 1) mineralogical model only available between 300K and 4200K
! 2) T_CMB = 4200K
k=1
do ir=1,nr+1
	do id=1,nd
		do ii=1,(nt+1)**2
			if(s2(k,pro)>4200.0) s2(k,pro)=4200.0
			if(s2(k,pro)<300.0) s2(k,pro)=300.0
			!if(not(s2(k,pro)>=300.0.or.s2(k,pro)<=4200.0)) then
			!	s2(k,pro)=meantemp(ind(ir))
			!endif	
	         k=k+1
		enddo
	enddo
enddo

 call read_MP_DATA

do ir=1,nr+1
	! 1 additional kilometer to fit the range of the min. models
	depth(ir)=rshl(1)+1000.0-rshl(ir)
	i = (ir-1)*nn+1
	j = ir*nn

	s2_T(i:j,pro)=s2(i:j,pro)
	call convert_T(depth(ir),1,s2(i:j,pro),nn,vp,vs,rho)
	s2_rho(i:j,pro)=rho(1,:)
enddo

 call rad_mean(s2_rho(:,pro),rmean_rho_loop)
 call rad_mean(s2_T(:,pro),rmean_T_loop)
rmean_rho=rmean_rho+rmean_rho_loop
rmean_T=rmean_T+rmean_T_loop

!call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!call system ('rm '//trim(gpath)//trim(lpath)//'/'//trim(fname)//'_'//iter)
!call MPI_BARRIER(MPI_COMM_WORLD,ierr)

enddo !proc

 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
f=0.0
f2=0.0
f3=0.0
f4=0.0
bin_cnt=0.0
rmax=0.0

do pro=1,nproc/nproc2

myid=myid2*nproc/nproc2+pro-1

!border_pos=5.0
border_neg=-0.5
depth_crit=150.0
!depth_cut=150.0

rmax_loop=0.0
rmax_loop_proc=0.0

k=1
do ir=1,nr+1
	do id=1,nd
		do ii=1,(nt+1)**2
			s2_rho(k,pro)=(s2_rho(k,pro)-rmean_rho(ir))/rmean_rho(ir)
			s2_T(k,pro)=(s2_T(k,pro)-rmean_T(ir))/rmean_T(ir)
			if(depth(ir)<=depth_cut*1000.0) then
				s2_rho(k,pro)=0.0
			else if(depth(ir)<=depth_crit*1000.0) then
				if(s2_rho(k,pro)*100<border_neg) then
					s2_rho(k,pro)=border_neg/100.0
				elseif(s2_rho(k,pro)*100>border_pos) then
					s2_rho(k,pro)=border_pos/100.0
				endif
			endif
			
			x=(ir-1)*3*nd*(nt+1)**2+(id-1)*(nt+1)**2+ii
			y=(ir-1)*3*nd*(nt+1)**2+nd*(nt+1)**2+(id-1)*(nt+1)**2+ii
			z=(ir-1)*3*nd*(nt+1)**2+2*nd*(nt+1)**2+(id-1)*(nt+1)**2+ii
			u_ampl(k,pro)=sqrt(u(x,pro)**2.0+u(y,pro)**2.0+u(z,pro)**2.0)
			
			! ######## histograms ########
!			val=s2_T(k,pro)*100.0
!			do ibin=1,nbin
!				indx=(ir-1)*nbin+ibin
!				if(val>=bin(ibin-1).and.val<bin(ibin)) then
!					bin_cnt(indx)=bin_cnt(indx)+1
!				endif
!			enddo
			! ########################
						
			if(s2_T(k,pro)>rmax_loop_proc(ir)) rmax_loop_proc(ir)=s2_T(k,pro)
			
			k=k+1
		enddo
	enddo
enddo

 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 call gridinit
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 !call MPI_REDUCE(bin_cnt,bin_cnt_tot,nbin*(nr+1),MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
 !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 
 call fetosh(s2_rho(:,pro),f_tmp,lmax)
 call fetosh(s2_T(:,pro),f2_tmp,lmax)
 !call fetosh(u_ampl(:,pro),f3_tmp,lmax)
 !call fetosh(shrh(:,pro),f4_tmp,lmax)
f=f+f_tmp
f2=f2+f2_tmp
!f3=f3+f3_tmp
!f4=f4+f4_tmp

 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 call MPI_REDUCE(rmax_loop_proc,rmax_loop,nr+1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(rmax_loop,(nr+1),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
! determine maximal variations in each layer
! and respective scaling factor to reduce variations in case of a too high temperature
do ir=1,nr+1
	if(rmax_loop(ir)>rmax(ir)) rmax(ir)=rmax_loop(ir)
enddo
 
enddo !proc

if(myid2==0) then
	!open(unit=314,file='/home/horbach/Work/NewPlates/'//trim(lpath)//'/'//cstage//'.dat',status='replace')
	!open(unit=314,file='/home/horbach/Work/TERRA/Adjoint_Test/'//cstage//'.dat',status='replace')
      !open(unit=314,file='/home/horbach/TERRA/Adjoint/Endfeld/'//trim(lpath)//'/'//cstage//'.dat',status='replace')
      !open(unit=314,file='/home/horbach/TERRA/Adjoint/runs/'//trim(lpath)//'/'//cstage//'.dat',status='replace')
	open(unit=314,file=trim(gpath)//'/dat-files/Iteration_'//iter_out//'/'//trim(cstage_out)//'_50.dat',status='replace')
	!open(unit=314,file=trim(gpath)//'/'//trim(cstage_out)//'_50.dat',status='replace')
	
     	write(314,*) lmax, nr+1, cname
      do i=1,nr+1
      	write(314,'(E20.12)') rshl(i)
      enddo
      do i=1,nr+1
      	write(314,'(E20.12)') rmean_rho(i)*1000.0
      enddo
      do i=1,nr+1
      	do j=1,(lmax+1)*(lmax+2)/2
      		write(314,'(2E18.10)') f(1,j,i), f(2,j,i)
      	enddo
      enddo
	do i=1,nr+1
		write(314,'(E20.12)') rmean_T(i)
	enddo
	do i=1,nr+1
		do j=1,(lmax+1)*(lmax+2)/2
			write(314,'(2E18.10)') f2(1,j,i), f2(2,j,i)
		enddo
	enddo
	!do i=1,nr+1
      	!do j=1,(lmax+1)*(lmax+2)/2
      	!	write(314,'(2E18.10)') f3(1,j,i), f3(2,j,i)
      	!enddo
      	!enddo
	!do i=1,nr+1
	!do j=1,(lmax+1)*(lmax+2)/2
	!	write(314,'(2E18.10)') f4(1,j,i), f4(2,j,i)
	!enddo
	!enddo
	close(314)

	! ####### write histogram data to file ###############
	!open(unit=315,file=trim(gpath)//trim(lpath)//'/hist'//trim(cstage_out)//'.dat',status='replace')
	!write(315,*) nbin, val_neg, val_pos
	!write(315,*) nr-1
      !do i=2,nr
      !	write(315,'(E20.12)') rshl(i)
      !enddo
      !do i=0,nbin-1
      !	write(315,'(E20.12)') bin(i)
      !enddo
      !do i=nbin+1,nbin*nr
      !	write(315,'(E20.12)') bin_cnt_tot(i)
      !enddo
      
	!close(315)
	
	! ####### write max-values to file ###############
	!open(unit=316,file=trim(gpath)//trim(lpath)//'/max_val',status='replace')
	!write(316,*) rmax
	!close(316)
endif
  
if(myid2==0) write(*,*) cstage_out, " done"

!call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!call system('rm '//trim(gpath)//trim(lpath)//'/'//trim(fname)//'_'//iter)
!call MPI_BARRIER(MPI_COMM_WORLD,ierr)

enddo		! multiple
enddo		! conversion
 
! leave parallel communication before exiting
if(myid2== 0) then
	write(*,*)
	write(*,*) "Servus!"
endif
 call MPI_FINALIZE(ierr)

end program SHconvert_tet


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


subroutine vecinunform44(u,nf)
implicit none
!...  This routine reads the nodal field u from logical unit nf
!     using unformated I/O and without any additional information.
 
      include 'size.h'
      include 'pcom.h'
      
      integer i, nf
      real u((nt+1)**2*nd*(nr+1))
 	
	read(nf) (u(i),i=1,(nt+1)**2*nd*(nr+1))
 	
end subroutine vecinunform44
      
      
subroutine vecin11(u,nf,nj)
implicit none
!...  This routine reads the nodal field u from logical unit nf
!...  using 1p,e10.3 format when ifmt = 0 and f10.3 when ifmt = 1.
 
      include 'size.h'
      include 'pcom.h'
      
      integer kr, kt, ii, nj, nf
      real u((nt+1)**2*nd*nj*(nr+1)), rshl(nr+1), propr(20)
      common /prty/ propr
      character*8 titl(4,4)
 
      read(nf,10) kr, kt
 10   format(2i5)
 
      if(kr.ne.nr .or. kt.ne.nt) then
         print *,'Array dimensions of data file and size.h inconsistent'
         ! LMU BS BS begin section added
         print *,'kr = ',kr,' and nr = ',nr
         print *,'kt = ',kt,' and nt = ',nt
         ! LMU BS BS end section added
         stop
      end if
 
      read(nf,20) titl
 20   format(4a8)
 
      read(nf,30) (rshl(ii),ii=1,nr+1)
      read(nf,30) propr
 30   format(1p10e15.8)
 
         read(nf,50) (u(ii),ii=1,(nt+1)**2*nd*nj*(nr+1))
50      format(15f10.3)

end subroutine vecin11

 
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
		MP_dir='./MP_DATA/MP_Stixrude/Pyrolite'
	elseif(comp=='PICLOGITE') then
		MP_dir='./MP_DATA/MP_Stixrude/Piclogite'
	else
		write(*,*) "Desired composition not available!"
		return
	endif

	open(unit=803,file=trim(MP_dir)//'/MP_Stixrude_Pyrolite_depth.dat',action='read',status='old')
	open(unit=804,file=trim(MP_dir)//'/MP_Stixrude_Pyrolite_T.dat',action='read',status='old')
	open(unit=805,file=trim(MP_dir)//'/MP_Stixrude_Pyrolite_vs.dat',action='read',status='old')
	open(unit=806,file=trim(MP_dir)//'/MP_Stixrude_Pyrolite_vp.dat',action='read',status='old')
	open(unit=807,file=trim(MP_dir)//'/MP_Stixrude_Pyrolite_rho.dat',action='read',status='old')

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
		MP_dir='./MP_DATA/MP_Antonio/Pyrolite'
	elseif(comp=='PICLOGITE') then
		MP_dir='./MP_DATA/MP_Antonio/Piclogite'
	else
		write(*,*) "Desired composition not available!"
		return
	endif

	open(unit=804,file=trim(MP_dir)//'/MP_Antonio_Pyrolite_Depth.dat',status='old')
	do i=1,nD_A
	read(804,*) depth_A(i)
	enddo
	close(unit=804)

	open(unit=804,file=trim(MP_dir)//'/MP_Antonio_Pyrolite_T.dat',status='old')
	read(804,*) (temp_A(i),i=1,nT_A)
	close(unit=804)

	open(unit=804,file=trim(MP_dir)//'/MP_Antonio_Pyrolite_G.dat',status='old')
	do j=1,nD_A
	read(804,*) (G_A(j,i),i=1,nT_A)
	enddo
	close(unit=804)

	open(unit=804,file=trim(MP_dir)//'/MP_Antonio_Pyrolite_K.dat',status='old')
	do j=1,nD_A
	read(804,*) (K_A(j,i),i=1,nT_A)
	enddo
	close(unit=804)

	open(unit=804,file=trim(MP_dir)//'/MP_Antonio_Pyrolite_rho.dat',status='old')
	do j=1,nD_A
	read(804,*) (rho_A(j,i),i=1,nT_A)
	enddo
	 close(unit=804)

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

!...  Generate nodal coordinate array xn for finest level. (This
!...  duplicates part of array xm, but is done for convenience
!...  because xn does not require use of a pointer.)
 
      mm = (mt + 1)**2
      m0 = mm*30 + 1
 
      call grdgen(a,mt)
 
      if(nproc .eq. 1) then
         call scopy((mt+1)**2*30, a, 1, xn, 1)
      else
         call subarray(a,xn,0,10,nd,mt,nt,3)
      endif
 
!...  Generate the nodal area arrays.
 
      call ndarea(a(m0),a(m0+mm),a(m0+2*mm),a(m0+3*mm),a(m0+4*mm),a,mt)
 
      if(nproc .eq. 1) then
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


blockdata bdconvect

	real meantemp(129), meantemp_bs(129)
      common /mtemp/ meantemp, meantemp_bs

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
     
!     new geotherm 200Ma forward simulation 08/11

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
   
   
   
end blockdata
    
