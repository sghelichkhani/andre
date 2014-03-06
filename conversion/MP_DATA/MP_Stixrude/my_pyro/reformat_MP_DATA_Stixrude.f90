program reformat_MP_DATA_Stixrude


implicit none

integer i,j

integer, parameter :: n_p=1401
integer, parameter :: n_T=471

double precision, dimension(1:n_p) :: p,d
double precision, dimension(1:n_T) :: T
double precision, dimension(1:n_p,1:n_T) :: vs,vp,vc,rho,vs_Q,vp_Q,Q_p,Q_s

character(len=4) T_str


p=0.
d=0.
T=0.
vs=0.
vp=0.
vc=0.
rho=0.
vs_Q=0.
vp_Q=0.

do i=1,n_T

   write(T_str,'(i4)') 300 + (i-1)*10

   write(*,*) T_str

   !if (T_str(1:1) /= '2') then
   
   open(unit=801,file='my_pyro.'//adjustl(T_str),status='old')

   do j=1,n_p

      read(801,15) p(j),d(j),T(i),rho(j,i),vc(j,i),vs(j,i),vp(j,i),vs_Q(j,i),vp_Q(j,i)!,Q_s(j,i),Q_p(j,i)

15 format(f6.2,f9.2,f9.2,6f11.5)
      
      !write(*,15) p(j),d(j),T(i),rho(j,i),vc(j,i),vs(j,i),vp(j,i),vs_Q(j,i),vp_Q(j,i)!,Q_s(j,i),Q_p(j,i)

   enddo

   close(801)

   !endif

enddo

open(unit=801,file='MP_Stixrude_P.dat',status='unknown')
open(unit=803,file='MP_Stixrude_depth.dat',status='unknown')
open(unit=804,file='MP_Stixrude_T.dat',status='unknown')
open(unit=805,file='MP_Stixrude_vs.dat',status='unknown')
open(unit=806,file='MP_Stixrude_vp.dat',status='unknown')
open(unit=807,file='MP_Stixrude_rho.dat',status='unknown')

open(unit=808,file='MP_Stixrude_vs_Q.dat',status='unknown')
open(unit=809,file='MP_Stixrude_vp_Q.dat',status='unknown')

write(804,'(471f9.2)') T

do j=1,n_p

  write(801,'(f6.2)') p(j)
  write(803,'(f9.2)') d(j)
  write(805,16) vs(j,:)
  write(806,16) vp(j,:)
  write(807,16) rho(j,:)

  write(808,16) vs_Q(j,:)
  write(809,16) vp_Q(j,:)

16 format(471f11.5)
  
enddo

close(801)
close(803)
close(804)
close(805)
close(806)
close(807)
close(808)
close(809)

end program reformat_MP_DATA_Stixrude

