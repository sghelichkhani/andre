program reformat

double precision, dimension(290,81) :: K

open(unit=801,file='G_f.dat')

do i=1,290
read(801,*) K(i,:)
enddo

close(801)

open(unit=801,file='G_new_f.dat')

do i=1,290
write(801,'(81f21.13)') K(i,:)
enddo

close(801)





end program reformat

