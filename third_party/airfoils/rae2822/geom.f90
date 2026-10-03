   program geom

!  Reads the text file containing the coordinates of the 
!  RAE 2822 airfoil and creates input files for ICEM CFD.

   implicit none

   character :: cline * 70 

   integer :: n

   real :: x, yl, yu

   open ( unit=7, file='geom.txt', form='formatted' )
   open ( unit=8, file='yl.pts', form='formatted' )
   open ( unit=9, file='yu.pts', form='formatted' )

   do  n = 1, 5
     read(7,'(a70)') cline
   enddo

   write(8,*) 65, 1
   write(9,*) 65, 1

   do  n = 1, 65
     read(7,*) x, yl, yu
     yl = - yl
     write(8,*) x, yl, 0.0
     write(9,*) x, yu, 0.0
   enddo

   stop
   end program geom
