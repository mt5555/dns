#include "macros.h"

subroutine abort(message)
use mpi
use params
implicit none
integer ierr
character*(*) message
character*15 :: pre="ASSERT FAILURE "
write(*,'(a)') pre // message

#ifdef IBM
call flush_(6)
#else
call flush(6)
#endif

#ifdef USE_MPI
   call MPI_abort(comm_3d,1,ierr)
#endif
stop
end subroutine


subroutine print_message(message)
use params
implicit none
character*(*) message
if (my_pe==io_pe) then
   write(*,'(a)') message
endif
end subroutine


subroutine wallclock(tmx)
use params
implicit none
integer count,count_rate,count_max
real*8 tmx
#ifdef MPI
real*8 MPI_Wtime
tmx = MPI_Wtime()
#else
call system_clock(count,count_rate,count_max)
tmx=real(count,r8kind)/real(count_rate,r8kind)
#endif

end subroutine





subroutine plotASCII(spectrum,n,title)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ASCII plot of spectrum(0:n) on a log-log scale
!
! Y-axis scale is  10^-10 ... 10^0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use params
implicit none
integer :: n
real*8  ::  spectrum(0:n)
character*(*) :: title


! local variables
integer,parameter :: numx=20,numy=13
integer i,j
character :: plot(0:numx,0:numy)
real*8 cspec(0:numx)
integer ix,iy


if (my_pe==io_pe) then
cspec=0
plot=" "

! bin all values
do i=0,n
   if (i==0) then
      ix=0
   else
      ix=1 + ( (numx-1)*log10(real(i))/log10(real(n)) )     ! ranges from 1..numx
   endif
   if (ix<0) ix=0
   if (ix>numx) ix=numx
   cspec(ix)=cspec(ix)+spectrum(i)
enddo

! scale from 0..numy, log scale
do i=0,numx
   iy = -(numy/10.0) * log10(1d-200+cspec(i))  ! ranges from 0..numy
   if (iy<0) iy=0
   if (iy>numy) iy=numy
   cspec(i)=iy
enddo

do i=0,numx
    j=cspec(i)
    plot(i,j)="*"
enddo

print *
print *,"1E0   |",plot(:,0),title
do i=1,numy-1
   print *,"      |",plot(:,i)
enddo
print *,"1E-10 |",plot(:,numy)
print *,"      +---------------------+"
write (*,'(a,i4)') "      k=0                 k=",n

endif
end subroutine






