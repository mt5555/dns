#include "macros.h"

subroutine abort(message)
use mpi
use params
implicit none
integer ierr
character(len=*) message
character(len=15) :: pre="ABORT"
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
character(len=*) message
if (my_pe==io_pe) then
   write(*,'(a)') message
endif
end subroutine


subroutine wallclock(tmx)

#ifdef USE_MPI
use params
use mpi
implicit none
real*8 tmx
tmx = MPI_Wtime()

#else
use params
implicit none
integer count,count_rate,count_max
real*8 tmx
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
character(len=*) :: title


! local variables
integer,parameter :: numx=20,numy=13
integer i,j
character :: plot(0:numx,0:numy)
real*8 cspec(0:numx)
real*8 :: delta,ireal
integer ix,iy


if (my_pe==io_pe) then
cspec=0
plot=" "

#if 0
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
#endif

! interpolate
cspec(0)=spectrum(0)
do ix=1,numx
   ireal = 10**( log10(real(n)) * (ix-1)/real(numx-1) )
   i=floor(ireal)
   delta=ireal-i
   if (i>=0 .and. i<=n-1) then
      cspec(ix)=spectrum(i)*(1-delta)+spectrum(i+1)*delta
   endif
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



subroutine random_data(buf,n)
implicit none
integer n
real*8 :: buf(n)
!call random_number(buf)
call gaussian(buf,n)
end subroutine




!.... Gaussian Random number generator ran1 from Numerical recipes
subroutine gaussian(buf,n)
implicit none
      
integer n
real*8 :: buf(n)
real*8 ::  fac,rsq,ran1(2)
integer i,j1,j2

! n=10:  i=1..5   j1max=9, j2max=10 
! n=11:  i=1..6   j1max=11, j2max=12(ignored)
ran1=2*ran1-1
do i=1,(n+1)/2
   j1=2*i-1
   j2=j1+1
   do 
      call random_number(ran1)
      ran1=2*ran1-1
      rsq = ran1(1)**2 + ran1(2)**2
      if (rsq<1 .and. rsq>0) exit
   enddo
   fac = Sqrt(-2.0 * Log(rsq)/rsq)
   buf(j1) = ran1(1) * fac
   if (j2<=n) buf(j2)= ran1(2) * fac
enddo

end subroutine

