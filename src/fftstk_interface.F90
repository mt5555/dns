#include "macros.h"
#if 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! our wrapper for stk
! 
!
! provides public interfaces:
!   fft_interface_init               call this before using any other routines
!   fft1                             fft along first dimension of 3D array
!   ifft1	                     ifft along first dimension of 3D array
! 
! Routines work on data of the form:  p(n1d,n2d,n3d)
! Size of the grid point data         p(1:n1,1:n2,1:n3)
! Size of fourier coefficients        p(1:n1,1:n2,1:n3)
!
! grid space data:    1 2 3 4 5 6 7 8 * 
! 
! fourier space:      * 0 1 1 2 2 3 3 4 
! rearranged:         0 4 1 1 2 2 3 3 *
!

!


FFT data representation:

sum over m=1..n/2:

   f = fhat(1)  +  2 fhat(2*m) cos(m*2pi*x) - 2*fhat(2m+1) sin(m*2pi*x)


     if isign = +1, and m coefficient vectors are supplied
     each containing the sequence:

     a(0),a(1),b(1),...,a(n/2)  (n values)

     then the result consists of m data vectors each
     containing the corresponding n gridpoint values:

     x(0), x(1), x(2),...,x(n-1)

     note: the fact that the gridpoint values x(j) are real
     implies that b(0)=b(n/2)=0.  for a call with isign=+1,
     it is not actually necessary to supply these zeros.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif


module fft_interface
use params
implicit none
integer, parameter ::  num_fftsizes=3

integer :: init=0
type fftdata_d
   CPOINTER :: ptrigs
   integer :: size
end type
type(fftdata_d) :: fftdata(num_fftsizes)


private :: fftinit, getindex
contains 



subroutine fft_interface_init(f,work,nx,ny,nz)
integer :: nx,ny,nz
real*8 :: f(nx,ny,nz)
real*8 :: work(nx,ny,nz)
integer :: i,n1,n1d,n2,n2d,n3,n3d,index

init=1
do i=1,num_fftsizes
   fftdata(i)%size = 0	
enddo
end subroutine


subroutine fft_get_mcord(mcord,n)
integer n,mcord(:)
integer i
do i=1,n
   mcord(i)=(i-1)/2	
   if (mod(i,2)==0) mcord(i)=-mcord(i)
   if (i==2) mcord(i)=n/2
enddo
end subroutine







subroutine fftinit(n,index)
integer n,index
character*80 message

if (init==0) call abortdns("fftstk_interface.F90: call fft_interface_init to initialize first!");
if (n>1000000) call abortdns("fftstk_interface.F90: n>1 million")

fftdata(index)%size=n

write(message,'(a,i6)') 'Initializing stk FFT of size n=',n
call print_message(message)

call rfft_init(n,fftdata(index)%ptrigs)

end subroutine




subroutine getindex(n1,index)
integer n1,index

character*80 message_str
integer i,k


i=0
do 
   i=i+1
   if (i>num_fftsizes) then
      write(message_str,'(a,i10)') "fftstk_interface.F90:  Failed initializing an fft of size =",n1
      call abortdns(message_str)
   endif

   if (fftdata(i)%size==0) then
      call fftinit(n1,i)      
      exit 
   endif
   if (n1==fftdata(i)%size) exit 
enddo
index=i
end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute 3D in-place iFFT of p
! FFT taken along first direction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ifft1(p,n1,n1d,n2,n2d,n3,n3d)
integer n1,n1d,n2,n2d,n3,n3d
real*8 p(n1d,n2d,n3d)

real*8 tmx1,tmx2
integer index,i,k,j

call wallclock(tmx1)
if (tims(19)==0) ncalls(19)=0  ! timer was reset, so reset counter too

if (n1==1) return
ASSERT("ifft1: dimension too small ",1+n1<=n1d);
call getindex(n1,index)

do k=1,n3
#if 0
   do j=1,n2
      call rfft_synthesis2(p(1,j,k),fftdata(index)%ptrigs)
   enddo
#else
   call rfft_synthesis_M2(p(1,1,k),fftdata(index)%ptrigs,n2,n1d)
#endif
enddo

call wallclock(tmx2) 
tims(19)=tims(19)+(tmx2-tmx1)          
ncalls(19)=ncalls(19)+1
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute 3D in-place FFT of p
! FFT taken along first direction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft1(p,n1,n1d,n2,n2d,n3,n3d)
integer n1,n1d,n2,n2d,n3,n3d
real*8 p(n1d,n2d,n3d)

integer index,i,j,k
real*8 tmx1,tmx2

call wallclock(tmx1)
if (tims(18)==0) ncalls(18)=0  ! timer was reset, so reset counter too

if (n1==1) return
ASSERT("fft1: dimension too small ",1+n1<=n1d);
call getindex(n1,index)


do k=1,n3
#if 0
   do j=1,n2   
      call rfft_analysis2(p(1,j,k),fftdata(index)%ptrigs)
      do i=1,n1
         p(i,j,k)=p(i,j,k)/n1
      enddo
   enddo
#endif
   call rfft_analysis_M2(p(1,1,k),fftdata(index)%ptrigs,n2,n1d)
   do j=1,n2   
      do i=1,n1
         p(i,j,k)=p(i,j,k)/n1
      enddo
   enddo
enddo
call wallclock(tmx2) 
tims(18)=tims(18)+(tmx2-tmx1)          
ncalls(18)=ncalls(18)+1
end subroutine








end ! module mod_fft_interface
