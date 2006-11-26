#include "macros.h"
#if 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! our wrapper for CXML FFT
! 
!
! provides public interfaces:
!   fft_interface_init               call this before using any other routines
!   fft1                             fft along first dimension of 3D array
!   ifft1	                     ifft along first dimension of 3D array
! 
! Routines work on data of the form:  p(n1d,n2d,n3d)
! Size of the grid point data         p(1:n1,1:n2,1:n3)
! Size of fourier coefficients        p(1:n1+2,1:n2+2,1:n3+2)
!

FFT data representation:

sum over m=1..n/2:

   f = fhat(1)  +  2 fhat(2*m+1) cos(m*2pi*x) - 2*fhat(2m+2) sin(m*2pi*x)


     if isign = +1, and m coefficient vectors are supplied
     each containing the sequence:

     a(0),b(0),a(1),b(1),...,a(n/2),b(n/2)  (n+2 values)

     then the result consists of m data vectors each
     containing the corresponding n+2 gridpoint values:

     x(0), x(1), x(2),...,x(n-1),0,0.

     note: the fact that the gridpoint values x(j) are real
     implies that b(0)=b(n/2)=0.  for a call with isign=+1,
     it is not actually necessary to supply these zeros.


In otherwords:

 grid space data:    1 2 3 4 5 6 7 8 * *
 
 fourier space:      0 0 1 1 2 2 3 3 4 4
 rearranged:         0 4 1 1 2 2 3 3 * *


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif


module fft_interface
implicit none
integer, parameter ::  num_fftsizes=3



integer :: init=0
type fftdata_d
   ! Real cpq data structure is about 160 words, but lets go 10x to be safe!
   ! we are supposed to pass in a named common block from a cmxl include
   ! file - but that does not allow multiple fft's!
   real*8 :: trigs(1600)
   integer :: size
end type
type(fftdata_d) :: fftdata(num_fftsizes)


private :: fftinit, getindex
contains 



subroutine fft_interface_init()
integer :: i

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
integer n,index,stat
character*80 message
#include <cxmldef.for>

if (init==0) call abortdns("fftcpq_interface.F90: call fft_interface_init to initialize first!");

write(message,'(a,i6)') 'Initializing CXML FFT of size n=',n
call print_message(message)

fftdata(index)%size=n
!stat=dfft_init_grp(n,fftdata(index)%trigs,.false.,n)
stat=dfft_init(n,fftdata(index)%trigs,.true.)




end subroutine




subroutine getindex(n1,index)
integer n1,index

character*80 message_str
integer i,k


i=0
do 
   i=i+1
   if (i>num_fftsizes) then
      write(message_str,'(a,i10)') "fft_interface.F90:  Failed initializing an fft of size =",n1
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
#include <cxmldef.for>
integer n1,n1d,n2,n2d,n3,n3d
real*8 p(n1d,n2d,n3d)

integer index,i,j,k,numffts
real*8 :: scale
if (n1==1) return
ASSERT("ifft1: dimension too small ",n1+2<=n1d)
call getindex(n1,index)
scale=n1


do k=1,n3
   do j=1,n2

      do i=1,n1
         p(i,j,k)=p(i,j,k)*scale
      enddo

!     move the last cosine mode back into correct location:
      p(n1+1,j,k)=p(2,j,k)
      p(2,j,k)=0             ! not needed?
      p(n1+2,j,k)=0          ! not needed?

      call dfft_apply('C','R','B',p(1,j,k),p(1,j,k),fftdata(index)%trigs,1)
   enddo

!   j=dfft_apply_grp('C','R','B',p(1,1,k),p(1,1,k),n2,1,&
!           fftdata(index)%trigs,1,n1d)
!   ASSERT("Error: ifft1() status <> 0",j==0)

enddo

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute 3D in-place FFT of p
! FFT taken along first direction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft1(p,n1,n1d,n2,n2d,n3,n3d)
#include <cxmldef.for>
integer n1,n1d,n2,n2d,n3,n3d
real*8 p(n1d,n2d,n3d)
real*8 :: scale

integer index,i,j,k,numffts
if (n1==1) return
ASSERT("fft1: dimension too small ",n1+2<=n1d);
call getindex(n1,index)

scale=n1
scale=1/scale

do k=1,n3

!   j=dfft_apply_grp('R','C','F',p(1,1,k),p(1,1,k),n2,1,&
!           fftdata(index)%trigs,1,n1d)
!   ASSERT("Error: fft1() status <> 0",j==0)


   do j=1,n2
      call dfft_apply('R','C','F',p(1,j,k),p(1,j,k),fftdata(index)%trigs,1)
!     move the last cosine mode into slot of first sine mode:
      p(2,j,k)=p(n1+1,j,k)
      do i=1,n1
         p(i,j,k)=p(i,j,k)*scale
      enddo
   enddo


enddo
end subroutine









end ! module mod_fft_interface


