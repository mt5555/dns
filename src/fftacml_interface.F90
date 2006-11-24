#include "macros.h"
#if 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! our wrapper for AMD ACML FFT 
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
ACML representation:

f=[fhat(1) + 2 fhat(2*m+1) cos(m*2pi*x) - 2*fhat(2m+2) sin(m*2pi*x) ] / sqrt(n)

X(i)   = real part of Z(i) i=0.. N/2
X(N-i) = imag part of Z(i) i=1.. (N-1)/2

 grid space data:    1 2 3 4 5 6 7 8 * *

 i                   1 2 3 4 5 6 7 8
 fourier space:      0 1 2 3 4 3 2 1 * *
 rearranged:         0 4 1 1 2 2 3 3 * *

xout(1)=x(1)
xout(2)=x(n/2+1)
do i=2,n/2
   xout(2*i-1) = x(i) 
   xout(2*i)   = x(n-i+2)
enddo


NOTE: ACML has multiple real FFT, but it doesnt have an "mode"
argument and it requires no padding (i.e. real*8 X(N,M))

integer :: info,n,mode
real*8 :: X(n)
real*8 :: comm(3*N+100)
mode=0:  default initialization
mode=1:  real transform is perfmed
mode=2:  init + transform
mode=100:  expensive init


   CALL DZFFT(0,N,X,COMM,INFO)  ! init
   CALL DZFFT(1,N,X,COMM,INFO)  ! forward transform
   DO 10 I = N/2+2, N           ! prep for inverse transform
      X(I) = -X(I)
10 CONTINUE
   CALL ZDFFT(2,N,X,COMM,INFO)  ! inverse transform


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
use params
implicit none
integer, parameter ::  num_fftsizes=3

integer :: init=0
type fftdata_d
   real*8,dimension(:),pointer :: trigs1
   real*8,dimension(:),pointer :: trigs2
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
!
!  i=1   0 cosine mode             mcord=0
!  i=2   n/2 cosine mode           mcord=  n/2
!  i=3   1 cosine mode             mcord =  1
!  i=4   1 sine mode               mcord = -1
!  i=5   2 cosine mode             mcord =  2
!  i=6   2 sine mode               mcord=  -2
!  etc...
!
integer n,mcord(:)
integer i
do i=1,n
   mcord(i)=(i-1)/2	
   if (mod(i,2)==0) mcord(i)=-mcord(i)
   if (i==2) mcord(i)=n/2
enddo
end subroutine







subroutine fftinit(n,index)
integer n,index,info
character(len=80) message

if (init==0) call abortdns("fftacml_interface.F90: call fft_interface_init to initialize first!")
if (n>1000000) call abortdns("fftacml_interface.F90: n>1 million")
if (n<0) call abortdns("Error: invalid value of n for fft")

fftdata(index)%size=n
allocate(fftdata(index)%trigs1(3*n + 100))
allocate(fftdata(index)%trigs2(3*n + 100))

write(message,'(a,i6)') 'Initializing ACML DZFFT of size n=',n
call print_message(message)

call DZFFTM(100,n,dummy,fftdata(index)%trigs1,info)

write(message,'(a,i6)') 'Initializing ACML ZDFFT of size n=',n
call print_message(message)

call ZDFFTM(100,n,dummy,fftdata(index)%trigs2,info)

write(message,'(a)') 'Initializing ACML done'
call print_message(message)


end subroutine




subroutine getindex(n1,index)
integer n1,index

character(len=80) message_str
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
integer n1,n1d,n2,n2d,n3,n3d
real*8 p(n1d,n2d,n3d)
real*8 w(n2*(n1+1))

real*8 :: scale=1,tmx1,tmx2
character(len=80) message_str
integer index,j,k

call wallclock(tmx1)
if (tims(19)==0) ncalls(19)=0  ! timer was reset, so reset counter too

if (n1==1) return
ASSERT("ifft1: dimension too small ",n1<=n1d)
call getindex(n1,index)



j=0  ! j=number of fft's computed for each k
do k=1,n3
   do j=1,n2

      ! move the last cosine mode back into correct location:
      p(n1+1,j,k)=p(2,j,k)
      !p(2,j,k)=0             ! not needed?
      !p(n1+2,j,k)=0          ! not needed?
   enddo

   call fft991(p(1,1,k),w,fftdata(index)%trigs,fftdata(index)%ifax,1,n1d,n1,n2,1)
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
real*8 :: scale
real*8 :: w(n2*(n1+1)) 
real*8 :: tmx1,tmx2
integer index,j,k

call wallclock(tmx1)
if (tims(18)==0) ncalls(18)=0  ! timer was reset, so reset counter too


if (n1==1) return
ASSERT("fft1: dimension too small ",n1<=n1d)
call getindex(n1,index)

scale=n1
scale=1/scale

do k=1,n3
   !   do j=1,n2
   !         p(n1+1,jj,k)=0
   !         p(n1+2,jj,k)=0
   !   enddo
   call fft991(p(1,1,k),w,fftdata(index)%trigs,fftdata(index)%ifax,1,n1d,n1,n2,-1)
   !     move the last cosine mode into slot of first sine mode:
   do j=1,n2
      p(2,j,k)=p(n1+1,j,k)
   enddo
   
enddo
call wallclock(tmx2) 
tims(18)=tims(18)+(tmx2-tmx1)          
ncalls(18)=ncalls(18)+1
end subroutine




end ! module mod_fft_interface


