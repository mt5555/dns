#include "macros.h"
#if 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! our wrapper for ECMWF FFT99 and SGI FFT
! (they both have the same output format)
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
   CPOINTER,dimension(:),pointer :: trigs
   integer :: ifax(13)
   integer :: size
end type
type(fftdata_d) :: fftdata(num_fftsizes)


integer, parameter ::  fftblocks=2000 ! do fft's in blocks of size fftblocks
                                      ! set very large to disable blocking and do 
                                      ! all fft's at once
private :: fftinit, getindex
contains 




subroutine fft_interface_init()

integer i
real*8 :: one=1

init=1
do i=1,num_fftsizes
   fftdata(i)%size = 0	
enddo
end subroutine




subroutine fft_get_mcord(mcord,n)
integer n,mcord(:)
integer i,m
do i=1,n
   m=(i-1)/2
   if (i==2) m=n/2   ! last cosine mode is stored at i=2
   mcord(i)=m	
enddo
end subroutine







subroutine fftinit(n,index)
integer n,index
character*80 message

if (init==0) call abort("fft99_interface.F90: call fft_interface_init to initialize first!");
if (n>1000000) call abort("fft99_interface.F90: n>1 million")

fftdata(index)%size=n
#ifdef USE_SGIFFT
allocate(fftdata(index)%trigs(n+15))
#else
allocate(fftdata(index)%trigs(3*n/2+1))
#endif

#ifdef USE_SGIFFT
write(message,'(a,i6)') 'Initializing SGI FFT of size n=',n
#else
write(message,'(a,i6)') 'Initializing FFT99 of size n=',n
#endif
call print_message(message)

#ifdef USE_SGIFFT
CALL DZFFTM (0, n, 1, 0, 0, 0, 0, 0, fftdata(index)%trigs, 0,0)
#else
call set99(fftdata(index)%trigs,fftdata(index)%ifax,n)
#endif
if (n<0) call abort("Error; invalid value of n for fft");

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
      call abort(message_str)
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
#ifdef USE_SGIFFT
real*8 w(n1+2)
#else
real*8 w(min(fftblocks,n2)*(n1+1))
#endif

real*8 :: scale=1
character*80 message_str

integer index,jj,j,k,numffts
if (n1==1) return
ASSERT("ifft1: dimension too small ",n1+2<=n1d);
call getindex(n1,index)



j=0  ! j=number of fft's computed for each k
do k=1,n3
   j=0  ! j=number of fft's computed for each k
   do while (j<n2)
      numffts=min(fftblocks,n2-j)	

!     move the last cosine mode back into correct location:
      do  jj=j+1,j+numffts
	p(n1+1,jj,k)=p(2,jj,k)
        !p(2,jj,k)=0             ! not needed?
        !p(n1+2,jj,k)=0          ! not needed?
      enddo     

#ifdef USE_SGIFFT
      CALL ZDFFTM (1, n1, numffts, scale, p(1,j+1,k), n1d/2, p(1,j+1,k), n1d,fftdata(index)%trigs, w,0)
#else
      call fft991(p(1,j+1,k),w,fftdata(index)%trigs,fftdata(index)%ifax,1,n1d,n1,numffts,1)
#endif

      j=j+numffts
   enddo
enddo

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
#ifdef USE_SGIFFT
real*8 w(n1+2)
#else
real*8 :: w(min(fftblocks,n2)*(n1+1)) 
#endif

integer index,jj,j,k,numffts
if (n1==1) return
ASSERT("fft1: dimension too small ",n1+2<=n1d);
call getindex(n1,index)

scale=n1
scale=1/scale

do k=1,n3
   j=0  ! j=number of fft's computed for each k
   do while (j<n2)
      numffts=min(fftblocks,n2-j)	
!      do jj=j+1,j+numffts
!         p(n1+1,jj,k)=0
!         p(n1+2,jj,k)=0
!      enddo

#ifdef USE_SGIFFT 
      CALL DZFFTM (-1, n1, numffts, scale, p(1,j+1,k), n1d, p(1,j+1,k), n1d/2,fftdata(index)%trigs, w,0)
#else
      call fft991(p(1,j+1,k),w,fftdata(index)%trigs,fftdata(index)%ifax,1,n1d,n1,numffts,-1)
#endif

!     move the last cosine mode into slot of first sine mode:
      do  jj=j+1,j+numffts
	p(2,jj,k)=p(n1+1,jj,k)
      enddo     


      j=j+numffts
   enddo
enddo
end subroutine









end ! module mod_fft_interface


