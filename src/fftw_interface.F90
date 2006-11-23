#include "macros.h"
#if 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! our wrapper for ECMWF FFT99 
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



fftw_plan fftw_plan_many_dft_r2c(int rank, const int *n, int howmany,
                                      double *in, const int *inembed,
                                      int istride, int idist,
                                      fftw_complex *out, const int *onembed,
                                      int ostride, int odist,
                                      unsigned flags);
fftw_plan fftw_plan_many_dft_c2r(int rank, const int *n, int howmany,
                                      fftw_complex *in, const int *inembed,
                                      int istride, int idist,
                                      double *out, const int *onembed,
                                      int ostride, int odist,
                                      unsigned flags);



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif


module fft_interface
use params
implicit none
include 'fftw3.f'

contains 




subroutine fft_interface_init()
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





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute 3D in-place iFFT of p
! FFT taken along first direction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ifft1(p,n1,n1d,n2,n2d,n3,n3d)
integer n1,n1d,n2,n2d,n3,n3d
real*8 p(n1d,n2d,n3d)

real*8 :: tmx1,tmx2
character(len=80) message_str
integer index,j,k
CPOINTER :: plan

call wallclock(tmx1)
if (tims(19)==0) ncalls(19)=0  ! timer was reset, so reset counter too

if (n1==1) return
ASSERT("ifft1: dimension too small ",n1+2<=n1d)
ASSERT("ifft1: dimension too small ",n2==n2d)

! number of FFTs to compute:  n2*n3
! stride:  n1d

do k=1,n3
   do j=1,n2
      ! move the last cosine mode back into correct location:
      p(n1+1,j,k)=p(2,j,k)
      p(2,j,k)=0             ! not needed?
      p(n1+2,j,k)=0          ! not needed?
   enddo
enddo

! FFTW_MEASURE
#if 1
do k=1,n3
   do j=1,n2
      call dfftw_plan_dft_c2r_1d(plan, n1,p(1,j,k),p(1,j,k),FFTW_ESTIMATE)
      call dfftw_execute(plan)
      call dfftw_destroy_plan(plan)
   enddo
enddo
#else
call dfftw_plan_many_dft_c2r(plan, 1, n1,n2*n3,p, n1d, &
                       1, n1d, p, n1d, 1, n1d, FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
#endif

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
real*8 :: p(n1d,n2d,n3d)
real*8 :: tmx1,tmx2
integer index,j,k
CPOINTER :: plan

call wallclock(tmx1)
if (tims(18)==0) ncalls(18)=0  ! timer was reset, so reset counter too


if (n1==1) return
ASSERT("fft1: dimension too small ",n1+2<=n1d)
ASSERT("fft1: dimension too small ",n2==n2d)

! number of FFTs to compute:  n2*n3
! stride:  n1d
#if 1
do k=1,n3
   do j=1,n2
      call dfftw_plan_dft_r2c_1d(plan, n1,p(1,j,k),p(1,j,k),FFTW_ESTIMATE)
      call dfftw_execute(plan)
      call dfftw_destroy_plan(plan)
   enddo
enddo


#else
call dfftw_plan_many_dft_r2c(plan, 1, n1,n2*n3,p, n1d, &
                       1, n1d, p, n1d, 1, n1d, FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
#endif


p=p/n1

do k=1,n3
   !   do j=1,n2
   !         p(n1+1,jj,k)=0
   !         p(n1+2,jj,k)=0
   !   enddo
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


