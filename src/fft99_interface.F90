#include "macros.h"
#if 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! our wrapper for ECMWF FFT99.
!
! provides public interfaces:
!   fft_interface_init               call this before using any other routines
!   fft
!   ifft
!   fft_derivatives
!   fft_laplace_inverse
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif


module fft_interface

implicit none
integer, parameter ::  num_fftsizes=3
real*8 :: pi2,pi2_squared

integer :: init=0
type fftdata_d
   real*8,dimension(:),pointer :: trigs
   integer :: ifax(13)
   integer :: size
end type
type(fftdata_d) :: fftdata(num_fftsizes)


integer, parameter ::  fftblocks=100   ! do fft's in blocks of size fftblocks
                                       ! set very large to disable
private :: fftinit, getindex
contains 




subroutine fft_interface_init()

integer i
real*8  :: pi,one=1

do i=1,num_fftsizes
   fftdata(i).size=0	
enddo
pi=4*atan(one)
pi2=2*pi
pi2_squared=4*pi*pi
init=1

end subroutine


subroutine fftinit(n,index)
integer n,index
character*80 message

if (init==0) call abort("fft99_interface.F90: call fft_interface_init to initialize first!");
if (n>1000000) call abort("fft99_interface.F90: n>1 million")

fftdata(index).size=n
allocate(fftdata(index).trigs(3*n/2+1))
call set99(fftdata(index).trigs,fftdata(index).ifax,n)


write(message,'(a,i6)') 'Initializing FFT of size n=',n
call print_message(message)
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

   if (fftdata(i).size==0) then
      call fftinit(n1,i)      
      exit 
   endif
   if (n1==fftdata(i).size) exit 
enddo
index=i
end subroutine






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute 3D in-place iFFT of p
! FFT taken along first direction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ifft(p,n1,n1d,n2,n2d,n3,n3d)
real*8 p(n1d,n2d,n3d)
real*8 w(min(fftblocks,n2)*(n1+1))
integer n1,n1d,n2,n2d,n3,n3d
character*80 message_str

integer index,jj,j,k,numffts

if (n1==1) return
call getindex(n1,index)


j=0  ! j=number of fft's computed for each k
do k=1,n3
   j=0  ! j=number of fft's computed for each k
   do while (j<n2)
      numffts=min(fftblocks,n2-j)	
!      do jj=j+1,j+numffts
!         p(2,jj,k)=0
!         p(n1+2,jj,k)=0
!      enddo
      call fft991(p(1,j+1,k),w,fftdata(index).trigs,fftdata(index).ifax,1,n1d,n1,numffts,1)
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
subroutine fft(p,n1,n1d,n2,n2d,n3,n3d)
integer n1,n1d,n2,n2d,n3,n3d
real*8 p(n1d,n2d,n3d)
real*8 :: w(min(fftblocks,n2)*(n1+1)) 

integer index,jj,j,k,numffts
if (n1==1) return
call getindex(n1,index)

do k=1,n3
   j=0  ! j=number of fft's computed for each k
   do while (j<n2)
      numffts=min(fftblocks,n2-j)	
!      do jj=j+1,j+numffts
!         p(n1+1,jj,k)=0
!         p(n1+2,jj,k)=0
!      enddo

      call fft991(p(1,j+1,k),w,fftdata(index).trigs,fftdata(index).ifax,1,n1d,n1,numffts,-1)
      j=j+numffts
   enddo
enddo


end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Compute FFT derivates.
!
!  input: px 
!  output:
!     if numder=1   return d/dx along first direction in px
!                   (and pxx is not accessed) 
!     if numder=2   return d2/dx2 along first direction in pxx
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_derivatives(px,pxx,numder,n1,n1d,n2,n2d,n3,n3d)
real*8 px(n1d,n2d,n3d)
real*8 pxx(n1d,n2d,n3d)
integer numder,n1,n1d,n2,n2d,n3,n3d

integer i,j,k,m
real*8 temp
ASSERT("fft99_interface.F90: numder<=2",numder<=2)

call fft(px,n1,n1d,n2,n2d,n3,n3d)

if (numder>=2) then
   do k=1,n3
   do j=1,n2
      do i=1,n1
         m=(i-1)/2
         pxx(i,j,k) = -m*m * pi2_squared * px(i,j,k)
      enddo
      ! tweak do that dxx = (dx)(dx)
      ! this is because we have a cos((n1/2) 2pi x) mode, but no sine!
      pxx(n1+1,j,k) = 0
      pxx(n1+2,j,k) = 0
   enddo
   enddo
   call ifft(pxx,n1,n1d,n2,n2d,n3,n3d)
endif

if (numder>=1) then
   do k=1,n3
   do j=1,n2
      do m = 0, n1/2-1
         i = 2*m+1
         temp =  pi2* m * px(i,j,k)
         px(i,j,k) = -pi2 *m * px(i+1,j,k)
         px(i+1,j,k) = temp
      enddo
      px(n1+1,j,k)=0
      px(n1+2,j,k)=0
   enddo
   enddo
   call ifft(px,n1,n1d,n2,n2d,n3,n3d)
endif
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Solve  [alpha + beta*Laplacian] p = rhs
!
! on input,  p = fourier coefficients of rhs
! on output, p = fourier coefficients of solution
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_laplace_inverse(p,n1,n1d,n2,n2d,n3,n3d,alpha,beta)
real*8 p(n1d,n2d,n3d)
real*8 alpha,beta
integer n1,n1d,n2,n2d,n3,n3d

integer i,j,k,im,jm,km
real*8 xfac

   do k=1,n3+2
      km=(k-1)/2
      do j=1,n2+2
         jm=(j-1)/2
         do i=1,n1+2
            im=(i-1)/2
            xfac= alpha + beta*(-im*im -km*km - jm*jm)*pi2_squared      
            if (xfac<>0) xfac = 1/xfac
            p(i,j,k)=p(i,j,k)*xfac
         enddo
      enddo
   enddo

end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Filter out the highest cosine mode
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_filter(p,n1,n1d,n2,n2d,n3,n3d)
real*8 p(n1d,n2d,n3d)
real*8 alpha,beta
integer n1,n1d,n2,n2d,n3,n3d

integer i,j,k,im,jm,km
real*8 xfac

   do k=1,n3+2
      km=(k-1)/2
      do j=1,n2+2
         jm=(j-1)/2
         do i=1,n1+2
            im=(i-1)/2
            if (km == n3/2) p(i,j,k)=0
            if (jm == n2/2) p(i,j,k)=0
            if (im == n1/2) p(i,j,k)=0
         enddo
      enddo
   enddo

end subroutine




end ! module mod_fft_interface
