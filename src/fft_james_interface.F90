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
!   fft_derivatives                  compute derivative along first dimension 
!                                    (input/ouput given in grid space)
! 
! Routines work on data of the form:  p(n1d,n2d,n3d)
! Size of the grid point data         p(1:n1,1:n2,1:n3)
! Size of fourier coefficients        p(1:n1,1:n2,1:n3)
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
implicit none
integer, parameter ::  num_fftsizes=3
real*8 :: pi2,pi2_squared

integer :: init=0
type fftdata_d
   real*8 :: ptrigs   ! used as a C pointer.  make sure at least 64bit
   integer :: size
end type
type(fftdata_d) :: fftdata(num_fftsizes)


private :: fftinit, getindex
contains 




subroutine fft_interface_init()

integer i
real*8  :: pi,one=1

do i=1,num_fftsizes
   fftdata(i)%size = 0	
enddo
pi=4*atan(one)
pi2=2*pi
pi2_squared=4*pi*pi
init=1
end subroutine


subroutine fft_get_mcord(mcord,n)
integer n,mcord(:)
integer i,m
do i=1,n
   m=i/2
   mcord(i)=m	
enddo
end subroutine







subroutine fftinit(n,index)
integer n,index
character*80 message

if (init==0) call abort("fft_james_interface.F90: call fft_interface_init to initialize first!");
if (n>1000000) call abort("fft_james_interface.F90: n>1 million")

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
      write(message_str,'(a,i10)') "fft_james_interface.F90:  Failed initializing an fft of size =",n1
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


integer index,k
if (n1==1) return
ASSERT("ifft1: dimension too small ",n1<=n1d);
call getindex(n1,index)

do k=1,n3
   call rfft_synthesis_m(p(1,1,k),fftdata(index)%ptrigs,n2,n1d)
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

integer index,i,j,k
if (n1==1) return
ASSERT("fft1: dimension too small ",n1<=n1d);
call getindex(n1,index)

scale=n1
scale=1/scale

do k=1,n3
   call rfft_analysis_m(p(1,1,k),fftdata(index)%ptrigs,n2,n1d)
   do j=1,n2   
   do i=1,n1
      p(i,j,k)=p(i,j,k)*scale
   enddo
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
ASSERT("fft99_interface.F90: numder>=1",numder>=1)

call fft1(px,n1,n1d,n2,n2d,n3,n3d)

if (numder>=2) then
   do k=1,n3
   do j=1,n2
      do i=1,n1-1
         m=i/2
         pxx(i,j,k) = -m*m * pi2_squared * px(i,j,k)
      enddo

      ! note: for i=n1, we are working with the cos(n1/2) mode
      ! but d/dx of this mode goes to sin(n1/2) = 0 on our grid, so we
      ! just take m=0 to make sure that dxx = dx dx
      pxx(n1,j,k)=0

   enddo
   enddo
   call ifft1(pxx,n1,n1d,n2,n2d,n3,n3d)
endif

   do k=1,n3
   do j=1,n2
      ! note: for i=2, m=0, we are actually working with the cos(n1/2) mode
      ! but d/dx of this mode goes to sin(n1/2) = 0, so just take m=0 

      px(1,j,k)=0              !m=0 cosine mode
      px(n1,j,k)=0             !m=n1/2 cosine mode
      do m = 1, n1/2-1
         i=2*m
         temp =  pi2* m * px(i,j,k)
         px(i,j,k) = -pi2 *m * px(i+1,j,k)
         px(i+1,j,k) = temp
      enddo
   enddo
   enddo
   call ifft1(px,n1,n1d,n2,n2d,n3,n3d)

end subroutine




subroutine divfree(u,p)
!
! make u divergence free
!    solve:  div(u) = laplacian(p)
!    then:   unew = u - grad(p)
!    
! 
!
use params
implicit none
real*8 u(nx,ny,nz,3)
real*8 p(nx,ny,nz)

!local variables
integer i,j,k
integer im,jm,km
real*8 :: xfac

ASSERT("divfree(): nslabx must be even ",mod(nslabx,2)==0)
ASSERT("divfree(): nslaby must be even ",mod(nslaby,2)==0)
ASSERT("divfree(): nslabz must be even ",(mod(nslabz,2)==0 .or. nslabz==1))

do i=1,3
   call fft3d(u(1,1,1,i),p)
enddo

   do k=nz1,nz2
      km=kmcord(k)
      if (km==g_nz/2) km=0
      do j=ny1,ny2
         jm=jmcord(j)
         if (jm==g_ny/2) jm=0
         do i=nx1,nx2
            im=imcord(i)
            if (im==g_nx/2) im=0

            ! compute the divergence
            p(i,j,k)=0
 
            if (im>0) then
            if (mod(i-nx1,2)==1) then
               p(i,j,k)=p(i,j,k) - im*u(i+1,j,k,1)
            else
               p(i,j,k)=p(i,j,k) + im*u(i-1,j,k,1)
            endif
            endif

            if (jm>0) then
            if (mod(j-ny1,2)==1) then
               p(i,j,k)=p(i,j,k) - jm*u(i,j+1,k,2)
            else
               p(i,j,k)=p(i,j,k) + jm*u(i,j-1,k,2)
            endif
            endif
            
            if (km>0) then
            if (mod(k-nz1,2)==1) then
               p(i,j,k)=p(i,j,k) - km*u(i,j,k+1,3)
            else
               p(i,j,k)=p(i,j,k) + km*u(i,j,k-1,3)
            endif
            endif

            ! compute laplacian inverse
            xfac= (im*im +km*km + jm*jm)
            if (xfac/=0) xfac = -1/xfac
            p(i,j,k)=xfac*p(i,j,k)


         enddo
      enddo
   enddo

   do k=nz1,nz2
      km=kmcord(k)
      if (km==g_nz/2) km=0
      do j=ny1,ny2
         jm=jmcord(j)
         if (jm==g_ny/2) jm=0
         do i=nx1,nx2
            im=imcord(i)
            if (im==g_nx/2) im=0

            ! compute u - dp/dx
            if (im>0) then
            if (mod(i-nx1,2)==1) then
               u(i,j,k,1)=u(i,j,k,1) + im*p(i+1,j,k)
            else
               u(i,j,k,1)=u(i,j,k,1) - im*p(i-1,j,k)
            endif
            endif

            if (jm>0) then
            if (mod(j-ny1,2)==1) then
               u(i,j,k,2)=u(i,j,k,2) + jm*p(i,j+1,k)
            else
               u(i,j,k,2)=u(i,j,k,2) - jm*p(i,j-1,k)
            endif
            endif

            if (km>0) then
            if (mod(k-nz1,2)==1) then
               u(i,j,k,3)=u(i,j,k,3) + km*p(i,j,k+1)
            else
               u(i,j,k,3)=u(i,j,k,3) - km*p(i,j,k-1)
            endif
            endif


         enddo
      enddo
   enddo

do i=1,3
   if (dealias) call fft_filter_dealias(u(1,1,1,i))
   call ifft3d(u(1,1,1,i),p)
enddo
end subroutine





end ! module mod_fft_interface
