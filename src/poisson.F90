#include "macros.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  compute derivative along index index=1,2 or 3
!  numder = 1  compute p_x, return in px.   (pxx is not accessed)
!  numder = 2  compute p_xx, return in pxx.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine der(p,px,pxx,pt,numder,index)
use params
use fft_interface
implicit none

!input:
integer numder,index
real*8 p(nx,ny,nz)    ! original data
real*8 pt(nx,ny,nz)   ! work array

!output:
real*8 pxx(nx,ny,nz)
real*8 px(nx,ny,nz)

integer n1,n1d,n2,n2d,n3,n3d

n1=nx2
n1d=nx
n2=ny2
n2d=ny
n3=nz2
n3d=nz


if (index==1) then

   px=p
   call fft_derivatives(px,pxx,numder,n1,n1d,n2,n2d,n3,n3d)

else if (index==2) then

   call transpose12(p,pt,SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)
   ! 1st derivative returned in pt, 2nd derivative returned in px
   call fft_derivatives(pt,px,numder,n1,n1d,n2,n2d,n3,n3d)
   if (numder==2) then
      call transpose12(px,pxx,DONT_SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)
   endif
   call transpose12(pt,px,SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)

else if (index==3) then
   call transpose13(p,pt,SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)
   call fft_derivatives(pt,px,numder,n1,n1d,n2,n2d,n3,n3d)
   if (numder==2) then
      call transpose13(px,pxx,DONT_SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)
   endif
   call transpose13(pt,px,SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)

endif

ASSERT("derivatives.F90 transpose error n1",n1==nx2)
ASSERT("derivatives.F90 transpose error n1d",n1d==nx)
ASSERT("derivatives.F90 transpose error n2",n2==ny2)
ASSERT("derivatives.F90 transpose error n2d",n2d==ny)
ASSERT("derivatives.F90 transpose error n3",n3==nz2)
ASSERT("derivatives.F90 transpose error n3d",n3d==nz)

end subroutine



subroutine divfree(u)
!
! make u divergence free
!    solve:  div(u) = laplacian(p)
!    then:   unew = u - grad(p)
!    
! 
!
use params
use fft_interface
implicit none
real*8 u(nx,ny,nz,3)
real*8 d1(nx,ny,nz)
real*8 work(nx,ny,nz)
real*8 p(nx,ny,nz)
real*8 :: dummy
real*8 :: alpha=0
real*8 :: beta=1

integer i,j,k
real*8 divu(nx,ny,nz)
real*8 px(nx,ny,nz)
real*8 pxx(nx,ny,nz)
real*8 lap(nx,ny,nz)
real*8 lap2(nx,ny,nz)

! compute p = div(u)
i=1
call der(u(1,1,1,i),d1,dummy,work,DX_ONLY,i)
p = d1
i=2
call der(u(1,1,1,i),d1,dummy,work,DX_ONLY,i)
p = p + d1
i=3
call der(u(1,1,1,i),d1,dummy,work,DX_ONLY,i)
p = p + d1

divu=p

! solve laplacian(p)=div(u)
call poisson(p,work,alpha,beta)


! compute u=u-grad(p)
do i=1,3
   call der(p,d1,dummy,work,DX_ONLY,i)
   u(:,:,:,i) = u(:,:,:,i) - d1
enddo


end subroutine








subroutine poisson(f,work,alpha,beta)
!
!  solve laplacian(p) = f
!  input:  f 
!  ouput:  f   will be overwritten with the solution p
!
use params
use fft_interface
implicit none
real*8 f(nx,ny,nz)    ! input/output
real*8 work(nx,ny,nz) ! work array
real*8 :: alpha
real*8 :: beta


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k

n1=nx2
n1d=nx
n2=ny2
n2d=ny
n3=nz2
n3d=nz

call fft1(f,n1,n1d,n2,n2d,n3,n3d)     
call transpose12(f,work,SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)  ! x,y,z -> y,x,z
call fft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose13(work,f,SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)  ! y,x,z -> z,x,y
call fft1(f,n1,n1d,n2,n2d,n3,n3d)

! solve [alpha + beta*Laplacian] p = f.  f overwritten with output  p
call fft_laplace_inverse(f,n1,n1d,n2,n2d,n3,n3d,alpha,beta)

call ifft1(f,n1,n1d,n2,n2d,n3,n3d)
call transpose13(f,work,SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)         ! z,x,y -> y,x,z
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose12(work,f,SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)         ! y,x,z -> x,y,z
call ifft1(f,n1,n1d,n2,n2d,n3,n3d)



ASSERT("poisson.F90 poisson: n1<>nx2",n1==nx2)
ASSERT("poisson.F90 poisson: n1d<>nx",n1d==nx)
ASSERT("poisson.F90 poisson: n2<>ny2",n2==ny2)
ASSERT("poisson.F90 poisson: n2d<>ny",n2d==ny)
ASSERT("poisson.F90 poisson: n3<>nz2",n3==nz2)
ASSERT("poisson.F90 poisson: n3d<>nz",n3d==nz)



end




subroutine fft3d(f,work,n1,n2,n3)
!
!  compute the spectrum, ouput in f
!  also output the dimensions of the fourier coefficints, n1,n2,n3
!  (since fft99 increases them by 2)
!
use params
use fft_interface
implicit none
real*8 f(nx,ny,nz)    ! input/output
real*8 work(nx,ny,nz) ! work array
integer n1,n1d,n2,n2d,n3,n3d

n1=nx2
n1d=nx
n2=ny2
n2d=ny
n3=nz2
n3d=nz

call fft1(f,n1,n1d,n2,n2d,n3,n3d)     
call transpose12(f,work,SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)  ! x,y,z -> y,x,z
call fft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose12(work,f,SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)  
call transpose13(f,work,SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)  ! x,y,z -> z,y,x
call fft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose13(work,f,SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)  


! dimension has grown by 2 because of fft99:
ASSERT("poisson.F90: fft3d: n1<>nx2",n1-2==nx2)
ASSERT("poisson.F90: fft3d: n1d<>nx",n1d==nx)
ASSERT("poisson.F90: fft3d: n2<>ny2",n2-2==ny2)
ASSERT("poisson.F90: fft3d: n2d<>ny",n2d==ny)
if (n3>1) then
   ASSERT("poisson.F90: fft3d: n3<>nz2",n3-2==nz2)
else
   ASSERT("poisson.F90: fft3d: n3<>nz2",n3==nz2)
endif
ASSERT("poisson.F90: fft3d: n3d<>nz",n3d==nz)



end




subroutine ifft3d(f,work,n1,n2,n3)
!
!  compute inverse fft 3d of f, return in f
!  n1,n2,n3 = size of fft coefficient array (which could be different
!  then grid point array)
!
use params
use fft_interface
implicit none
real*8 f(nx,ny,nz)    ! input/output
real*8 work(nx,ny,nz) ! work array


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k


n1d=nx
n2d=ny
n3d=nz


call transpose13(f,work,SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)         ! x,y,z -> z,y,x
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose13(work,f,SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)       
call transpose12(f,work,SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)        ! x,y,z -> y,x,z
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose12(work,f,SWAP_INDEX,n1,n1d,n2,n2d,n3,n3d)       
call ifft1(f,n1,n1d,n2,n2d,n3,n3d)


ASSERT("poisson.F90: ifft3d: n1<>nx2",n1==nx2)
ASSERT("poisson.F90: ifft3d: n1d<>nx",n1d==nx)
ASSERT("poisson.F90: ifft3d: n2<>ny2",n2==ny2)
ASSERT("poisson.F90: ifft3d: n2d<>ny",n2d==ny)
ASSERT("poisson.F90: ifft3d: n3<>nz2",n3==nz2)
ASSERT("poisson.F90: ifft3d: n3d<>nz",n3d==nz)



end

