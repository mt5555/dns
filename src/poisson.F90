#include "macros.h"

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
call der(u(1,1,1,i),d1,dummy,work,1,i)
p = d1
i=2
call der(u(1,1,1,i),d1,dummy,work,1,i)
p = p + d1
i=3
call der(u(1,1,1,i),d1,dummy,work,1,i)
p = p + d1

divu=p

! solve laplacian(p)=div(u)
call poisson(p,work,alpha,beta)


! compute u=u-grad(p)
do i=1,3
   call der(p,d1,dummy,work,1,i)
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

print *,'first fft dims: ',n1,n1d,n2,n2d,n3,n3d
call fft(f,n1,n1d,n2,n2d,n3,n3d)     
n1=n1+2

call transpose12(f,work,1,n1,n1d,n2,n2d,n3,n3d)  ! x,y,z -> y,x,z
call fft(work,n1,n1d,n2,n2d,n3,n3d)
n1=n1+2

call transpose13(work,f,1,n1,n1d,n2,n2d,n3,n3d)  ! y,x,z -> z,x,y
call fft(f,n1,n1d,n2,n2d,n3,n3d)
n1=n1+2


! solve [alpha + beta*Laplacian] p = f.  f overwritten with output p
call fft_laplace_inverse(f,n1-2,n1d,n2-2,n2d,n3-2,n3d,alpha,beta)


n1=n1-2
call ifft(f,n1,n1d,n2,n2d,n3,n3d)
call transpose13(f,work,1,n1,n1d,n2,n2d,n3,n3d)         ! z,x,y -> y,x,z

n1=n1-2
call ifft(work,n1,n1d,n2,n2d,n3,n3d)
call transpose12(work,f,1,n1,n1d,n2,n2d,n3,n3d)         ! y,x,z -> x,y,z

n1=n1-2
call ifft(f,n1,n1d,n2,n2d,n3,n3d)



ASSERT("poisson.F90: n1<>nx2",n1==nx2)
ASSERT("poisson.F90: n1d<>nx",n1d==nx)
ASSERT("poisson.F90: n2<>ny2",n2==ny2)
ASSERT("poisson.F90: n2d<>ny",n2d==ny)
ASSERT("poisson.F90: n3<>nz2",n3==nz2)
ASSERT("poisson.F90: n3d<>nz",n3d==nz)



end



subroutine filter(f,work)
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


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k

n1=nx2
n1d=nx
n2=ny2
n2d=ny
n3=nz2
n3d=nz

print *,'first fft dims: ',n1,n1d,n2,n2d,n3,n3d
call fft(f,n1,n1d,n2,n2d,n3,n3d)     
n1=n1+2

call transpose12(f,work,1,n1,n1d,n2,n2d,n3,n3d)  ! x,y,z -> y,x,z
call fft(work,n1,n1d,n2,n2d,n3,n3d)
n1=n1+2

call transpose13(work,f,1,n1,n1d,n2,n2d,n3,n3d)  ! y,x,z -> z,x,y
call fft(f,n1,n1d,n2,n2d,n3,n3d)
n1=n1+2


! solve [alpha + beta*Laplacian] p = f.  f overwritten with output p
call fft_filter(f,n1-2,n1d,n2-2,n2d,n3-2,n3d)


n1=n1-2
call ifft(f,n1,n1d,n2,n2d,n3,n3d)
call transpose13(f,work,1,n1,n1d,n2,n2d,n3,n3d)         ! z,x,y -> y,x,z

n1=n1-2
call ifft(work,n1,n1d,n2,n2d,n3,n3d)
call transpose12(work,f,1,n1,n1d,n2,n2d,n3,n3d)         ! y,x,z -> x,y,z

n1=n1-2
call ifft(f,n1,n1d,n2,n2d,n3,n3d)



ASSERT("poisson.F90: n1<>nx2",n1==nx2)
ASSERT("poisson.F90: n1d<>nx",n1d==nx)
ASSERT("poisson.F90: n2<>ny2",n2==ny2)
ASSERT("poisson.F90: n2d<>ny",n2d==ny)
ASSERT("poisson.F90: n3<>nz2",n3==nz2)
ASSERT("poisson.F90: n3d<>nz",n3d==nz)



end

