#include "macros.h"

subroutine poisson(f,work)
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
real*8 :: alpha = 0
real*8 :: beta = 1
!local
integer n1,n1d,n2,n2d,n3,n3d

n1=nx2
n1d=nx
n2=ny2
n2d=ny
n3=nz2
n3d=nz

call fft(f,n1,n1d,n2,n2d,n3,n3d)     
n1=n1+2

print *,'before: ',n1,n1d,n2,n2d,n3,n3d
call transpose12(f,work,1,n1,n1d,n2,n2d,n3,n3d)  ! x,y,z -> y,x,z
print *,'after: ',n1,n1d,n2,n2d,n3,n3d
call fft(work,n1,n1d,n2,n2d,n3,n3d)
n1=n1+2
#if 0


call transpose13(work,f,1,n1,n1d,n2,n2d,n3,n3d)  ! y,x,z -> z,x,y
call fft(f,n1,n1d,n2,n2d,n3,n3d)
n1=n1+2

! solve [alpha + beta*Laplacian] p = f.  f overwritten with output p
!call fft_laplace_inverse(f,n1-2,n1d,n2-2,n2d,n3-2,n3d,alpha,beta)

n1=n1-2
call ifft(f,n1,n1d,n2,n2d,n3,n3d)
call transpose13(f,work,1,n1,n1d,n2,n2d,n3,n3d)         ! z,x,y -> y,x,z
#endif

n1=n1-2
call ifft(work,n1,n1d,n2,n2d,n3,n3d)
print *,'before: ',n1,n1d,n2,n2d,n3,n3d
call transpose12(work,f,1,n1,n1d,n2,n2d,n3,n3d)         ! y,x,z -> x,y,z
print *,'after: ',n1,n1d,n2,n2d,n3,n3d


n1=n1-2
call ifft(f,n1,n1d,n2,n2d,n3,n3d)



ASSERT("poisson.F90: n1<>nx2",n1==nx2)
ASSERT("poisson.F90: n1d<>nx",n1d==nx)
ASSERT("poisson.F90: n2<>ny2",n2==ny2)
ASSERT("poisson.F90: n2d<>ny",n2d==ny)
ASSERT("poisson.F90: n3<>nz2",n3==nz2)
ASSERT("poisson.F90: n3d<>nz",n3d==nz)

end
