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

!local
integer n1,n1d,n2,n2d,n3,n3d

n1=nx2
n1d=nx
n2=ny2
n2d=ny
n3=nz2
n3d=nz

call fft(f,n1,n1d,n2,n2d,n3,n3d)          
call transpose12(f,work,n1,n1d,n2,n2d,n3,n3d)  ! x,y,z -> y,x,z
call fft(work,n1,n1d,n2,n2d,n3,n3d)                        
call transpose13(work,f,n1,n1d,n2,n2d,n3,n3d)  ! y,x,z -> z,x,y
call fft(f,n1,n1d,n2,n2d,n3,n3d)

call fft_laplace_inverse3d(f,n1,n1d,n2,n2d,n3,n3d)

call ifft(f,n1,n1d,n2,n2d,n3,n3d)
call transpose13(f,work,n1,n1d,n2,n2d,n3,n3d)         ! z,x,y -> y,x,z
call ifft(work,n1,n1d,n2,n2d,n3,n3d)
call transpose12(work,f,n1,n1d,n2,n2d,n3,n3d)         ! y,x,z -> x,y,z
call ifft(f,n1,n1d,n2,n2d,n3,n3d)


ASSERT("poisson.F90",n1==nx2)
ASSERT("poisson.F90",n1d==nx)
ASSERT("poisson.F90",n2==ny2)
ASSERT("poisson.F90",n2d==ny)
ASSERT("poisson.F90",n3==nz2)
ASSERT("poisson.F90",n3d==nz)

end
