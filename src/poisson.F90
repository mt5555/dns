poisson(f,work)
!
!  solve laplacian(p) = f
!  overrite f with the solution p 
!
implicit none
use params
! input/output
real*8 f(nxd,nyd,nzd)
! work array
real*8 work(nxd,nyd,nzd)

!local
integer n1,n1d,n2,n2d,n3,n3d

n1=nx
n1d=nxd
n2=ny
n2d=nyd
n3=nz
n3d=nzd

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


ASSERT(n1==nx)
ASSERT(n1d==nxd)
ASSERT(n2==ny)
ASSERT(n2d==nyd)
ASSERT(n3==nz)
ASSERT(n3d==nzd)

end
