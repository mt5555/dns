poisson(fhat,f)
!
!  solve laplacian(fhat) = f
!
implicit none
use params
! input
real*8 f(nxd,nyd,nzd)
! output
real*8 fhat(nxd,nyd,nzd)

!local
integer n1,n1d,n2,n2d,n3,n3d

n1=nx
n1d=nxd
n2=ny
n2d=nyd
n3=nz
n3d=nzd

call fft(f,fhat,n1,n1d,n2,n2d,n3,n3d)          
call transpose12(fhat,f,n1,n1d,n2,n2d,n3,n3d)  ! x,y,z -> y,x,z
call fft(f,fhat,n1,n1d,n2,n2d,n3,n3d)                        
call transpose13(fhat,f,n1,n1d,n2,n2d,n3,n3d)  ! y,x,z -> z,x,y
call fft(f,fhat,n1,n1d,n2,n2d,n3,n3d)

fhat /= (l**2 + m**2 + n**2)

call ifft(fhat,f,n1,n1d,n2,n2d,n3,n3d)
call transpose13(f,fhat,n1,n1d,n2,n2d,n3,n3d)         ! z,x,y -> y,x,z
call ifft(fhat,f,n1,n1d,n2,n2d,n3,n3d)
call transpose12(f,fhat,n1,n1d,n2,n2d,n3,n3d)         ! y,x,z -> x,y,z
call ifft(fhat,f,n1,n1d,n2,n2d,n3,n3d)

ASSERT(n1==nx)
ASSERT(n1d==nxd)
ASSERT(n2==ny)
ASSERT(n2d==nyd)
ASSERT(n3==nz)
ASSERT(n3d==nzd)

end
