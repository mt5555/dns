

der(p,px,pxx,n,index)
!
!  compute derivative along index index=1,2 or 3
!  n = number of derivatives to compute
!

implicit none
use params

!input:
integer n,index
real*8 p(nxd,nyd,nzd)
!output:
real*8 px(nxd,nyd,nzd)
real*8 pxx(nxd,nyd,nzd)

!local
real*8 pt(nxd*nyd*nzd)
real*8 tmp(nxd*nyd*nzd)

integer n1,n1d,n2,n2d,n3,n3d

n1=nx
n1d=nxd
n2=ny
n2d=nyd
n3=nz
n3d=nzd


if (index==1) then
   fft_derivatives(p,px,pxx,n1,n1d,n2,n2d,n3,n3d)
else if (index==2) then
   transpose12(p,ptmp,n1,n1d,n2,n2d,n3,n3d)
   ! tmp 1st derivative, px=2nd derivative
   fft_derivatives(pt,tmp,px,n,n1,n1d,n2,n2d,n3,n3d)
   if (n==2) transpose12(px,pxx,n1,n1d,n2,n2d,n3,n3d)
   transpose12(pxt,px,n1,n1d,n2,n2d,n3,n3d)
else if (index==3) then
   transpose13(p,pt,n1,n1d,n2,n2d,n3,n3d0
   fft_derivatives(pt,tmp,px,1,n1,n1d,n2,n2d,n3,n3d)
   if (n==2) transpose13(px,pxx,n1,n1d,n2,n2d,n3,n3d)
   transpose13(tmp,px,n1,n1d,n2,n2d,n3,n3d)
endif

ASSERT(n1==nx)
ASSERT(n1d==nxd)
ASSERT(n2==ny)
ASSERT(n2d==nyd)
ASSERT(n3==nz)
ASSERT(n3d==nzd)

end


