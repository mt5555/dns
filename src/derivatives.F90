der(p,px,pxx,pt,n,index)
!
!  compute derivative along index index=1,2 or 3
!  n = number of derivatives to compute
!

implicit none
use params

!input:
integer n,index
real*8 p(nxd,nyd,nzd)    ! original data
real*8 pt(nxd,nyd,nzd)   ! work array

!output:
real*8 pxx(nxd,nyd,nzd)
real*8 px(nxd,nyd,nzd)

integer n1,n1d,n2,n2d,n3,n3d

n1=nx
n1d=nxd
n2=ny
n2d=nyd
n3=nz
n3d=nzd


if (index==1) then
   px=p
   fft_derivatives(px,pxx,n,n1,n1d,n2,n2d,n3,n3d)
else if (index==2) then
   transpose12(p,pt,n1,n1d,n2,n2d,n3,n3d)
   ! 1st derivative returned in pt, 2nd derivative returned in px
   fft_derivatives(pt,px,n,n1,n1d,n2,n2d,n3,n3d)
   if (n==2) transpose12(px,pxx,n1,n1d,n2,n2d,n3,n3d)
   transpose12(pt,px,n1,n1d,n2,n2d,n3,n3d)
else if (index==3) then
   transpose13(p,pt,n1,n1d,n2,n2d,n3,n3d0
   fft_derivatives(pt,px,1,n1,n1d,n2,n2d,n3,n3d)
   if (n==2) transpose13(px,pxx,n1,n1d,n2,n2d,n3,n3d)
   transpose13(pt,px,n1,n1d,n2,n2d,n3,n3d)
endif

ASSERT(n1==nx)
ASSERT(n1d==nxd)
ASSERT(n2==ny)
ASSERT(n2d==nyd)
ASSERT(n3==nz)
ASSERT(n3d==nzd)

end









