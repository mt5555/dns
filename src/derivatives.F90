#include "macros.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  compute derivative along index index=1,2 or 3
!  n = 1  compute p_x, return in px.   (pxx is not accessed)
!  n = 2  compute p_xx, return in pxx.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine der(p,px,pxx,pt,n,index)
use params
use fft_interface
implicit none

!input:
integer n,index
real*8 p(nx,ny,nz)    ! original data
real*8 pt(nx,ny,nz)   ! work array

!output:
real*8 pxx(nx,ny,nz)
real*8 px(nx,ny,nz)

integer n1,n1d,n2,n2d,n3,n3d

n1=nx2
n1d=nx
n2=ny
n2d=ny2
n3=nz
n3d=nz2


if (index==1) then
   px=p
   call fft_derivatives(px,pxx,n,n1,n1d,n2,n2d,n3,n3d)
else if (index==2) then
   call transpose12(p,pt,1,n1,n1d,n2,n2d,n3,n3d)
   ! 1st derivative returned in pt, 2nd derivative returned in px
   call fft_derivatives(pt,px,n,n1,n1d,n2,n2d,n3,n3d)
   if (n==2) call transpose12(px,pxx,0,n1,n1d,n2,n2d,n3,n3d)
   call transpose12(pt,px,1,n1,n1d,n2,n2d,n3,n3d)
else if (index==3) then
   call transpose13(p,pt,1,n1,n1d,n2,n2d,n3,n3d)
   call fft_derivatives(pt,px,1,n1,n1d,n2,n2d,n3,n3d)
   if (n==2) call transpose13(px,pxx,0,n1,n1d,n2,n2d,n3,n3d)
   call transpose13(pt,px,1,n1,n1d,n2,n2d,n3,n3d)
endif

ASSERT("derivatives.F90 transpose error",n1==nx2)
ASSERT("derivatives.F90 transpose error",n1d==nx)
ASSERT("derivatives.F90 transpose error",n2==ny2)
ASSERT("derivatives.F90 transpose error",n2d==ny)
ASSERT("derivatives.F90 transpose error",n3==nz2)
ASSERT("derivatives.F90 transpose error",n3d==nz)

end subroutine









