#include "macros.h"

subroutine ns3d(rhs,Q,time)
!
! vor(1) = w_y - v_z
! vor(2) = u_z - w_x 
! vor(3) = v_x - u_y
!
use params
implicit none


! input
real*8 Q(nx,ny,nz,n_var)
real*8 time
! output
real*8 rhs(nx,ny,nz,n_var)
!local
real*8 d1(nx,ny,nz)
real*8 d2(nx,ny,nz)
real*8 work(nx,ny,nz)
real*8 dummy
real*8,external :: norm_divergence 
integer i,j,k




rhs=0
! compute viscous terms (in rhs) and vorticity
do i=1,3
   ! compute u_x, u_xx
   call der(Q(1,1,1,i),d1,d2,work,DX_AND_DXX,1)
   rhs(:,:,:,i) = rhs(:,:,:,i) + mu*d2 - Q(:,:,:,1)*d1
!   if (i==3) vor(:,:,:,2) = vor(:,:,:,2) - d1
!   if (i==2) vor(:,:,:,3) = vor(:,:,:,3) + d1

   ! compute u_y, u_yy
   call der(Q(1,1,1,i),d1,d2,work,DX_AND_DXX,2)
   rhs(:,:,:,i) = rhs(:,:,:,i) + mu*d2 - Q(:,:,:,2)*d1
!   if (i==3) vor(:,:,:,1) = vor(:,:,:,1) + d1
!   if (i==1) vor(:,:,:,3) = vor(:,:,:,3) -d1



   ! compute u_z, u_zz
   call der(Q(1,1,1,i),d1,d2,work,DX_AND_DXX,3)
   rhs(:,:,:,i) = rhs(:,:,:,i) + mu*d2 - Q(:,:,:,3)*d1
!   if (i==2) vor(:,:,:,1) = vor(:,:,:,1) -d1
!   if (i==1) vor(:,:,:,2) = vor(:,:,:,2) +d1

enddo


if (dealias) then
do i=1,3
   call fft3d(rhs(1,1,1,i),d1)
   call fft_filter_dealias(rhs(1,1,1,i))
   call ifft3d(rhs(1,1,1,i),d1)
enddo
endif

! apply b.c. to rhs:
call bc_rhs(rhs)



! make rhs divergence free:
call divfree(rhs,d2,d1,work)


end











