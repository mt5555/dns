#include "macros.h"

subroutine ns3d(rhs,Q,time,ke_diss)
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
real*8 :: ke_diss


! output
real*8 rhs(nx,ny,nz,n_var)

!local
real*8 d1(nx,ny,nz)
real*8 d2(nx,ny,nz)
real*8 work(nx,ny,nz)
real*8 dummy
integer i,j,k


ke_diss=0
rhs=0
! compute viscous terms (in rhs) and vorticity
do i=1,3

   ! compute u_x, u_xx
   call der(Q(1,1,1,i),d1,d2,work,DX_AND_DXX,1)
   rhs(:,:,:,i) = rhs(:,:,:,i) + mu*d2 - Q(:,:,:,1)*d1
!   if (i==3) vor(:,:,:,2) = vor(:,:,:,2) - d1
!   if (i==2) vor(:,:,:,3) = vor(:,:,:,3) + d1
   ke_diss=ke_diss + mu*sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,i)*d2(nx1:nx2,ny1:ny2,nz1:nz2))

   ! compute u_y, u_yy
   call der(Q(1,1,1,i),d1,d2,work,DX_AND_DXX,2)
   rhs(:,:,:,i) = rhs(:,:,:,i) + mu*d2 - Q(:,:,:,2)*d1
!   if (i==3) vor(:,:,:,1) = vor(:,:,:,1) + d1
!   if (i==1) vor(:,:,:,3) = vor(:,:,:,3) -d1
   ke_diss=ke_diss + mu*sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,i)*d2(nx1:nx2,ny1:ny2,nz1:nz2))



   ! compute u_z, u_zz
   call der(Q(1,1,1,i),d1,d2,work,DX_AND_DXX,3)
   rhs(:,:,:,i) = rhs(:,:,:,i) + mu*d2 - Q(:,:,:,3)*d1
!   if (i==2) vor(:,:,:,1) = vor(:,:,:,1) -d1
!   if (i==1) vor(:,:,:,2) = vor(:,:,:,2) +d1
   ke_diss=ke_diss + mu*sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,i)*d2(nx1:nx2,ny1:ny2,nz1:nz2))


enddo



! apply b.c. to rhs:
call bc_rhs(rhs)


end











