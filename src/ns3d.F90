#include "macros.h"

subroutine ns3d(Q,rhs,time)
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
real*8 vor(nx,ny,nz,3)     ! could be removed and data accumulated in rhs
real*8 dummy
real*8 :: alpha=0
real*8 :: beta=1
integer i,j,k


rhs=0
vor=0
! compute viscous terms (in rhs) and vorticity
do i=1,3
   ! compute u_x, u_xx
   call der(Q(1,1,1,i),d1,d2,work,DX_AND_DXX,1)
   rhs(:,:,:,i) = rhs(:,:,:,i) + d2
   if (i==3) vor(:,:,:,2) = vor(:,:,:,2) - d1
   if (i==2) vor(:,:,:,3) = vor(:,:,:,3) + d1

   ! compute u_y, u_yy
   call der(Q(1,1,1,i),d1,d2,work,DX_AND_DXX,2)
   rhs(:,:,:,i) = rhs(:,:,:,i) + d2
   if (i==3) vor(:,:,:,1) = vor(:,:,:,1) + d1
   if (i==1) vor(:,:,:,3) = vor(:,:,:,3) -d1

   ! compute u_z, u_zz
   call der(Q(1,1,1,i),d1,d2,work,DX_AND_DXX,3)
   rhs(:,:,:,i) = rhs(:,:,:,i) + d2
   if (i==2) vor(:,:,:,1) = vor(:,:,:,1) -d1
   if (i==1) vor(:,:,:,2) = vor(:,:,:,2) +d1
enddo



! build up the rhs
! rhs = mu*rhs + Q cross vor 
! vor = Q cross vor

rhs(:,:,:,1)=mu*rhs(:,:,:,1) + Q(:,:,:,2)*vor(:,:,:,3) - Q(:,:,:,3)*vor(:,:,:,2)
rhs(:,:,:,2)=mu*rhs(:,:,:,2) + Q(:,:,:,3)*vor(:,:,:,1) - Q(:,:,:,1)*vor(:,:,:,3)
rhs(:,:,:,3)=mu*rhs(:,:,:,2) + Q(:,:,:,1)*vor(:,:,:,2) - Q(:,:,:,2)*vor(:,:,:,1)


! compute d2 = div(q cross vor), or use the full rhs
i=1
call der(rhs(1,1,1,i),d1,dummy,work,DX_ONLY,i)
d2 = d1
i=2
call der(rhs(1,1,1,i),d1,dummy,work,DX_ONLY,i)
d2 = d2+d1
i=3
call der(rhs(1,1,1,i),d1,dummy,work,DX_ONLY,i)
d2 = d2+d1



! solve laplacian p = div(q cross vor).  d2 is overritten with p. 
call poisson(d2,work,alpha,beta)

! add grad p to the RHS  (p still stored in d2)
do i=1,3
   call der(d2,d1,dummy,work,DX_ONLY,i)
   rhs(:,:,:,i) = rhs(:,:,:,i) + d1
enddo

end











