#include "macros.h"

subroutine ns3d(rhs,Q,time,compute_ints,ints)
!
! evaluate RHS of N.S. equations:   -u dot grad(u) + mu * laplacian(u)
! if compute_ints==1, then we also return the following integrals:
!       ints(3)           ke disspaation from diffusion
!       ints(4)           vorticity z-component
!       ints(5)           helicity
!     
!
! vor(1) = w_y - v_z
! vor(2) = u_z - w_x 
! vor(3) = v_x - u_y
!
! hel = u (w_y-v_z) + v (u_z - w_x)  + w (v_x - u_y)
!
use params
implicit none

! input
real*8 Q(nx,ny,nz,n_var)
real*8 time
integer compute_ints

! output
real*8 rhs(nx,ny,nz,n_var)
real*8 ints(nints)

!local
real*8 d1(nx,ny,nz)
real*8 d2(nx,ny,nz)
real*8 work(nx,ny,nz)
real*8 dummy
real*8 :: ke_diss,vor,hel
integer i,j,k,numder


ke_diss=0
vor=0
hel=0
numder=DX_ONLY
if (mu>0) numder=DX_AND_DXX

rhs=0
! compute viscous terms (in rhs) and vorticity
do i=1,3

   ! compute u_x, u_xx
   call der(Q(1,1,1,i),d1,d2,work,numder,1)
!   d1=0
   rhs(:,:,:,i) = rhs(:,:,:,i) + mu*d2 - Q(:,:,:,1)*d1

   ke_diss=ke_diss + mu*sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,i)*d2(nx1:nx2,ny1:ny2,nz1:nz2))
   if (i==2) then  ! dv/dx, part of vor(3)
      vor=vor + sum(d1(nx1:nx2,ny1:ny2,nz1:nz2))
      hel=hel + sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,3)*d1(nx1:nx2,ny1:ny2,nz1:nz2))
   endif
   if (i==3) then  ! dw/dx, part of vor(2)
      hel=hel - sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,2)*d1(nx1:nx2,ny1:ny2,nz1:nz2))
   endif


   ! compute u_y, u_yy
   call der(Q(1,1,1,i),d1,d2,work,numder,2)
!   d1=0
   rhs(:,:,:,i) = rhs(:,:,:,i) + mu*d2 - Q(:,:,:,2)*d1


   ke_diss=ke_diss + mu*sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,i)*d2(nx1:nx2,ny1:ny2,nz1:nz2))
   if (i==1) then  ! du/dy part of vor(3)
      vor=vor - sum(d1(nx1:nx2,ny1:ny2,nz1:nz2))
      hel=hel - sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,3)*d1(nx1:nx2,ny1:ny2,nz1:nz2))
   endif
   if (i==3) then  ! dw/dy part of vor(1)
      hel=hel + sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,1)*d1(nx1:nx2,ny1:ny2,nz1:nz2))
   endif



   ! compute u_z, u_zz
   call der(Q(1,1,1,i),d1,d2,work,numder,3)
!   d1=0
   rhs(:,:,:,i) = rhs(:,:,:,i) + mu*d2 - Q(:,:,:,3)*d1

   ke_diss=ke_diss + mu*sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,i)*d2(nx1:nx2,ny1:ny2,nz1:nz2))
   if (i==1) then  ! du/dz part of vor(2)
      hel=hel + sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,2)*d1(nx1:nx2,ny1:ny2,nz1:nz2))
   endif
   if (i==2) then  ! dv/dz part of vor(1)
      hel=hel - sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,1)*d1(nx1:nx2,ny1:ny2,nz1:nz2))
   endif


enddo



! apply b.c. to rhs:
call bc_rhs(rhs)

if (compute_ints==1) then
   ints(3)=ke_diss
   ints(4)=vor
   ints(5)=hel
endif

end











