#include "macros.h"

subroutine ns3d(rhs,Q,time,compute_ints)
!
! evaluate RHS of N.S. equations:   -u dot grad(u) + mu * laplacian(u)
!
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

!local
real*8 d1(nx,ny,nz)
real*8 d2(nx,ny,nz)
real*8 work(nx,ny,nz)
real*8 dummy,tmx1,tmx2
real*8 :: ke_diss,ke_diss2,vor,hel
integer n,i,j,k,numder

call wallclock(tmx1)

ke_diss=0
ke_diss2=0
vor=0
hel=0
numder=DX_ONLY
if (mu>0) numder=DX_AND_DXX   ! only compute 2nd derivatives if mu>0


rhs=0
! compute viscous terms (in rhs) and vorticity
do n=1,3


   ! compute u_x, u_xx
   call der(Q(1,1,1,n),d1,d2,work,numder,1)

   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2

      rhs(i,j,k,n) = rhs(i,j,k,n) +  mu*d2(i,j,k) - Q(i,j,k,1)*d1(i,j,k) 

      ! Q,d1,d2 are in cache, so we can do these sums for free?
      ke_diss=ke_diss + mu*Q(i,j,k,n)*d2(i,j,k)
      ke_diss2=ke_diss2 - mu*d1(i,j,k)**2
      if (n==2) then  ! dv/dx, part of vor(3)
         vor=vor + d1(i,j,k)
         hel=hel + Q(i,j,k,3)*d1(i,j,k)
      endif
      if (n==3) then  ! dw/dx, part of vor(2)
         hel=hel - Q(i,j,k,2)*d1(i,j,k)
      endif
   enddo
   enddo
   enddo



   ! compute u_y, u_yy
   call der(Q(1,1,1,n),d1,d2,work,numder,2)

   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2

      rhs(i,j,k,n) = rhs(i,j,k,n) +  mu*d2(i,j,k) - Q(i,j,k,2)*d1(i,j,k) 

      ke_diss=ke_diss + mu*Q(i,j,k,n)*d2(i,j,k)
      ke_diss2=ke_diss2 - mu*d1(i,j,k)**2
      if (n==1) then  ! du/dy part of vor(3)
         vor=vor - d1(i,j,k)
         hel=hel - Q(i,j,k,3)*d1(i,j,k)
      endif
      if (n==3) then  ! dw/dy part of vor(1)
         hel=hel + Q(i,j,k,1)*d1(i,j,k)
      endif

   enddo
   enddo
   enddo



   ! compute u_z, u_zz
   call der(Q(1,1,1,n),d1,d2,work,numder,3)

   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2

      rhs(i,j,k,n) = rhs(i,j,k,n) +  mu*d2(i,j,k) - Q(i,j,k,3)*d1(i,j,k) 

      ke_diss=ke_diss + mu*Q(i,j,k,n)*d2(i,j,k)
      ke_diss2=ke_diss2 - mu*d1(i,j,k)**2
      if (n==1) then  ! du/dz part of vor(2)
         hel=hel + Q(i,j,k,2)*d1(i,j,k)
      endif
      if (n==2) then  ! dv/dz part of vor(1)
         hel=hel - Q(i,j,k,1)*d1(i,j,k)
      endif

   enddo
   enddo
   enddo



enddo



! apply b.c. to rhs:
call bc_rhs(rhs)

if (compute_ints==1) then
   ints_timeDU=time
   ints(3)=ke_diss2
   ints(4)=vor
   ints(5)=hel
endif

call wallclock(tmx2)
tims(5)=tims(5)+(tmx2-tmx1)
end











