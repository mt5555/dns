#include "macros.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine to take one Runge-Kutta 4th order time step
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk4(time,Q,Q2,Q3,Q4,rhs,work1,work2)
use params
implicit none
real*8 :: time
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Q2(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: Q3(nx,ny,nz,n_var)
real*8 :: Q4(nx,ny,nz,n_var)
real*8 :: rhs(nx,ny,nz,n_var)

! never use work1,work2,rhs,q4
call rk4reshape(time,Q,Q2(1,1,1,1),Q2(1,1,1,2),Q3(1,1,1,1),Q3(1,1,1,2))
end subroutine


subroutine rk4reshape(time,Q,rhs,w_old,w_tmp,psi)
use params
implicit none
real*8 :: time
real*8 :: Q(nx,ny,n_var)
real*8 :: w_tmp(nx,ny)
real*8 :: w_old(nx,ny)
real*8 :: psi(nx,ny)
real*8 :: rhs(nx,ny)


! local variables
real*8 :: vel,u,v
real*8 :: mone=-1,zero=0
integer i,j,k,n,ierr
logical,save :: firstcall=.true.


if (firstcall) then
   firstcall=.false.
   if (equations/=2) then
      call print_message("Error: psi-vor model can only run equations=2 (psi-vor)")
      call abort("initial conditions are probably incorrect.")
   endif
   if (ndim/=2) then
      call abort("Error: psi-vor model cannot run in 3D")
   endif
   if (nz/=1) then
      call abort("Error: psi-vor model cannot run in 3D")
   endif
   if (alpha_value/=0) then
      call abort("Error: dnsgrid cannot handle alpha>0.")
   endif
   Q(:,:,2)=Q(:,:,1)
   call helmholtz_inv(Q(1,1,2),rhs,zero,mone)

endif




!
!  Q(:,:,:,1) = vorticity
!  Q(:,:,:,2) = stream function
!  Q(:,:,:,3) = ?
!


w_old=Q(:,:,1)


! stage 1
call ns3D(rhs,Q,Q(1,1,2),time,1)
Q(:,:,1)=Q(:,:,1)+delt*rhs/6.0

! stage 2
w_tmp = w_old + delt*rhs/2.0
psi=w_tmp
call helmholtz_inv(psi,rhs,zero,mone)
call ns3D(rhs,w_tmp,psi,time+delt/2.0,0)
Q(:,:,1)=Q(:,:,1)+delt*rhs/3.0

! stage 3
w_tmp = w_old + delt*rhs/2.0
psi=w_tmp
call helmholtz_inv(psi,rhs,zero,mone)
call ns3D(rhs,w_tmp,psi,time+delt/2.0,0)
Q(:,:,1)=Q(:,:,1)+delt*rhs/3.0

! stage 4
w_tmp = w_old + delt*rhs
psi=w_tmp
call helmholtz_inv(psi,rhs,zero,mone)
call ns3D(rhs,w_tmp,psi,time+delt,0)
Q(:,:,1)=Q(:,:,1)+delt*rhs/6.0

Q(:,:,2)=Q(:,:,1)
call helmholtz_inv(Q(1,1,2),rhs,zero,mone)


time = time + delt


call ghost_update_x_reshape(Q(1,1,2),1)
call ghost_update_y_reshape(Q(1,1,2),1)

! compute KE, max U  
maxs(1:4)=0
do j=ny1,ny2
do i=nx1,nx2
   u=-( 2*(Q(i,j+1,2)-Q(i,j-1,2))/3 -  &
        (Q(i,j+2,2)-Q(i,j-2,2))/12          )/dely

   v=( 2*(Q(i+1,j,2)-Q(i-1,j,2))/3 -  &
        (Q(i+2,j,2)-Q(i-2,j,2))/12          )/delx

   maxs(1)=max(maxs(1),abs(u))
   maxs(2)=max(maxs(1),abs(v))
   maxs(3)=0
   vel = abs(u)/delx + abs(v)/dely 
   maxs(4)=max(maxs(4),vel)
enddo
enddo


end subroutine




subroutine ns3D(rhs,w,psi,time,compute_ints)
!
! evaluate RHS of N.S. equations:   -u dot grad(u) + mu * laplacian(u)
!
! if compute_ints==1, then we also return the following integrals:
!       ints(3)           ke disspaation from diffusion
!       ints(4)           vorticity z-component
!     
!
use params
use fft_interface
implicit none

! input
real*8 w(nx,ny)
real*8 psi(nx,ny)
real*8 time
integer compute_ints

! output
real*8 rhs(nx,ny)


!local
real*8 dummy,tmx1,tmx2
real*8 :: ke,ke_diss,ke_diss2,vor,dx,dxx,dy,dyy,ensave,u,v
integer n,i,j,k

call wallclock(tmx1)

ke=0
ke_diss=0
ke_diss2=0
vor=0
ensave=0



rhs=0

call ghost_update_x_reshape(psi,1)
call ghost_update_y_reshape(psi,1)
call ghost_update_x_reshape(w,1)
call ghost_update_y_reshape(w,1)

do j=ny1,ny2
do i=nx1,nx2
   u=-( 2*(psi(i,j+1)-psi(i,j-1))/3 -  &
        (psi(i,j+2)-psi(i,j-2))/12          )/dely

   v=( 2*(psi(i+1,j)-psi(i-1,j))/3 -  &
        (psi(i+2,j)-psi(i-2,j))/12          )/delx

   dx=( 2*(w(i+1,j)-w(i-1,j))/3 -  &
        (w(i+2,j)-w(i-2,j))/12          )/delx
   dxx=(-w(i+2,j) + 16*w(i+1,j) - 30*w(i,j) + &
        16*w(i-1,j) - w(i-2,j)) / (12*delx*delx)
   
   dy=( 2*(w(i,j+1)-w(i,j-1))/3 -  &
        (w(i,j+2)-w(i,j-2))/12          )/dely
   dyy=(-w(i,j+2) + 16*w(i,j+1) - 30*w(i,j) + &
        16*w(i,j-1) - w(i,j-2)) / (12*dely*dely)
   
   rhs(i,j) = rhs(i,j) +  mu*(dxx+dyy) - u*dx - v*dy
   
   ke = ke + .5*(u**2 + v**2)

   vor=vor + w(i,j)
   ensave = ensave + w(i,j)**2
enddo
enddo




! apply b.c. to rhs:
!call bc_rhs(rhs)


if (compute_ints==1) then
   ints(2)=ke_diss2/g_nx/g_ny     ! gradu dot gradu
   !ints(3) = forcing terms
   ints(4)=vor/g_nx/g_ny
   ints(6)=ke/g_nx/g_ny
   ints(7)=ensave /g_nx/g_ny
   ! ints(8) = < u,div(tau)' >   (alpha model only)
   ! ints(9)  = < u,f >  (alpha model only)
   ints(10)=ke_diss/g_nx/g_ny     ! u dot laplacian u

endif

call wallclock(tmx2)
tims(5)=tims(5)+(tmx2-tmx1)
end











subroutine ghost_update_x_reshape(psi,n)
use params
use ghost
implicit none
integer :: n
real*8 :: psi(nx,ny,nz,n)
call ghost_update_x(psi,n)
end
subroutine ghost_update_y_reshape(psi,n)
use params
use ghost
implicit none
integer :: n
real*8 :: psi(nx,ny,nz,n)
call ghost_update_y(psi,n)
end
