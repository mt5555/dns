#include "macros.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine to take one Runge-Kutta 4th order time step
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk4(time,Q,rhs,Q_old,Q_tmp,q4,work1,p)
use params
implicit none
real*8 :: time
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: rhs(nx,ny,nz,n_var)
real*8 :: Q_tmp(nx,ny,nz,n_var)
real*8 :: Q_old(nx,ny,nz,n_var)
real*8 :: q4(nx,ny,nz,n_var)     ! not used
real*8 :: work1(nx,ny,nz)
real*8 :: p(nx,ny,nz)

! local variables
real*8 :: vel
integer i,j,k,n,ierr
logical,save :: firstcall=.true.


if (firstcall) then
   firstcall=.false.
   if (alpha_value/=0) then
      call abort("Error: dnsgrid cannot handle alpha>0.")
   endif
endif






#define USE_RK4
#ifdef USE_RK4

Q_old=Q


! stage 1
call ns3D(rhs,Q,time,1)
call divfree_ghost(rhs,p,work1)
Q=Q+delt*rhs/6.0

! stage 2
Q_tmp = Q_old + delt*rhs/2.0
call ns3D(rhs,Q_tmp,time+delt/2.0,0)
call divfree_ghost(rhs,p,work1)
Q=Q+delt*rhs/3.0

! stage 3
Q_tmp = Q_old + delt*rhs/2.0
call ns3D(rhs,Q_tmp,time+delt/2.0,0)
call divfree_ghost(rhs,p,work1)
Q=Q+delt*rhs/3.0

! stage 4
Q_tmp = Q_old + delt*rhs
call ns3D(rhs,Q_tmp,time+delt,0)
Q=Q+delt*rhs/6.0
call divfree_ghost(Q,p,work1)





#else

! stage 1
call ns3D(rhs,Q,time,1)
call divfree_ghost(rhs,p,work1)
Q=Q+delt*rhs/3

! stage 2
Q_tmp = rhs
call ns3D(rhs,Q,time+delt/3,0)
call divfree_ghost(rhs,p,work1)
rhs = -5*Q_tmp/9 + rhs
Q=Q + 15*delt*rhs/16


! stage 3
Q_tmp=rhs
call ns3D(rhs,Q,time+3*delt/4,0)
rhs = -153*Q_tmp/128 + rhs
Q=Q+8*delt*rhs/15
call divfree_ghost(Q,p,work1)


#endif




time = time + delt


! compute KE, max U  
maxs(1:4)=0
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   do n=1,3
      maxs(n)=max(maxs(n),abs(Q(i,j,k,n)))   ! max u,v,w
   enddo
   if (ndim==3) then
      vel = abs(Q(i,j,k,1))/delx + abs(Q(i,j,k,2))/dely + abs(Q(i,j,k,3))/delz
   else
      vel = abs(Q(i,j,k,1))/delx + abs(Q(i,j,k,2))/dely 
   endif
   maxs(4)=max(maxs(4),vel)
enddo
enddo
enddo


end subroutine rk4  





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
use fft_interface
use ghost
implicit none

! input
real*8 Q(nx,ny,nz,n_var)
real*8 time
integer compute_ints

! output
real*8 rhs(nx,ny,nz,n_var)


!local
real*8 dummy,tmx1,tmx2
real*8 :: ke,ke_diss,ke_diss2,vor,hel,gradu_diss,d1,d2
integer n,i,j,k

call wallclock(tmx1)

ke=0
ke_diss=0
ke_diss2=0
gradu_diss=0
vor=0
hel=0








rhs=0

call ghost_update_x(Q,ndim)
do n=1,ndim

   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      ke = ke + .5*Q(i,j,k,n)**2

      d1=( 2*(Q(i+1,j,k,n)-Q(i-1,j,k,n))/3 -  &
                 (Q(i+2,j,k,n)-Q(i-2,j,k,n))/12          )/delx
      d2=(-Q(i+2,j,k,n) + 16*Q(i+1,j,k,n) - 30*Q(i,j,k,n) + &
                   16*Q(i-1,j,k,n) - Q(i-2,j,k,n)) / (12*delx*delx)

      rhs(i,j,k,n) = rhs(i,j,k,n) +  mu*d2 - Q(i,j,k,1)*d1

      ! Q,d1,d2 are in cache, so we can do these sums for free?
      ke_diss=ke_diss + Q(i,j,k,n)*d2
      ke_diss2=ke_diss2 + d1**2
      gradu_diss=gradu_diss + d2**2
      if (n==2) then  ! dv/dx, part of vor(3)
         vor=vor + d1
         hel=hel + Q(i,j,k,3)*d1
      endif
      if (n==3) then  ! dw/dx, part of vor(2)
         hel=hel - Q(i,j,k,2)*d1
      endif
   enddo
   enddo
   enddo
enddo

call ghost_update_y(Q,ndim)
do n=1,ndim

   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2

      d1=( 2*(Q(i,j+1,k,n)-Q(i,j-1,k,n))/3 -  &
                  (Q(i,j+2,k,n)-Q(i,j-2,k,n))/12          )/dely
      d2=(-Q(i,j+2,k,n) + 16*Q(i,j+1,k,n) - 30*Q(i,j,k,n) + &
                   16*Q(i,j-1,k,n) - Q(i,j-2,k,n)) / (12*dely*dely)
      rhs(i,j,k,n) = rhs(i,j,k,n) +  mu*d2 - Q(i,j,k,2)*d1

      ke_diss=ke_diss + Q(i,j,k,n)*d2
      ke_diss2=ke_diss2  + d1**2
      gradu_diss=gradu_diss + d2**2
      if (n==1) then  ! du/dy part of vor(3)
         vor=vor - d1
         hel=hel - Q(i,j,k,3)*d1
      endif
      if (n==3) then  ! dw/dy part of vor(1)
         hel=hel + Q(i,j,k,1)*d1
      endif

   enddo
   enddo
   enddo
enddo


if (ndim==3) then
call ghost_update_z(Q,ndim)
do n=1,ndim


   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2

      d1=( 2*(Q(i,j,k+1,n)-Q(i,j,k-1,n))/3 -  &
                  (Q(i,j,k+2,n)-Q(i,j,k-2,n))/12          )/delz
      d2=(-Q(i,j,k+2,n) + 16*Q(i,j,k+1,n) - 30*Q(i,j,k,n) + &
                   16*Q(i,j,k-1,n) - Q(i,j,k-2,n)) / (12*delz*delz)
      rhs(i,j,k,n) = rhs(i,j,k,n) +  mu*d2 - Q(i,j,k,3)*d1

      ke_diss=ke_diss + Q(i,j,k,n)*d2
      ke_diss2=ke_diss2 + d1**2
      gradu_diss=gradu_diss + d2**2
      if (n==1) then  ! du/dz part of vor(2)
         hel=hel + Q(i,j,k,2)*d1
      endif
      if (n==2) then  ! dv/dz part of vor(1)
         hel=hel - Q(i,j,k,1)*d1
      endif

   enddo
   enddo
   enddo
enddo
endif



if (compute_ints==1) then
   ints(1)=gradu_diss/g_nx/g_ny/g_nz

   ints(2)=ke_diss2/g_nx/g_ny/g_nz     ! gradu dot gradu

   !ints(3) = forcing terms
   ints(4)=vor/g_nx/g_ny/g_nz
   ints(5)=hel/g_nx/g_ny/g_nz
   ints(6)=ke/g_nx/g_ny/g_nz
   ints(7)=ints(2)  ! this is only true for periodic incompressible case
   ! ints(8) = < u,div(tau)' >   (alpha model only)
   ! ints(9)  = < u,f >  (alpha model only)
   ints(10)=mu*ke_diss/g_nx/g_ny/g_nz     ! u dot laplacian u

endif

call wallclock(tmx2)
tims(5)=tims(5)+(tmx2-tmx1)
end











subroutine divfree_ghost(u,p,work)
!
! make u divergence free
!    solve:  div(u) = laplacian(p)
!    then:   unew = u - grad(p)
!    
! 
!
use params
use fft_interface
use ghost
implicit none
real*8 :: u(nx,ny,nz,3)
real*8 :: p(nx,ny,nz)
real*8 :: work(nx,ny,nz)

!local
real*8 :: dummy(1),tol
real*8 :: alpha=0
real*8 :: beta=1
integer i,j,k,n
external helmholtz_periodic,helmholtz_dirichlet


! solve laplacian(p)=div(u)

if (bdy_x1==PERIODIC .and. bdy_y1==PERIODIC .and. bdy_z1==PERIODIC) then

   ! p = ux+vy+wz
   call ghost_update_x(u(1,1,1,1),1)
   call ghost_update_y(u(1,1,1,2),1)
   if (ndim==3) call ghost_update_z(u(1,1,1,3),1)
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      p(i,j,k)=  &
        ( 2*(u(i+1,j,k,1)-u(i-1,j,k,1))/3 -  (u(i+2,j,k,1)-u(i-2,j,k,1))/12 )/delx
   
      p(i,j,k)=p(i,j,k)+&
        ( 2*(u(i,j+1,k,2)-u(i,j-1,k,2))/3 -  (u(i,j+2,k,2)-u(i,j-2,k,2))/12 )/dely
   
   if (ndim==3) then
      p(i,j,k)=p(i,j,k)+&
           ( 2*(u(i,j,k+1,3)-u(i,j,k-1,3))/3 -  (u(i,j,k+2,3)-u(i,j,k-2,3))/12 )/delz
   endif

   enddo
   enddo
   enddo


   call helmholtz_periodic_inv(p,work,alpha,beta)
   !work=p  ! RHS
   !p=0  ! initial guess
   !tol=1e-10
   !call cgsolver(p,work,alpha,beta,tol,work2,helmholtz_periodic,.false.)

else
   stop 'divfree_ghost: only supports periodic case'
endif




! compute u=u-grad(p)
call ghost_update_x(p,1)
call ghost_update_y(p,1)
if (ndim==3) call ghost_update_z(p,1)

do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   u(i,j,k,1)=u(i,j,k,1)-   &
        ( 2*(p(i+1,j,k)-p(i-1,j,k))/3 -  (p(i+2,j,k)-p(i-2,j,k))/12 )/delx
   
   u(i,j,k,2)=u(i,j,k,2)-   &
        ( 2*(p(i,j+1,k)-p(i,j-1,k))/3 -  (p(i,j+2,k)-p(i,j-2,k))/12 )/dely
   
   if (ndim==3) then
      u(i,j,k,3)=u(i,j,k,3)-   &
           ( 2*(p(i,j,k+1)-p(i,j,k-1))/3 -  (p(i,j,k+2)-p(i,j,k-2))/12 )/delz
   endif
enddo
enddo
enddo


end subroutine
