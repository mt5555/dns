#include "macros.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine to take one Runge-Kutta 4th order time step
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk4(time,Q,rhs,Q_old,Q_tmp,q4,work1,work2)
use params
implicit none
real*8 :: time
real*8 :: Q(nx,ny,n_var)
real*8 :: rhs(nx,ny,n_var)
real*8 :: Q_tmp(nx,ny,n_var)
real*8 :: Q_old(nx,ny,n_var)
real*8 :: q4(nx,ny,n_var)
real*8 :: work1(nx,ny)
real*8 :: work2(nx,ny)

! local variables
real*8 :: ints_buf(nints),vel
integer i,j,k,n,ierr
logical,save :: firstcall=.true.


if (firstcall) then
   firstcall=.false.
   if (ndim/=2) then
      call abort("Error: shallow water model cannot run in 3D")
   endif
   if (nz/=1) then
      call abort("Error: shallow water model cannot run in 3D")
   endif
   if (alpha_value/=0) then
      call abort("Error: shallow cannot yet handle alpha>0.")
   endif
endif



#define USE_RK4
#ifdef USE_RK4

Q_old=Q


! stage 1
call ns3D(rhs,Q,time,1,work1,work2,q4)
if (dealias) call dealias_gridspace(rhs,work1)
Q=Q+delt*rhs/6.0

! stage 2
Q_tmp = Q_old + delt*rhs/2.0
call ns3D(rhs,Q_tmp,time+delt/2.0,0,work1,work2,q4)
if (dealias) call dealias_gridspace(rhs,work1)
Q=Q+delt*rhs/3.0

! stage 3
Q_tmp = Q_old + delt*rhs/2.0
call ns3D(rhs,Q_tmp,time+delt/2.0,0,work1,work2,q4)
if (dealias) call dealias_gridspace(rhs,work1)
Q=Q+delt*rhs/3.0

! stage 4
Q_tmp = Q_old + delt*rhs
call ns3D(rhs,Q_tmp,time+delt,0,work1,work2,q4)
Q=Q+delt*rhs/6.0
if (dealias) call dealias_gridspace(Q,work1)




#else

! stage 1
call ns3D(rhs,Q,time,1,work1,work2,q4)
if (dealias) call dealias_gridspace(rhs,work1)
Q=Q+delt*rhs/3

! stage 2
Q_tmp = rhs
call ns3D(rhs,Q,time+delt/3,0,work1,work2,q4)
if (dealias) call dealias_gridspace(rhs,work1)
rhs = -5*Q_tmp/9 + rhs
Q=Q + 15*delt*rhs/16


! stage 3
Q_tmp=rhs
call ns3D(rhs,Q,time+3*delt/4,0,work1,work2,q4)
rhs = -153*Q_tmp/128 + rhs
Q=Q+8*delt*rhs/15
if (dealias) call dealias_gridspace(Q,work1)


#endif



time = time + delt


! compute KE, max U  
maxs(1:4)=0
do j=ny1,ny2
do i=nx1,nx2
   do n=1,3
      maxs(n)=max(maxs(n),abs(Q(i,j,n)))   ! max u,v,w
   enddo
   vel = abs(Q(i,j,1))/delx + abs(Q(i,j,2))/dely 
   maxs(4)=max(maxs(4),vel)
enddo
enddo


end subroutine rk4  





subroutine ns3d(rhs,Q,time,compute_ints,d1,d2,work)
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
implicit none

! input
real*8 Q(nx,ny,n_var)
real*8 time
integer compute_ints

! output
real*8 rhs(nx,ny,n_var)

!work
real*8 d1(nx,ny)
real*8 d2(nx,ny)
real*8 work(nx,ny,n_var)

!local
real*8 dummy,tmx1,tmx2
real*8 :: pe,ke,ke_diss,ke_diss2,vor,hel,gradu_diss
integer n,i,j,k,numder

call wallclock(tmx1)

ke=0
pe=0
ke_diss=0
ke_diss2=0
gradu_diss=0
vor=0
hel=0
numder=DX_ONLY
if (mu>0) numder=DX_AND_DXX   ! only compute 2nd derivatives if mu>0



! compute divergence ( h u)
do j=ny1,ny2
do i=nx1,nx2
   rhs(i,j,1)=Q(i,j,1)*Q(i,j,3)
   rhs(i,j,2)=Q(i,j,2)*Q(i,j,3)
enddo
enddo
call divergence(d1,rhs,d2,work)
do j=ny1,ny2
do i=nx1,nx2
   rhs(i,j,3)=-d1(i,j)
enddo
enddo



! coriolis force 
do j=ny1,ny2
do i=nx1,nx2
   rhs(i,j,1)= fcor*Q(i,j,2)
   rhs(i,j,2)=-fcor*Q(i,j,1)
enddo
enddo




! advection and viscous terms
do n=1,2


   ! compute u_x, u_xx
   call der(Q(1,1,n),d1,d2,work,numder,1)

   do j=ny1,ny2
   do i=nx1,nx2
      ke = ke + .5*Q(i,j,3)*Q(i,j,n)**2

      rhs(i,j,n) = rhs(i,j,n) +  mu*d2(i,j) - Q(i,j,1)*d1(i,j) 

      ! Q,d1,d2 are in cache, so we can do these sums for free?
      ke_diss=ke_diss - Q(i,j,3)*Q(i,j,n)*d2(i,j)
      !ke_diss2=ke_diss2 + d1(i,j)**2
      gradu_diss=gradu_diss + d2(i,j)**2
      if (n==2) then  ! dv/dx, part of vor(3)
         vor=vor + d1(i,j)
         hel=hel + Q(i,j,3)*d1(i,j)
      endif
      if (n==3) then  ! dw/dx, part of vor(2)
         hel=hel - Q(i,j,2)*d1(i,j)
      endif
   enddo
   enddo



   ! compute u_y, u_yy
   call der(Q(1,1,n),d1,d2,work,numder,2)

   do j=ny1,ny2
   do i=nx1,nx2

      rhs(i,j,n) = rhs(i,j,n) +  mu*d2(i,j) - Q(i,j,2)*d1(i,j) 

      ke_diss=ke_diss - Q(i,j,3)*Q(i,j,n)*d2(i,j)
      !ke_diss2=ke_diss2  + d1(i,j)**2
      gradu_diss=gradu_diss + d2(i,j)**2
      if (n==1) then  ! du/dy part of vor(3)
         vor=vor - d1(i,j)
         hel=hel - Q(i,j,3)*d1(i,j)
      endif
      if (n==3) then  ! dw/dy part of vor(1)
         hel=hel + Q(i,j,1)*d1(i,j)
      endif

   enddo
   enddo


enddo


! g grad (h+h_s)
do j=ny1,ny2
do i=nx1,nx2
   d2(i,j)=Q(i,j,3) ! + h_s(i,j)
   pe=pe+.5*grav*Q(i,j,3)**2
enddo
enddo
do n=1,2
   call der(d2,d1,dummy,work,DX_ONLY,n)
   do j=ny1,ny2
   do i=nx1,nx2
      rhs(i,j,n)=rhs(i,j,n)-grav*d1(i,j)
   enddo
   enddo
enddo



! apply b.c. to rhs:
call bc_rhs(rhs)


if (compute_ints==1) then
   ints(2)=ke_diss/g_nx/g_ny     ! u dot laplacian u
   !ints(2)=ke_diss2/g_nx/g_ny     ! gradu dot gradu

   !ints(3) = forcing terms

   ints(4)=vor/g_nx/g_ny
   ints(5)=hel/g_nx/g_ny
   ints(6)=(ke+pe)/g_nx/g_ny
   ints(10)=ke/g_nx/g_ny

   ints(7)=ints(2)  ! this is only true for periodic incompressible case
   ! ints(8) = < u,div(tau)' >   (alpha model only)
   ! ints(9)  = < u,f >  (alpha model only)
   ints(1)=gradu_diss/g_nx/g_ny

endif

call wallclock(tmx2)
tims(5)=tims(5)+(tmx2-tmx1)
end











