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
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: rhs(nx,ny,nz,n_var)
real*8 :: Q_tmp(nx,ny,nz,n_var)
real*8 :: Q_old(nx,ny,nz,n_var)
real*8 :: q4(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

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
call ns3D(rhs,Q,time,1,work1,work2,q4)
call divfree_gridspace(rhs,q4,work1,work2)
Q=Q+delt*rhs/6.0

! stage 2
Q_tmp = Q_old + delt*rhs/2.0
call ns3D(rhs,Q_tmp,time+delt/2.0,0,work1,work2,q4)
call divfree_gridspace(rhs,q4,work1,work2)
Q=Q+delt*rhs/3.0

! stage 3
Q_tmp = Q_old + delt*rhs/2.0
call ns3D(rhs,Q_tmp,time+delt/2.0,0,work1,work2,q4)
call divfree_gridspace(rhs,q4,work1,work2)
Q=Q+delt*rhs/3.0

! stage 4
Q_tmp = Q_old + delt*rhs
call ns3D(rhs,Q_tmp,time+delt,0,work1,work2,q4)
Q=Q+delt*rhs/6.0
call divfree_gridspace(Q,q4,work1,work2)





#else

! stage 1
call ns3D(rhs,Q,time,1,work1,work2,q4)
call divfree_gridspace(rhs,q4,work1,work2)
Q=Q+delt*rhs/3

! stage 2
Q_tmp = rhs
call ns3D(rhs,Q,time+delt/3,0,work1,work2,q4)
call divfree_gridspace(rhs,q4,work1,work2)
rhs = -5*Q_tmp/9 + rhs
Q=Q + 15*delt*rhs/16


! stage 3
Q_tmp=rhs
call ns3D(rhs,Q,time+3*delt/4,0,work1,work2,q4)
rhs = -153*Q_tmp/128 + rhs
Q=Q+8*delt*rhs/15
call divfree_gridspace(Q,q4,work1,work2)


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
real*8 Q(nx,ny,nz,n_var)
real*8 time
integer compute_ints

! output
real*8 rhs(nx,ny,nz,n_var)

!work
real*8 d1(nx,ny,nz)
real*8 d2(nx,ny,nz)
real*8 work(nx,ny,nz)

!local
real*8 dummy,tmx1,tmx2
real*8 :: ke,ke_diss,ke_diss2,vor,hel,gradu_diss
integer n,i,j,k,numder

call wallclock(tmx1)

ke=0
ke_diss=0
ke_diss2=0
gradu_diss=0
vor=0
hel=0
numder=DX_ONLY
if (mu>0) numder=DX_AND_DXX   ! only compute 2nd derivatives if mu>0


rhs=0
! compute viscous terms (in rhs) and vorticity
do n=1,ndim


   ! compute u_x, u_xx
   call der(Q(1,1,1,n),d1,d2,work,numder,1)

   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      ke = ke + .5*Q(i,j,k,n)**2

      rhs(i,j,k,n) = rhs(i,j,k,n) +  mu*d2(i,j,k) - Q(i,j,k,1)*d1(i,j,k) 

      ! Q,d1,d2 are in cache, so we can do these sums for free?
      ke_diss=ke_diss + Q(i,j,k,n)*d2(i,j,k)
      ke_diss2=ke_diss2 + d1(i,j,k)**2
      gradu_diss=gradu_diss + d2(i,j,k)**2
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

      ke_diss=ke_diss + Q(i,j,k,n)*d2(i,j,k)
      ke_diss2=ke_diss2  + d1(i,j,k)**2
      gradu_diss=gradu_diss + d2(i,j,k)**2
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


   if (ndim==3) then
   ! compute u_z, u_zz
   call der(Q(1,1,1,n),d1,d2,work,numder,3)

   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2

      rhs(i,j,k,n) = rhs(i,j,k,n) +  mu*d2(i,j,k) - Q(i,j,k,3)*d1(i,j,k) 

      ke_diss=ke_diss + Q(i,j,k,n)*d2(i,j,k)
      ke_diss2=ke_diss2 + d1(i,j,k)**2
      gradu_diss=gradu_diss + d2(i,j,k)**2
      if (n==1) then  ! du/dz part of vor(2)
         hel=hel + Q(i,j,k,2)*d1(i,j,k)
      endif
      if (n==2) then  ! dv/dz part of vor(1)
         hel=hel - Q(i,j,k,1)*d1(i,j,k)
      endif

   enddo
   enddo
   enddo
   endif


enddo



! apply b.c. to rhs:
call bc_rhs(rhs)
!call divfree_gridspace(rhs,work,d1,d2)

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
   ints(10)=ke_diss/g_nx/g_ny/g_nz     ! u dot laplacian u

endif

call wallclock(tmx2)
tims(5)=tims(5)+(tmx2-tmx1)
end











