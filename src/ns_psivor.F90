#include "macros.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine to take one Runge-Kutta 4th order time step
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if 0

outline of algorithm for infinite domain computation.

w = vorticity
PSI = stream function

1. w initialized on interior and boundary
2. call gost_update on w
3. compute_psi 
4. compute (u,v) from grad PSI
5. compute grad(w) and w_xx, boundary terms hand coded in loop:
      go to 2nd order centered near boundary


6. advance w in time:   w_t + u dot grad(w) = viscosity
7. set w on boundary (needed for compute_PSI) and call ghost_update(w)
      inflow: 0
      outflow: chose so 2nd order will give one sided difference?
goto #3


compute_PSI:
  1. solve for PSI on boundary with B-S
  2  solve for PSI on interior
  3. set ghost cells for PSI so that 4th order formula because
     2nd order near boundary
  4. call ghost update on PSI

#endif



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

! not used: work1, rhs, Q3,Q4
call rk4reshape(time,Q,Q(1,1,1,3),Q2(1,1,1,1),Q2(1,1,1,2),Q2(1,1,1,3),work1)
end subroutine


subroutine rk4reshape(time,Q,rhs,w_old,w_tmp,psi,work)
use params
implicit none
real*8 :: time
real*8 :: Q(nx,ny,n_var)
real*8 :: w_tmp(nx,ny)
real*8 :: w_old(nx,ny)
real*8 :: psi(nx,ny)
real*8 :: rhs(nx,ny)
real*8 :: work(nx,ny)


! local variables
real*8 :: vel,u,v
integer i,j,k,n,ierr
logical,save :: firstcall=.true.
integer,save :: btype


if (firstcall) then
   firstcall=.false.
   if (equations/=NS_PSIVOR) then
      call print_message("Error: psi-vor model can only run equations=NS_PSIVOR")
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

   if (bdy_x1==PERIODIC .and. bdy_y1==PERIODIC) then
      ! periodic
      btype=0
   else if  ( &
        ! not periodic, but can still do with all ghostcells:
        ((bdy_x1==PERIODIC) .or. (bdy_x1==REFLECT) .or. (bdy_x1==REFLECT_ODD)) .and. &
        ((bdy_x2==PERIODIC) .or. (bdy_x2==REFLECT) .or. (bdy_x2==REFLECT_ODD)) .and. &
        ((bdy_y1==PERIODIC) .or. (bdy_y1==REFLECT) .or. (bdy_y1==REFLECT_ODD)) .and. &
        ((bdy_y2==PERIODIC) .or. (bdy_y2==REFLECT) .or. (bdy_y2==REFLECT_ODD))  ) then
        btype=1
   else
      btype=2
   endif

   ! initial vorticity should have been set on the boundary, so
   ! we can compute PSI right now:
   call compute_psi(Q(1,1,1),Q(1,1,2),rhs,work,btype)
   call bc_impose(Q(1,1,1),Q(1,1,2))

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
call bc_impose(w_tmp,Q(1,1,2))
call compute_psi(w_tmp,psi,rhs,work,btype)
call ns3D(rhs,w_tmp,psi,time+delt/2.0,0)
Q(:,:,1)=Q(:,:,1)+delt*rhs/3.0



! stage 3
w_tmp = w_old + delt*rhs/2.0
call bc_impose(w_tmp,Q(1,1,2))
call compute_psi(w_tmp,psi,rhs,work,btype)
call ns3D(rhs,w_tmp,psi,time+delt/2.0,0)
Q(:,:,1)=Q(:,:,1)+delt*rhs/3.0

! stage 4
w_tmp = w_old + delt*rhs
call bc_impose(w_tmp,Q(1,1,2))
call compute_psi(w_tmp,psi,rhs,work,btype)
call ns3D(rhs,w_tmp,psi,time+delt,0)
Q(:,:,1)=Q(:,:,1)+delt*rhs/6.0
call bc_impose(Q(1,1,1),Q(1,1,2))

call compute_psi(Q(1,1,1),Q(1,1,2),rhs,work,btype)
time = time + delt





! compute KE, max U  
maxs(1:4)=0
do j=ny1,ny2
do i=nx1,nx2
   u=( 2*(Q(i,j+1,2)-Q(i,j-1,2))/3 -  &
        (Q(i,j+2,2)-Q(i,j-2,2))/12          )/dely

   v=-( 2*(Q(i+1,j,2)-Q(i-1,j,2))/3 -  &
        (Q(i+2,j,2)-Q(i-2,j,2))/12          )/delx

   maxs(1)=max(maxs(1),abs(u))
   maxs(2)=max(maxs(1),abs(v))
   maxs(3)=0
   vel = abs(u)/delx + abs(v)/dely 
   maxs(4)=max(maxs(4),vel)
enddo
enddo


end subroutine




subroutine compute_psi(w,psi,b,work,btype)
!
!  btype=0   all periodic  use periodic FFT solver 
!  btype=1   all periodic,reflect, reflect-odd
!                use ghost cell CG solver (no boundaries)
!  btype=2   use biotsavart for boundary and dirichlet solver
!
use params
implicit none
real*8 w(nx,ny)
real*8 psi(nx,ny)
real*8 work(nx,ny)
real*8 b(nx,ny)
integer :: btype

!local
real*8 :: one=1,zero=0,tol=1e-8
external helmholtz_dirichlet,helmholtz_periodic,helmholtz_periodic_ghost



! if all b.c. periodic:
if (btype==0) then
   psi=-w
   call helmholtz_periodic_inv(psi,work,zero,one)

   !psi=0  ! initial guess
   !b=-w
   !call cgsolver(psi,b,zero,one,tol,work,helmholtz_periodic,.true.)

else if (btype==1) then
   ! this code can handle any compination of PERIODIC, REFLECT, REFLECT-ODD:

   b=-w  ! be sure to copy ghost cell data also!
   ! apply compact correction to b:
   call helmholtz_dirichlet_setup(b,psi,work,0)
   call cgsolver(psi,b,zero,one,tol,work,helmholtz_periodic_ghost,.false.)

else
   psi=0  ! initial guess
   call bc_biotsavart(w,psi)    !update PSI on boundary using biot-savart law

   b=-w  ! be sure to copy ghost cell data also!

   ! copy b.c. from psi into 'b', and apply compact correction to b:
   call helmholtz_dirichlet_setup(b,psi,work,1)
   call cgsolver(psi,b,zero,one,tol,work,helmholtz_dirichlet,.false.)


   !update PSI 1st row of ghost cells so that our 4th order differences
   !near the boundary look like 2nd order centered
   call bc_onesided(psi)
endif

call ghost_update_x_reshape(psi,1)
call ghost_update_y_reshape(psi,1)

end subroutine







subroutine bc_impose(w,psi)
!
! on non-periodic or non-reflective boundarys:
!
! set w=0 on boundary for inflow
! interpolate for outflow
!  
use params
implicit none
real*8 w(nx,ny)
real*8 psi(nx,ny)


!local
integer i,j
real*8 :: u,v

if (my_x==0 .and. bdy_x1==INFLOW0_ONESIDED) then
   !             nx1   nx1+1   nx1+2    nx1+3
   ! stencil (   -1             1              ) /2h
   !                    -3      4        1     ) /2h
   !      
   !
   !         -(nx1) + (nx1+2)  = -3(nx1+1) + 4(nx1+2) + (nx1+3)
   !          nx1 = 3(nx1+1) - 3(nx1+2) - (nx1+3)
   do j=ny1,ny2

      i=nx1+1
      u=( 2*(psi(i,j+1)-psi(i,j-1))/3 -  &
           (psi(i,j+2)-psi(i,j-2))/12          )

      if (u>=0) then !inflow
         w(nx1,j)=0
      else
         w(nx1,j)= 3*w(nx1+1,j)-3*w(nx1+2,j)-w(nx1+3,j)
      endif
   enddo
endif
if (my_x==ncpu_x-1 .and. bdy_x2==INFLOW0_ONESIDED) then
   do j=ny1,ny2

      i=nx2-1
      u=( 2*(psi(i,j+1)-psi(i,j-1))/3 -  &
           (psi(i,j+2)-psi(i,j-2))/12          )


      if (u<=0) then
         w(nx2,j)=0
      else
         w(nx2,j)= 3*w(nx2-1,j)-3*w(nx2-2,j)-w(nx2-3,j)
      endif
   enddo
endif


if (my_y==0 .and. bdy_y1==INFLOW0_ONESIDED) then
   do i=nx1,nx2

      j=ny1+1
      v=-( 2*(psi(i+1,j)-psi(i-1,j))/3 -  &
           (psi(i+2,j)-psi(i-2,j))/12          )

      if (v>=0) then
         w(i,ny1)=0
      else
         w(i,ny1)= 3*w(i,ny1+1)  - 3*w(i,ny1+2) - w(i,ny1+3)
      endif
   enddo
endif

if (my_y==ncpu_y-1 .and. bdy_y2==INFLOW0_ONESIDED) then
   do i=nx1,nx2

      j=ny2-1
      v=-( 2*(psi(i+1,j)-psi(i-1,j))/3 -  &
           (psi(i+2,j)-psi(i-2,j))/12          )

      if (v<=0) then
         w(i,ny2)=0
      else
         w(i,ny2)=  3*w(i,ny2-1) - 3*w(i,ny2-2) - w(i,ny2-3)
      endif
   enddo
endif




call ghost_update_x_reshape(w,1)
call ghost_update_y_reshape(w,1)

end subroutine






subroutine bc_onesided(w)
! on non-periodic or non-reflective boundarys:
!
! fiddle first ghost cell so that 4th order derivative at 1st interior point
! will be the same as a 2nd order centered scheme.  
! (assuming boundary values already set)

use params
implicit none
real*8 w(nx,ny)


!local
integer i,j

if (my_x==0 .and. bdy_x1==INFLOW0_ONESIDED) then
   !           nx1-1     nx1   nx1+1   nx1+2    nx1+3
   ! stencil ( 1/12     -2/3      0     2/3     -1/12)     /h
   !                    -1/2            1/2                /h
   ! multiply both sided by 12h:
   !           nx1-1     nx1   nx1+1   nx1+2    nx1+3
   ! stencil    1       -8      0       8       -1
   !                    -6              6         

   do j=ny1,ny2
      w(nx1-1,j)= 2*w(nx1,j)  - 2*w(nx1+2,j) +  w(nx1+3,j)
   enddo
endif
if (my_x==ncpu_x-1 .and. bdy_x2==INFLOW0_ONESIDED) then
   !           nx2-3   nx2-2   nx2-1   nx2     nx2+1
   ! stencil    1       -8      0       8       -1
   !                    -6              6         
   do j=ny1,ny2
      w(nx2+1,j)= 2*w(nx2,j)  -2*w(nx2-2,j)  +  w(nx2-3,j)
   enddo
endif


if (my_y==0 .and. bdy_y1==INFLOW0_ONESIDED) then
   do i=nx1,nx2
      w(i,ny1-1)= 2*w(i,ny1)  - 2*w(i,ny1+2) +  w(i,ny1+3)
   enddo
endif

if (my_y==ncpu_y-1 .and. bdy_y2==INFLOW0_ONESIDED) then
   do i=nx1,nx2
      w(i,ny2+1)=  2*w(i,ny2)  -2*w(i,ny2-2)  +  w(i,ny2-3)
   enddo
endif



end subroutine






subroutine bc_biotsavart(w,psi)
! on non-periodic or non-reflective boundarys:
! use biot-savar law to compute boundary data for PSI.
!
! except for y1 boundary, where we set PSI=0
!
!     Finds psi on boundary from w using Biot-Savart
!
!     psi = -(1/4pi)*sum w(i,j)*logterm(i,j, 0,k,w)*Delta A - ubar*y
!
!     is the streamfunction in a reference frame moving with the 
!     predetermined velocity ubar
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

use params
use mpi
implicit none
real*8 w(nx,ny)
real*8 psi(nx,ny)

! local
integer i,j,k,l,ierr
real*8 :: dela
real*8,external :: logterm
real*8 :: ubar=0
real*8 :: eps=1e-8
real*8 :: psi_b(max(g_ny,g_nx),2,2)   ! psi_b(k,1,1) = x1 boundary
                                      ! psi_b(k,1,2) = x2 boundary
                                      ! psi_b(k,2,1) = y1 boundary
                                      ! psi_b(k,2,2) = y2 boundary
real*8 :: temp(max(g_ny,g_nx),2,2)


! init to zero on boundary
psi_b=0

do j=ny1,ny2
do i=nx1,nx2
   if (abs(w(i,j)).ge.eps) then
      do k=1,g_nx  !,10
         psi_b(k,2,2) = psi_b(k,2,2) - w(i,j)*logterm(i,j,k,g_ny)
      enddo
      do k=2,g_ny  !,10
         psi_b(k,1,1) = psi_b(k,1,1) - w(i,j)*logterm(i,j,1,k)
         psi_b(k,1,2) = psi_b(k,1,2) - w(i,j)*logterm(i,j,g_nx,k)
      enddo
   endif
enddo
enddo

#if 0
C     INTERPOLATE TO INTERMEDIATE POINTS
      call intpsi(psi)
#endif


#ifdef USE_MPI
temp=psi_b
k=max(g_ny,g_nx)*4
call MPI_allreduce(temp,psi_b,k,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif


! add ubar correction
! only needs to be done on processer which owns boundary data
! but lets have everyone do it for now...
dela = delx*dely/(4*pi)
do k=1,g_nx
   ! y2 boundary
   psi_b(k,2,2) = psi_b(k,2,2)*dela - ubar*g_ycord(g_ny)
enddo
do k=2,g_ny
   ! x1,x2 boundary
   psi_b(k,1,1) = psi_b(k,1,1)*dela  - ubar*g_ycord(k)
   psi_b(k,1,2) = psi_b(k,1,2)*dela  - ubar*g_ycord(k)
enddo


! set PSI boundary data:

if (my_x==0 .and. bdy_x1==INFLOW0_ONESIDED) then
   do j=ny1,ny2
      l = j-ny1+1 + nslaby*my_y
      psi(nx1,j)= psi_b(L,1,1)
   enddo
endif
if (my_x==ncpu_x-1 .and. bdy_x2==INFLOW0_ONESIDED) then
   do j=ny1,ny2
      l = j-ny1+1 + nslaby*my_y
      psi(nx2,j)= psi_b(L,1,2)
   enddo
endif

if (my_y==0 .and. bdy_y1==INFLOW0_ONESIDED) then
   do i=nx1,nx2
      l = i-nx1+1 + nslabx*my_x
      psi(i,ny1)= psi_b(L,2,1)
   enddo
endif

if (my_y==ncpu_y-1   .and. bdy_y2==INFLOW0_ONESIDED) then
   do i=nx1,nx2
      l = i-nx1+1 + nslabx*my_x
      psi(i,ny2)= psi_b(L,2,2)
   enddo
endif






return
end subroutine



real*8 FUNCTION logterm(i,j,k,l)
!
!     Finds streamfunction at (k,l) induced by filament pair at (i,j)
!
use params
implicit none
integer i,j,k,l
real*8 difx,dify,sumy,denom1,denom2

difx = xcord(i) - g_xcord(k)
dify = ycord(j) - g_ycord(l)
sumy = ycord(j) + g_ycord(l)

denom1 = difx**2 + dify**2
denom2 = difx**2 + sumy**2

logterm = log(denom1/denom2)
end function 





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
real*8 :: ke,ens_diss,vor,dx,dxx,dy,dyy,ensave,u,v
integer n,i,j,k

call wallclock(tmx1)

ke=0
vor=0
ensave=0
ens_diss=0



rhs=0


do j=ny1,ny2
do i=nx1,nx2
   u=( 2*(psi(i,j+1)-psi(i,j-1))/3 -  &
        (psi(i,j+2)-psi(i,j-2))/12          )/dely

   v=-( 2*(psi(i+1,j)-psi(i-1,j))/3 -  &
        (psi(i+2,j)-psi(i-2,j))/12          )/delx


   if ((my_x==0 .and. bdy_x1==INFLOW0_ONESIDED .and. i<=nx1+2) .or. &
       (my_x==ncpu_x-1 .and. bdy_x2==INFLOW0_ONESIDED .and. i>=nx2-2)) then
      ! centered, 2nd order
      dx=( w(i+1,j)-w(i-1,j) ) / (2*delx)
      dxx=(w(i+1,j)-2*w(i,j)+w(i-1,j) ) / (delx*delx)
   else
      dx=( 2*(w(i+1,j)-w(i-1,j))/3 -  &
           (w(i+2,j)-w(i-2,j))/12          )/delx
      dxx=(-w(i+2,j) + 16*w(i+1,j) - 30*w(i,j) + &
           16*w(i-1,j) - w(i-2,j)) / (12*delx*delx)
   endif

   if ((my_y==0 .and. bdy_y1==INFLOW0_ONESIDED .and. j<=ny1+2) .or. &
       (my_y==ncpu_y-1 .and. bdy_y2==INFLOW0_ONESIDED .and. j>=ny2-2)) then
      ! centered, 2nd order
      dy=( w(i,j+1)-w(i,j-1) ) / (2*dely)
      dyy=(w(i,j+1)-2*w(i,j)+w(i,j-1) ) / (dely*dely)
    else
       dy=( 2*(w(i,j+1)-w(i,j-1))/3 -  &
            (w(i,j+2)-w(i,j-2))/12          )/dely
       dyy=(-w(i,j+2) + 16*w(i,j+1) - 30*w(i,j) + &
            16*w(i,j-1) - w(i,j-2)) / (12*dely*dely)
   endif



   rhs(i,j) = rhs(i,j) +  mu*(dxx+dyy) - u*dx - v*dy
   
   ke = ke + .5*(u**2 + v**2)

   vor=vor + w(i,j)
   ensave = ensave + w(i,j)**2
   ens_diss=ens_diss + w(i,j)*(dxx+dyy)
enddo
enddo






if (compute_ints==1) then
   !ints(3) = forcing terms
   ints(4)=vor/g_nx/g_ny
   ints(5)=ens_diss/g_nx/g_ny
   ints(6)=ke/g_nx/g_ny
   ints(7)=ensave /g_nx/g_ny
   ! ints(8) = < u,div(tau)' >   (alpha model only)
   ! ints(9)  = < u,f >  (alpha model only)
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
