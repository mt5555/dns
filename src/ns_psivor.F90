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



subroutine rk4(time,Q,Qsave,Q2,Q3,rhs,work1,work2)
use params
implicit none
real*8 :: time
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qsave(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: Q2(nx,ny,nz,n_var)
real*8 :: Q3(nx,ny,nz,n_var)
real*8 :: rhs(nx,ny,nz,n_var)

call rk4reshape(time,Q,Qsave(1,1,1,1),Qsave(1,1,1,2),Q2(1,1,1,1),Q2(1,1,1,2),&
Q3(1,1,1,1),work1,work2)
end subroutine


subroutine rk4reshape(time,Q,w0,psi0,rhs,w_old,w_tmp,psi,work)
use params
use bc
use tracers
use ellipse
implicit none
real*8 :: time
real*8 :: Q(nx,ny,n_var)
real*8 :: w0(nx,ny)
real*8 :: psi0(nx,ny)
real*8 :: w_tmp(nx,ny)
real*8 :: w_old(nx,ny)
real*8 :: psi(nx,ny)
real*8 :: rhs(nx,ny)
real*8 :: work(nx,ny)


! local variables
real*8 :: vel,u,v
integer i,j,k,n,ierr
integer,save :: ncall=0
integer,save :: btype

integer :: comp_psi0 = 0       ! compute psi on boundary at time=0 
integer :: comp_psi_rk13=0     ! compute psi on boundary after RK stages 1..3?
integer :: comp_psi_rk4=0      ! compute psi on boundary after RK stage 4?

ncall=ncall+1


if (biotsavart_apply>0) then
   comp_psi0=1
   comp_psi_rk4=0
   if (mod(ncall,biotsavart_apply)==0) comp_psi_rk4=1
else
   ! do not update psi on boundary
endif




if (ncall==1) then
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

   ! set boundary data
   call bcw_impose(w0)
   w_tmp=0
   psi0=0
   ! initialize PSI on interior.  
   ! recompute PSI on boundary if comp_psi0==1
   call compute_psi(w0,psi0,rhs,work,w_tmp,comp_psi0)
endif




!
!  Q(:,:,:,1) = u
!  Q(:,:,:,2) = v
!  Q(:,:,:,3) = vor
!


! find the location of peak vorticity (needed by tracer_advance)
call comp_ellipse(w0,0,1)

w_old=w0


! stage 1
call ns2D(rhs,w0,psi0,time,1); call tracer_advance(psi0,Q,1,time)
w0=w0+delt*rhs/6.0 


! stage 2
w_tmp = w_old + delt*rhs/2.0
call bcw_impose(w_tmp)
call compute_psi(w_tmp,psi,rhs,work,psi0,comp_psi_rk13)
call ns2D(rhs,w_tmp,psi,time+delt/2.0,0); call tracer_advance(psi,Q,2,time)
w0=w0+delt*rhs/3.0


! stage 3
w_tmp = w_old + delt*rhs/2.0
call bcw_impose(w_tmp)
call compute_psi(w_tmp,psi,rhs,work,psi0,comp_psi_rk13)
call ns2D(rhs,w_tmp,psi,time+delt/2.0,0) ; call tracer_advance(psi,Q,3,time)
w0=w0+delt*rhs/3.0



! stage 4
w_tmp = w_old + delt*rhs
call bcw_impose(w_tmp)
call compute_psi(w_tmp,psi,rhs,work,psi0,comp_psi_rk13)
call ns2D(rhs,w_tmp,psi,time+delt,0) ; call tracer_advance(psi,Q,4,time)
w0=w0+delt*rhs/6.0
call bcw_impose(w0)



call compute_psi(w0,psi0,rhs,work,psi,comp_psi_rk4)
time = time + delt



! compute KE, max U  
maxs(1:4)=0
do j=inty1,inty2
do i=intx1,intx2
      
      u=( 2*(psi0(i,j+1)-psi0(i,j-1))/3 -  &
           (psi0(i,j+2)-psi0(i,j-2))/12          )/dely
      
      v=-( 2*(psi0(i+1,j)-psi0(i-1,j))/3 -  &
           (psi0(i+2,j)-psi0(i-2,j))/12          )/delx
      
      maxs(1)=max(maxs(1),abs(u))
      maxs(2)=max(maxs(2),abs(v))
      maxs(3)=0
      vel = abs(u)/delx + abs(v)/dely 
      maxs(4)=max(maxs(4),vel)
      Q(i,j,1)=u
      Q(i,j,2)=v
enddo
enddo



end subroutine










subroutine ns2D(rhs,w,psi,time,compute_ints)
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
real*8 :: ke,ens_diss,vor,dx,dxx,dy,dyy,ensave,u,v,maxvor
integer n,i,j,k

call wallclock(tmx1)

ke=0
vor=0
ensave=0
ens_diss=0
maxvor=0


rhs=0


do j=inty1,inty2
do i=intx1,intx2
   u=( 2*(psi(i,j+1)-psi(i,j-1))/3 -  &
        (psi(i,j+2)-psi(i,j-2))/12          )/dely

   v=-( 2*(psi(i+1,j)-psi(i-1,j))/3 -  &
        (psi(i+2,j)-psi(i-2,j))/12          )/delx


   if ( (REALBOUNDARY(bdy_x1) .and. i<=bx1+2) .or. &
        (REALBOUNDARY(bdy_x2) .and. i>=bx2-2) ) then
      ! centered, 2nd order
      dx=( w(i+1,j)-w(i-1,j) ) / (2*delx)
      dxx=(w(i+1,j)-2*w(i,j)+w(i-1,j) ) / (delx*delx)
   else
      dx=( 2*(w(i+1,j)-w(i-1,j))/3 -  &
           (w(i+2,j)-w(i-2,j))/12          )/delx
      dxx=(-w(i+2,j) + 16*w(i+1,j) - 30*w(i,j) + &
           16*w(i-1,j) - w(i-2,j)) / (12*delx*delx)
   endif


   if ( (REALBOUNDARY(bdy_y1) .and. j<=by1+2) .or. &
        (REALBOUNDARY(bdy_y2) .and. j>=by2-2)) then
      ! centered, 2nd order
      dy=( w(i,j+1)-w(i,j-1) ) / (2*dely)
      dyy=(w(i,j+1)-2*w(i,j)+w(i,j-1) ) / (dely*dely)
    else
       dy=( 2*(w(i,j+1)-w(i,j-1))/3 -  &
            (w(i,j+2)-w(i,j-2))/12          )/dely
       dyy=(-w(i,j+2) + 16*w(i,j+1) - 30*w(i,j) + &
            16*w(i,j-1) - w(i,j-2)) / (12*dely*dely)
   endif

   rhs(i,j) = rhs(i,j) - u*dx - v*dy + mu*(dxx + dyy)

   if (compute_ints==1) then
         maxvor=max(maxvor,abs(w(i,j)))
         ke = ke + .5*(u**2 + v**2)
         vor=vor + w(i,j)
         ensave = ensave + w(i,j)**2
         ens_diss=ens_diss + w(i,j)*mu*(dxx+dyy)
   endif

enddo
enddo


if (compute_ints==1) then
   !ints(3) = forcing terms
   ints(4)=vor/o_nx/o_ny
   ints(5)=ens_diss/o_nx/o_ny
   ints(6)=ke/o_nx/o_ny
   ints(7)=ensave /o_nx/o_ny
   ! ints(8) = < u,div(tau)' >   (alpha model only)
   ! ints(9)  = < u,f >  (alpha model only)

   maxs(5)=maxvor
endif



call wallclock(tmx2)
tims(5)=tims(5)+(tmx2-tmx1)
end














