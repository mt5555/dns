#include "macros.h"
!
! NOTE: before using an iterative solver, see the notes
! at the top of fftops.F90 about how to set
! DX4DX4 and HINV_DX4DX4
!
!
subroutine cgsolver(sol,b,a1,a2,tol,h,matvec,precond)
!
! CG for solving:
!   matvec(a1,a2,h) * sol = b
!
! Normally, matvec is the operator:  [a1 + a2 Laplacian]
! but for the shallow water equations, matvec is the operator:
!                 [ h a1 + a2 div h grad ] = b
! (note that we have multiplied through by h, so if you want to solve:
!     [a1 + 1/h a2 div h grad ] sol = b, then multiply 'b' by 'h'
!     before calling cgsolver)
!
!
! input parameters:
! b   = rhs
! sol = starting approximate guess and final solution
!       For dirichelt problems, sol and b most both be set
!       with the the dirichlet data on the boundary.
!       this is done by calling helmholtz_dirichlet_setup()
!
! tol = tolerance for convergence (10e-6 is sufficient)
!
! if equations=SHALLOW  (shallow water equations), then
! h = height field.  (and helmholtz = [a1 + a2 (1/h) div h grad] )
!
! routines called
!
! itermax = maximum number of allowed iterations
! matvec(x,y,a1,a2,h)  performs  y=Ax
! ddot(x,y)            returns <x,y>
!
! preconditioner:   if precond=.true., use a periodic left preconditioner.
! instead of Ax=b, we solve  HAx=Hb  where H = helmholz_periodic_inv
!
use params
implicit none
real*8 :: h(nx,ny,nz)
real*8 :: sol(nx,ny,nz)
real*8 :: b(nx,ny,nz)
real*8 :: tol,a1,a2
logical :: precond
external matvec

!local 
real*8 :: P(nx,ny,nz)
real*8 :: R(nx,ny,nz)
real*8 :: AP(nx,ny,nz)
real*8 :: work(nx,ny,nz)

real*8,external  :: ddot

integer i,j
real*8 alpha,sigma,tau,beta,err
real*8 res_init,tol1
integer itqmr
integer itermax

if (precond) then
   call helmholtz_periodic_inv(b,work,a1,a2)
   ! use preconditioned value as initial guess also:
   sol=b
endif

itermax=5000
res_init=sqrt(ddot(b,b))          
if (tol>1) then
   itermax=tol
   tol1=1e-12
else
   if (res_init<1) then
      tol1=tol
   else
      tol1=tol*res_init
   endif
endif
err=1d99

itqmr=0
call matvec(sol,R,a1,a2,h)
if (precond) call helmholtz_periodic_inv(R,work,a1,a2)
R=b-R
P=R
alpha=ddot(R,R)
err=sqrt(alpha)
if (err<tol1) goto 200

100 continue
   itqmr = itqmr+1
   sigma = alpha

   call matvec(P,AP,a1,a2,h)
   if (precond) call helmholtz_periodic_inv(AP,work,a1,a2)

   tau = sigma/ ddot(P,AP)

   sol=sol+P*tau
   R = R - AP*tau
   
   alpha = ddot(R,R)
   beta = alpha/sigma
   P=R + P*beta
   
   err=sqrt(alpha)
   !write(*,'(a,i4,a,e11.6,a,e11.6)')  'CG iter=',itqmr,' ||RHS||=',res_init,' residual: ',err

   if (err.gt.tol1 .and. itqmr.lt.itermax) goto 100
200 continue

if (itqmr>25 .and. (my_pe==io_pe)) then
   write(*,'(a,i4,a,e11.6,a,e11.6)')  'Final CG iter=',itqmr,' || RHS ||=',res_init,' residual: ',err
endif
   
if (itqmr>=itermax) then
   stop 'CG iteration failure'   
endif
end subroutine







subroutine jacobi(sol,b,a1,a2,tol,h,matvec,precond)
!
! Jacobi iteration.  See documentation for cgsolver() above
!
use params
implicit none
real*8 :: sol(nx,ny,nz)
real*8 :: b(nx,ny,nz)
real*8 :: h(nx,ny,nz)
real*8 :: tol,a1,a2
logical :: precond
external matvec

!local 
real*8 :: R(nx,ny,nz)
real*8 :: work(nx,ny,nz)

real*8,external  :: ddot

! local
integer i,j,k
real*8 err
real*8 res_init,tol1,w
integer itqmr
integer itermax

itermax=50000
w = a1 - a2*(pi2/(3*min(delx,dely,delz)))**2
w=2.5*w  ! safety factor
w = 1/w

if (precond) then
   call helmholtz_periodic_inv(b,work,a1,a2)
   ! use preconditioned value as initial guess also:
   sol=b
   w=1
endif


call matvec(sol,R,a1,a2,h)
if (precond) call helmholtz_periodic_inv(R,work,a1,a2)

R=b-R
sol = sol + w*R

res_init = sqrt(ddot(b,b))


if (tol>1) then
   itermax=tol
   tol1=1e-12
else
   if (res_init<1) then
      tol1=tol
   else
      tol1=tol*res_init
   endif
endif
   

itqmr=0

do 
   itqmr = itqmr+1
   
   call matvec(sol,R,a1,a2,h)
   if (precond) call helmholtz_periodic_inv(R,work,a1,a2)

   R=b-R
   sol = sol + w*R
   
   err=sqrt(ddot(R,R))
   if (mod(itqmr,25)==0) &
      write(*,444) itqmr,res_init,err
   if (err <= tol1 .or. itqmr>=itermax) exit
enddo

      write(*,444) itqmr,res_init,err
444   format('Jacobi iter=',i8,' ||rhs||, residual=',e11.6,'  ',e11.6)

if (itqmr>=itermax) then
   print *,'Jacobi iteration failure...'   
endif

end subroutine

















subroutine compute_psi(w,psi,b,work,psi0,runbs)
!
!  btype=0   all periodic  use periodic FFT solver 
!  btype=1   all periodic,reflect, reflect-odd
!                use ghost cell CG solver (no boundaries)
!  btype=2   use biotsavart for boundary and dirichlet solver
!
use params
use ghost
use bc
implicit none
real*8 w(nx,ny,nz)
real*8 psi(nx,ny,nz)
real*8 work(nx,ny,nz)
real*8 b(nx,ny,nz)
real*8 psi0(nx,ny,nz)   ! initial guess, if needed
integer :: runbs 

!local
real*8 :: one=1,zero=0,tol=1e-10
integer,save :: btype=-1
real*8 :: tmx1,tmx2



external helmholtz_dirichlet,helmholtz_periodic,helmholtz_periodic_ghost
integer i,j

call wallclock(tmx1)

! check global boundary conditions and see what kind of solver we need
if (btype==-1) then
   if (g_bdy_x1==PERIODIC .and. g_bdy_y1==PERIODIC) then
      ! periodic
      btype=0
   else if  ( &  ! not periodic, but can still do with all ghostcells:
        ((g_bdy_x1==PERIODIC) .or. (g_bdy_x1==REFLECT) .or. (g_bdy_x1==REFLECT_ODD)) .and. &
        ((g_bdy_x2==PERIODIC) .or. (g_bdy_x2==REFLECT) .or. (g_bdy_x2==REFLECT_ODD)) .and. &
        ((g_bdy_y1==PERIODIC) .or. (g_bdy_y1==REFLECT) .or. (g_bdy_y1==REFLECT_ODD)) .and. &
        ((g_bdy_y2==PERIODIC) .or. (g_bdy_y2==REFLECT) .or. (g_bdy_y2==REFLECT_ODD))  ) then
        btype=1
   else if (REALBOUNDARY(g_bdy_x1)  .and. REALBOUNDARY(g_bdy_x2) .and. &
            g_bdy_y1==REFLECT_ODD .and. REALBOUNDARY(g_bdy_y2)) then
      btype=2
   else
      ! other types of b.c. not supported, mostly because bc_biotsavart
      ! is hard-coded for btype=2
      call abortdns("compute_psi(): specified boundery conditions not supported")
   endif
endif



! if all b.c. periodic:
if (btype==0) then
   psi=-w
   call helmholtz_periodic_inv(psi,work,zero,one)

   !psi=0  ! initial guess
   !b=-w
   !call cgsolver(psi,b,zero,one,tol,work,helmholtz_periodic,.true.)

else if (btype==1) then
   ! this code can handle any compination of PERIODIC, REFLECT, REFLECT-ODD:

   psi=psi0  ! initial guess
   b=-w  ! be sure to copy ghost cell data also!
   ! apply compact correction to b:
   call helmholtz_dirichlet_setup(b,psi,work,0)
   call cgsolver(psi,b,zero,one,tol,work,helmholtz_periodic_ghost,.false.)

else if (btype==2) then

   psi=0  ! initial guess
   call bc_biotsavart(w,psi,runbs)    !update PSI on boundary using biot-savart law
   b=-w  ! be sure to copy ghost cell data also!

   ! apply compact correction to 'b', then set b.c. of b:
   call helmholtz_dirichlet_setup(b,psi,work,1)
   psi=b; call helmholtz_dirichlet_inv(psi,work,zero,one) 
   !call cgsolver(psi,b,zero,one,tol,work,helmholtz_dirichlet,.false.)

   !update PSI 1st row of ghost cells so that our 4th order differences
   !near the boundary look like 2nd order centered
   call bc_onesided(psi)
endif



call ghost_update_x(psi,1)
call ghost_update_y(psi,1)

      
call wallclock(tmx2)
tims(15)=tims(15)+(tmx2-tmx1)


end subroutine












real*8 function ddot(a,b)
use params
use mpi
implicit none
real*8 :: a(nx,ny,nz)
real*8 :: b(nx,ny,nz)

! local
integer :: i,j,k ,ierr
real*8 :: sm

ddot=0
do k=bz1,bz2
do j=by1,by2
do i=bx1,bx2
   ddot=ddot+a(i,j,k)*b(i,j,k)
enddo
enddo
enddo

#ifdef USE_MPI
sm=ddot
call mpi_allreduce(sm,ddot,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif

return
end function





