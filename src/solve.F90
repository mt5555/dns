#include "macros.h"
!
! NOTE: before using an iterative solver, see the nots
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
external :: matvec

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


itermax=1000
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
   write(*,'(a,i4,a,e11.6,a,e11.6)')  'CG iter=',itqmr,' ||RHS||=',res_init,' residual: ',err

if (err.gt.tol1 .and. itqmr.lt.itermax) goto 100
200 continue

if (itqmr>25) then
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
external :: matvec

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

itermax=2000
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
   write(*,444) itqmr,res_init,err
   if (err <= tol1 .or. itqmr>=itermax) exit
enddo

      write(*,444) itqmr,res_init,err
444   format('Jacobi iter=',i4,' ||rhs||, residual=',e11.6,'  ',e11.6)

if (itqmr>=itermax) stop 'Jacobi iteration failure'   

end subroutine


















real*8 function ddot(a,b)
use params
use mpi
implicit none
real*8 :: a(nx,ny,nz)
real*8 :: b(nx,ny,nz)
real*8 :: sm

! local
integer :: i,j,k ,ierr

ddot=0
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   ddot=ddot+a(i,j,k)*b(i,j,k)
enddo
enddo
enddo

#ifdef USE_MPI
sm=ddot
call MPI_allreduce(sm,ddot,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif

return
end function





