#include "macros.h"
subroutine jacobi(sol,b,a1,a2,tol)
!
! Jacobi iteration for solving   [a1 + a2 Laplacian] x = b
!
! input parameters:
! b   = rhs
! sol = starting approximate guess and final soln
! tol = tolerance for convergence (10e-6 is sufficient)
!
! routines called
!
! itermax = maximum number of allowed iterations
! helmholtz(x,y) :  performs  y=Ax
! ddot(x,y)       returns <x,y>
!
use params
implicit none
real*8 :: sol(nx,ny,nz)
real*8 :: b(nx,ny,nz)
real*8 :: tol,a1,a2

!local 
real*8 :: R(nx,ny,nz)

real*8,external  :: ddot

! local
integer i,j
real*8 err
real*8 res_init,tol1,w
integer itqmr
integer itermax

itermax=1000
w = a1 - a2*(pi2/(3*min(delx,dely,delz)))**2
w=2.5*w  ! safter factor
w = 1/w



call helmholtz(sol,R,a1,a2)
R=b-R
sol = sol + w*R


res_init = sqrt(ddot(R,R))
if (tol>1) then
   itermax=tol
   tol1=1e-12
else
   tol1=tol*res_init
endif
   
if (res_init<1e-20) return


itqmr=0

do 
   itqmr = itqmr+1
   
   call helmholtz(sol,R,a1,a2)
   R=b-R
   sol = sol + w*R
   
   err=sqrt(ddot(R,R))
   !write(*,444) itqmr,res_init,err
   if (err <= tol1 .or. itqmr>=itermax) exit
enddo

      write(*,444) itqmr,res_init,err
444   format('Jacobi iter=',i4,' res_init, residual=',e11.6,'  ',e11.6)

if (itqmr>=itermax) stop 'Jacobi iteration failure'   

end subroutine






subroutine cg(sol,b,a1,a2,tol)
!
! CG for solving   [a1 + a2 Laplacian] x = b
!
! input parameters:
! b   = rhs
! sol = starting approximate guess and final soln
! tol = tolerance for convergence (10e-6 is sufficient)
!
! routines called
!
! itermax = maximum number of allowed iterations
! helmholtz(x,y) :  performs  y=Ax
! ddot(x,y)       returns <x,y>
!
use params
implicit none
real*8 :: sol(nx,ny,nz)
real*8 :: b(nx,ny,nz)
real*8 :: tol,a1,a2

!local 
real*8 :: P(nx,ny,nz)
real*8 :: R(nx,ny,nz)
real*8 :: AP(nx,ny,nz)

real*8,external  :: ddot

! local
integer i,j
real*8 alpha,sigma,tau,beta,err
real*8 res_init,tol1
integer itqmr
integer itermax

itermax=1000


      call helmholtz(sol,R,a1,a2)
      R=b-R
      P=R

      alpha=ddot(R,R)
      res_init = sqrt(alpha)
      if (tol>1) then
         itermax=tol
         tol1=1e-12
      else
         tol1=tol*res_init
      endif


      if (res_init<1e-20) return

      itqmr=0
100   continue
      itqmr = itqmr+1
      sigma = alpha

      call helmholtz(P,AP,a1,a2)
      tau = sigma/ ddot(P,AP)

      sol=sol+P*tau
      R = R - AP*tau

      alpha = ddot(R,R)
      beta = alpha/sigma
      P=R + P*beta

      err=sqrt(alpha)
      if (err.gt.tol1 .and. itqmr.lt.itermax) goto 100

      write(*,444) itqmr,res_init,err
444   format('CG iter=',i4,' res_init, residual=',e11.6,'  ',e11.6)

if (itqmr>=itermax) stop 'CG iteration failure'   
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





