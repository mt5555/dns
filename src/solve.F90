#include "macros.h"

#define PRECOND


subroutine jacobi(sol,b,a1,a2,tol,h)
!
! Jacobi iteration for solving   [a1 + a2 Laplacian] x = b
!
! for shallow water equations, we are solving:
!                 [ a1 + 1/h a2 div h grad ] = b
! so we first multiply through by h to get a SPD system:
!                 [ h a1 + a2 div h grad ] = h b
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
real*8 :: h(nx,ny,nz)
real*8 :: tol,a1,a2

!local 
real*8 :: R(nx,ny,nz)
real*8 :: work(nx,ny,nz)

real*8,external  :: ddot

! local
integer i,j
real*8 err
real*8 res_init,tol1,w
integer itqmr
integer itermax

itermax=1000
w = a1 - a2*(pi2/(3*min(delx,dely,delz)))**2
w=2.5*w  ! safety factor
w = 1/w

if (equations==1) b=h*b

#ifdef PRECOND
call helmholtz_inv(b,work,a1,a2)
w=1
#endif



call helmholtz(sol,R,a1,a2,h)
#ifdef PRECOND
call helmholtz_inv(R,work,a1,a2)
#endif
R=b-R
sol = sol + w*R


res_init = sqrt(ddot(b,b))
if (tol>1) then
   itermax=tol
   tol1=1e-12
else
   tol1=tol*res_init
endif
   

itqmr=0

do 
   itqmr = itqmr+1
   
   call helmholtz(sol,R,a1,a2,h)
#ifdef PRECOND
   call helmholtz_inv(R,work,a1,a2)
#endif

   R=b-R
   sol = sol + w*R
   
   err=sqrt(ddot(R,R))
   !write(*,444) itqmr,res_init,err
   if (err <= tol1 .or. itqmr>=itermax) exit
enddo

      write(*,444) itqmr,res_init,err
444   format('Jacobi iter=',i4,' ||rhs||, residual=',e11.6,'  ',e11.6)

if (itqmr>=itermax) stop 'Jacobi iteration failure'   

end subroutine







subroutine cg(sol,b,a1,a2,tol)
real*8 :: dummy
call cg_shallow(sol,b,a1,a2,tol,dummy)
end subroutine






subroutine cg_shallow(sol,b,a1,a2,tol,h)
!
! CG for solving:
!   equations=0:    [a1 + a2 Laplacian] x = b
!
! for shallow water equations, we are solving:
!                 [ a1 + 1/h a2 div h grad ] = b
! so we first multiply through by h to get a SPD system:
!                 [ h a1 + a2 div h grad ] = h b
!
!
! input parameters:
! b   = rhs
! sol = starting approximate guess and final soln
! tol = tolerance for convergence (10e-6 is sufficient)
!
! if equations=1  (shallow water equations), then
! h = height field.  (and helmholtz = [a1 + a2 (1/h) div h grad] )
!
! routines called
!
! itermax = maximum number of allowed iterations
! helmholtz(x,y) :  performs  y=Ax
! ddot(x,y)       returns <x,y>
!
! preconditioner:   instead of Ax=b
! we solve  HAx=Hb  where H = helmholz_inv
!
use params
implicit none
real*8 :: h(nx,ny,nz)
real*8 :: sol(nx,ny,nz)
real*8 :: b(nx,ny,nz)
real*8 :: tol,a1,a2

!local 
real*8 :: P(nx,ny,nz)
real*8 :: R(nx,ny,nz)
real*8 :: AP(nx,ny,nz)
real*8 :: work(nx,ny,nz)

real*8,external  :: ddot

integer i,j
real*8 alpha,sigma,tau,beta,err,err_old
real*8 res_init,tol1,pa1,pa2
integer itqmr
integer itermax

if (equations==1) b=h*b


pa1=1
pa2=0
#ifdef PRECOND
pa1=a1
pa2=a2
call helmholtz_inv(b,work,pa1,pa2)
! use preconditioned value as initial guess also:
sol=b
#endif


itermax=1000
res_init=sqrt(ddot(sol,sol))          
if (tol>1) then
   itermax=tol
   tol1=1e-12
else
   tol1=tol*res_init
endif
err=1d99


itqmr=0

50    continue
      call helmholtz(sol,R,a1,a2,h)
#ifdef PRECOND
      call helmholtz_inv(R,work,pa1,pa2)
#endif
      R=b-R
      P=R
      alpha=ddot(R,R)

100   continue
      itqmr = itqmr+1
      sigma = alpha

      call helmholtz(P,AP,a1,a2,h)
#ifdef PRECOND
      call helmholtz_inv(AP,work,pa1,pa2)
#endif
      tau = sigma/ ddot(P,AP)

      sol=sol+P*tau
      R = R - AP*tau

      alpha = ddot(R,R)
      beta = alpha/sigma
      P=R + P*beta

      err_old=err
      err=sqrt(alpha)

      if (err.gt.tol1 .and. itqmr.lt.itermax) goto 100

      write(*,444) itqmr,res_init,err
444   format('CG iter=',i4,' ||RHS||, residual=',e11.6,'  ',e11.6)

if (itqmr>=itermax) then
   print *,'min/max',minval(h),maxval(h)
   stop 'CG iteration failure'   
endif
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







subroutine helmholtz(f,lf,alpha,beta,h)
!
!  input: f
!  output: lf
!
!  lf = [alpha + beta*laplacian](f)
!
use params
use fft_interface
use transpose
implicit none
real*8 f(nx,ny,nz)    ! input
real*8 lf(nx,ny,nz)    ! output
real*8 h(nx,ny,nz)    ! height field (used only for shallow water equations)
real*8 :: alpha
real*8 :: beta

!local
real*8 work(nx,ny,nz)    ! work array
real*8 gradf(nx,ny,nz,2) ! work array
real*8 fxx(nx,ny,nz)     ! work array
real*8 dummy
integer n

if (equations==1) then
   do n=1,2
      call der(f,gradf(1,1,1,n),dummy,work,DX_ONLY,n)
      gradf(:,:,:,n)=h*gradf(:,:,:,n)
   enddo
   call divergence(lf,gradf,fxx,work)
   lf=alpha*f*h + beta*lf
else
   ! linear helmholtz operator
   lf=alpha*f
   do n=1,3
      call der(f,gradf,fxx,work,DX_AND_DXX,n)
      lf=lf+beta*fxx
   enddo
endif

end subroutine

