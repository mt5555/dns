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

















subroutine compute_psi(w,psi,b,work,psi0)
!
!  btype=0   all periodic  use periodic FFT solver 
!  btype=1   all periodic,reflect, reflect-odd
!                use ghost cell CG solver (no boundaries)
!  btype=2   use biotsavart for boundary and dirichlet solver
!
use params
use ghost
implicit none
real*8 w(nx,ny,nz)
real*8 psi(nx,ny,nz)
real*8 work(nx,ny,nz)
real*8 b(nx,ny,nz)
real*8 psi0(nx,ny,nz)   ! initial guess, if needed

!local
real*8 :: one=1,zero=0,tol=1e-6
integer,save :: btype=-1

external helmholtz_dirichlet,helmholtz_periodic,helmholtz_periodic_ghost
integer i,j


if (btype==-1) then
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
   else if (bdy_x1==INFLOW0_ONESIDED .and. bdy_x2==INFLOW0_ONESIDED .and. &
            bdy_y1==REFLECT_ODD .and. bdy_y2==INFLOW0_ONESIDED) then
      btype=2
   else
      call abort("compute_psi(): boundery conditions not supported")
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
   call bc_biotsavart(w,psi)    !update PSI on boundary using biot-savart law

   b=-w  ! be sure to copy ghost cell data also!

   ! copy b.c. from psi into 'b', and apply compact correction to b:
   call helmholtz_dirichlet_setup(b,psi,work,1)
   call cgsolver(psi,b,zero,one,tol,work,helmholtz_dirichlet,.false.)


   !update PSI 1st row of ghost cells so that our 4th order differences
   !near the boundary look like 2nd order centered
   call bc_onesided(psi)
endif

call ghost_update_x(psi,1)
call ghost_update_y(psi,1)

end subroutine




subroutine helmholtz_dirichlet_inv(f,work,alpha,beta)
!
!  solve [alpha + beta*laplacian](p) = f
!  input:  f 
!  ouput:  f   will be overwritten with the solution p
!  b.c. for f specified in f. 
!
use params
use fft_interface
use transpose
implicit none
real*8 f(nx,ny,nz)    ! input/output
real*8 work(nx,ny,nz) ! work array
real*8 :: alpha
real*8 :: beta


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k
real*8 phi(nx,ny,nz)    
real*8 lphi(nx,ny,nz)    


!NOTE: we dont need phi, lphi.  when this code is debugged, replace by:
! save boundary data in single 1D array
! set f = 0 on boundary
! using boundary data alone: compute f = f - lphi 
!      along points just inside boundary
!
! solve as before
! restor boundary conditions in f from 1D array.




! phi = f on boundary, zero inside
phi=f
call zero_boundary(phi)
phi=f-phi
call helmholtz_dirichlet(phi,lphi,alpha,beta,work)


! convert problem to zero b.c. problem
f=f-lphi
call zero_boundary(f)
! solve Helm(f) = b with 0 b.c.
f=f+phi


end subroutine




subroutine helmholtz_periodic_inv(f,work,alpha,beta)
!
!  solve [alpha + beta*laplacian](p) = f
!  input:  f 
!  ouput:  f   will be overwritten with the solution p
!
use params
use fft_interface
use transpose
implicit none
real*8 f(nx,ny,nz)    ! input/output
real*8 work(nx,ny,nz) ! work array
real*8 :: alpha
real*8 :: beta


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k

if (beta==0) then
   f=f/alpha
   return
endif


call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d) 
call fft1(work,n1,n1d,n2,n2d,n3,n3d)     
call transpose_from_x(work,f,n1,n1d,n2,n2d,n3,n3d) 


call transpose_to_y(f,work,n1,n1d,n2,n2d,n3,n3d)  ! x,y,z -> y,x,z
call fft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(work,f,n1,n1d,n2,n2d,n3,n3d) 

call transpose_to_z(f,work,n1,n1d,n2,n2d,n3,n3d)  
call fft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(work,f,n1,n1d,n2,n2d,n3,n3d)  

! solve [alpha + beta*Laplacian] p = f.  f overwritten with output  p
call fft_laplace_inverse(f,alpha,beta)

call transpose_to_z(f,work,n1,n1d,n2,n2d,n3,n3d)       
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(work,f,n1,n1d,n2,n2d,n3,n3d)       

call transpose_to_y(f,work,n1,n1d,n2,n2d,n3,n3d)       
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(work,f,n1,n1d,n2,n2d,n3,n3d)         ! y,x,z -> x,y,z

call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d) 
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_x(work,f,n1,n1d,n2,n2d,n3,n3d )



end














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





