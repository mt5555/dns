#include "macros.h"
!
!  diagnostics for psi-vor model
!
subroutine output_model(doit_model,time,Q,Qhat,q1,q2,q3,work1,work2)
use params
use ellipse
use tracers
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: q2(nx,ny,nz,n_var)
real*8 :: q3(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time
logical :: doit_model

if (.not.doit_model) return

call tracers_save(io_pe,time)
call comp_ellipse(Qhat(1,1,1,1),0,0)
call ellipse_output(time)


end subroutine




