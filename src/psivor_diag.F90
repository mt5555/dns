#include "macros.h"
!
!  diagnostics for psi-vor model
!
subroutine output_model(time,Q,Qhat,q1,q2,q3,work1,work2)
use params
use ellipse
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: q2(nx,ny,nz,n_var)
real*8 :: q3(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time

call ellipse_output()


end subroutine




subroutine comp_ellipse_reshape(w,setmax)
use params
use ellipse
implicit none
real*8 :: w(nx,ny)
integer :: setmax
call comp_ellipse(w,setmax)
end
