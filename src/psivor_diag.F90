#include "macros.h"

logical function call_output_model()
!
!  The subroutine output_model() below is called every timestep
!  This routine will return .true. if it needs to be called,
!  and .false. if output_model() will do nothing if called
!
!  The point of this routine is for the super efficient ns_xpencil.F90 when
!  using an x-pencil decompostion.  That decomposition must be converted
!  back to our reference decompostion before output and these diagnostics
!  can be called.
!
!  failsafe version of call_output_model():  always return .true.
!  No extra costs for most models in the DNS code: only when using 
!  ns_xpencil.F90.   (and ns_xpencil.F90 is only necessary when running
!  with a pencil decomposition instead of slabs)
!  
call_output_model=.true.
return
end function


!
!  diagnostics for psi-vor model
!
subroutine output_model(doit_model,doit_diag,time,Q,Qhat,q1,q2,q3,work1,work2)
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
logical :: doit_model,doit_diag

if (.not.doit_model) return

call tracers_save(io_pe,time)

if (init_cond_subtype<100) then
   call comp_ellipse(Qhat(1,1,1,1),0,0)
   call ellipse_output(time)
endif

end subroutine




