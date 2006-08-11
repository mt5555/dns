#include "macros.h"


function call_output_model()
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
return .true.
end function


!
!  diagnostics for shallow water model
!
subroutine output_model(doit_model,doit_diag,time,Q,Qhat,q1,q2,q3,work1,work2)
use params
use spectrum
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(*)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: q2(nx,ny,nz,n_var)
real*8 :: q3(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time
logical :: doit_model,doit_diag

real*8 :: tsave

if (compute_transfer) then
   ! spec_r computed last time step
   ! spec_diff,spec_model were computed in RHS computation at the
   ! beginning of this flag (becuase compute_transfer flag was set)
   ! So they are all known at time_old. Now compute spec_r_new 
   ! (used to compute edot_r at mid-time level)
   call compute_Edotspec_shallow(time,Q,q1,work1,work2)

   ! e_dot_r    known at (time + time_old)/2
   ! spec_diff,spec_model  known at time_old
   ! compute spec_diff at 'time' be re-calling getrhs with compute_transfer
   ! flag still set = .true. and compute_ints==1

   spec_diff_new=spec_diff  
   spec_model_new=spec_model
   tsave=transfer_comp_time
   call getrhs(4,q1,Qhat,Q,time,1,work1,work2)
   transfer_comp_time=tsave

   spec_diff = (spec_diff_new+spec_diff)/2
   spec_model = (spec_model_new+spec_model)/2
   
   ! output all the spectrum:
   call output_tran(time,Q,q1,q2,q3,work1,work2)


   compute_transfer=.false.
endif


if (doit_model .or. time==time_initial) then
if ( g_bdy_x1==PERIODIC .and. &
     g_bdy_y1==PERIODIC .and. &
     g_bdy_z1==PERIODIC) then
   call compute_spec_shallow(time,Q,q1,work1,work2)
   call output_spec(time,time_initial)
   
   !set this flag so that for next timestep, we will compute and save
   !spectral transfer functions:
   compute_transfer=.true.
endif
endif


end subroutine
