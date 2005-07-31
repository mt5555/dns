#include "macros.h"
subroutine output_model(doit_model,doit_diag,time,Q,Qhat,q1,q2,q3,work1,work2)
use params
use spectrum
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(g_nz2,nslabx,ny_2dz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: q2(nx,ny,nz,n_var)
real*8 :: q3(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time
logical :: doit_model,doit_diag

! local variables
integer,parameter :: nints_e=49,npints_e=51
real*8 :: ints_e(nints_e)
real*8 :: pints_e(npints_e,n_var)
real*8 :: x,zero_len
real*8 :: divx,divi
integer i,j,k,n,ierr,csig
integer :: n1,n1d,n2,n2d,n3,n3d
character(len=80) :: message
CPOINTER fid,fidj,fidS,fidcore




if (compute_transfer) then
   compute_transfer=.false.
   ! spec_r computed last time step
   ! spec_diff, spec_f, spec_rhs were computed in RHS computation at the
   ! beginning of this flag (becuase compute_transfer flag was set)
   ! So they are all known at time_old. Now compute spec_r_new 
   ! (used to compute edot_r)
   call compute_Edotspec(time,Q,q1,work1,work2)
   ! output all the spectrum:
   call output_tran(time,Q,q1,q2,q3,work1,work2)
endif


! compute spectrum
! always compute at first timestep because transfer cannot be computed
! on last timestep.   
if (doit_model .or. time==time_initial) then
if ( g_bdy_x1==PERIODIC .and. &
     g_bdy_y1==PERIODIC .and. &
     g_bdy_z1==PERIODIC) then


   call compute_spec(time,Q,q1,work1,work2)
   call compute_spec_2d(time,Q,q1,work1,work2)
   call output_spec(time,time_initial)
   call output_helicity_spec(time,time_initial)  ! put all hel spec in same file
   call output_2d_spec(time,time_initial)  

   !set this flag so that for next timestep, we will compute and save
   !spectral transfer functions:
   compute_transfer=.true.


   ! for incompressible equations, print divergence as diagnostic:
   if (equations==NS_UVW) then
      call compute_div(Q,q1,work1,work2,divx,divi)
      write(message,'(3(a,e12.5))') 'max(div)=',divx
      call print_message(message)	
   endif


endif
endif

! do PDF's and scalars if doit_model=.true., OR if this is a restart
! but we have computed new passive scalars.
if ((compute_passive_on_restart .and. time==time_initial) .or. &
    doit_model) then
   ! do the rest of this suburoutine
else
   return
endif






end subroutine



