#include "macros.h"
subroutine init_data(Q,Qhat,work1,work2)
use params
use mpi
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
if (restart==1) then
#ifdef BYTESWAP_INPUT
   call set_byteswap_input(1);
#endif
   if (byteswap_input) call set_byteswap_input(1);

   ! initialize some constants, if needed on restart runs:
   if (init_cond==3) call init_data_sht(Q,Qhat,work1,work2,0)      ! set grav, fcor
   if (init_cond==4) call init_data_vxpair(Q,Qhat,work1,work2,0) ! set xscale, yscale... 

   call init_data_restart(Q,Qhat,work1,work2)

   if (init_cond==9) call init_data_decay(Q,Qhat,work1,work2,2,0,0)

else
   if (init_cond==0) call init_data_khblob(Q,Qhat,work1,work2)
   if (init_cond==1) call init_data_kh(Q,Qhat,work1,work2)
   if (init_cond==2) call init_data_lwisotropic(Q,Qhat,work1,work2,1,0)
   if (init_cond==3) call init_data_sht(Q,Qhat,work1,work2,1)
   if (init_cond==4) call init_data_vxpair(Q,Qhat,work1,work2,1)
   if (init_cond==5) call init_data_lwisotropic(Q,Qhat,work1,work2,1,1)
   if (init_cond==6) call init_data_zero(Q,Qhat,work1,work2)
   if (init_cond==7) call init_data_decay(Q,Qhat,work1,work2,1,0,0)
   if (init_cond==8) call init_data_decay(Q,Qhat,work1,work2,1,1,0)
endif
end subroutine


subroutine init_data_restart(Q,Qhat,work1,work2)
!
! low wave number, quasi isotropic initial condition
!
use params
use mpi
use fft_interface
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

!local
character(len=80) message
character(len=80) fname
integer :: n

Q=0
time_initial=-1
call input_uvw(time_initial,Q,Qhat,work1,work2)

write(message,'(a,f10.4)') "restart time=",time_initial
call print_message(message)
end subroutine









