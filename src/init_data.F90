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
   !call set_byteswap_input(1);
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
if (equations==NS_UVW) then

if (udm_input) then
   call udm_read_uvw(Q,Qhat,work1,work2)
else
   if (rw_spec) then
   call print_message("Restarting from file restart.[us,vs,ws]")
   fname = rundir(1:len_trim(rundir)) // "restart.us"
   call singlefile_io2(time_initial,Q(1,1,1,1),fname,work1,work2,1,io_pe,rw_spec)
   fname = rundir(1:len_trim(rundir)) // "restart.vs"
   call singlefile_io2(time_initial,Q(1,1,1,2),fname,work1,work2,1,io_pe,rw_spec)
   if (n_var==3) then
      fname = rundir(1:len_trim(rundir)) // "restart.ws"
      call singlefile_io2(time_initial,Q(1,1,1,n_var),fname,work1,work2,1,io_pe,rw_spec)
   endif
   do n=1,ndim
      call ifft3d(Q(1,1,1,n),work1)
   enddo
   else
   call print_message("Restarting from file restart.[uvw]")
   fname = rundir(1:len_trim(rundir)) // "restart.u"
   call singlefile_io2(time_initial,Q(1,1,1,1),fname,work1,work2,1,io_pe,rw_spec)
   fname = rundir(1:len_trim(rundir)) // "restart.v"
   call singlefile_io2(time_initial,Q(1,1,1,2),fname,work1,work2,1,io_pe,rw_spec)
   if (n_var==3) then
      fname = rundir(1:len_trim(rundir)) // "restart.w"
      call singlefile_io2(time_initial,Q(1,1,1,n_var),fname,work1,work2,1,io_pe,rw_spec)
   endif
   endif
endif

else if (equations==SHALLOW) then
   call print_message("Restarting from file restart.[uvh]")
   fname = rundir(1:len_trim(rundir)) // "restart.u"
   call singlefile_io(time_initial,Q(1,1,1,1),fname,work1,work2,1,io_pe)
   fname = rundir(1:len_trim(rundir)) // "restart.v"
   call singlefile_io(time_initial,Q(1,1,1,2),fname,work1,work2,1,io_pe)
   if (n_var==3) then
      fname = rundir(1:len_trim(rundir)) // "restart.h"
      call singlefile_io(time_initial,Q(1,1,1,n_var),fname,work1,work2,1,io_pe)
   endif
else if (equations==NS_PSIVOR) then
   call print_message("Restarting from file restart.vor")
   call abort("restart error: check it is ok that we setup I.C. before reading data")
   fname = rundir(1:len_trim(rundir)) // "restart.vor"
   call singlefile_io(time_initial,Qhat(1,1,1,1),fname,work1,work2,1,io_pe)
   !fname = rundir(1:len_trim(rundir)) // "restart.psi"
   !call singlefile_io(time_initial,Qhat,fname,work1,work2,1,io_pe)
endif



write(message,'(a,f10.4)') "restart time=",time_initial
call print_message(message)
end subroutine









