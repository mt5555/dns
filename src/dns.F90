#include "macros.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Fully resolved compressible Navier Stokes code
! copyright 2001 by Hal Marshal, Mark Taylor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program DNS
use params
use mpi
implicit none
real*8,save :: Q(nx,ny,nz,n_var)
real*8,save :: Qhat(nx,ny,nz,n_var)
real*8,save :: q1(nx,ny,nz,n_var)
real*8,save :: work1(nx,ny,nz)
real*8,save :: work2(nx,ny,nz)
character(len=80) message
integer ierr,i,j,k,n
real*8 tmx1,tmx2,tims_max(ntimers),tims_ave(ntimers)


call init_mpi       
call wallclock(tmx1)  ! wallclock may use MPI timers, call after init_mpi
call init_mpi_comm3d()

call init_model()

!write(message,'(a)') 'Running some tests'
!call print_message(message)
!call test           ! optional testing  routines go here
!stop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Set initial data, make it divergence free
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(message,'(a)') 'Initial data'
call print_message(message)
Q=0
call init_data(Q,Qhat,work1,work2)

if (equations==NS_UVW) then
   call print_message('Projecting initial data...')
   call divfree_gridspace(Q,work1,work2,q1) 
else if (equations==SHALLOW .or. equations==NS_PSIVOR) then
   if (dealias)  then
      call print_message('Dealiasing initial data...')
      call dealias_gridspace(Q,work1)
   endif
endif


#ifdef USE_MPI
call MPI_Barrier(comm_3d,ierr)
#endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Solve the equations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call wallclock(tmx2)
write(message,'(a)') 'Time stepping...'
call print_message(message)
! set all counters to zero (so they dont include initialization)
tims=0
! tims(1) times the total initialization
tims(1)=tmx2-tmx1

call dns_solve(Q,Qhat,q1,work1,work2)






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Print timing information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
tims_max=tims
tims_ave=tims
#ifdef USE_MPI
   call MPI_allreduce(tims,tims_max,ntimers,MPI_REAL8,MPI_MAX,comm_3d,ierr)
   call MPI_allreduce(tims,tims_ave,ntimers,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   tims_ave=tims_ave/initial_live_procs
#endif
tims_max=tims_max/60
tims_ave=tims_ave/60



call print_message('CPU times (min):   (avg/max)')
write(message,'(a,2f9.2,a)') 'initialization: ',tims_ave(1),tims_max(1)
call print_message(message)
write(message,'(a,2f9.2,a)') 'dns_solve:      ',tims_ave(2),tims_max(2)
call print_message(message)
write(message,'(a,2f9.2,a,f4.3,a)') '   time_control       ',tims_ave(3),tims_max(3)
call print_message(message)
write(message,'(a,2f9.2,a,f4.3,a)') '   RHS                ',tims_ave(5),tims_max(5)
call print_message(message)
if (tims_max(12)>0) then
   write(message,'(a,2f9.2,a,f4.3,a)') '   compute_pdf        ',tims_ave(12),tims_max(12)
   call print_message(message)
endif
if (tims_max(13)>0) then
   write(message,'(a,2f9.2,a,f4.3,a)') '   ghost_update       ',tims_ave(13),tims_max(13)
   call print_message(message)
endif
if (maxval(tims_max(6:11))>0) then
   write(message,'(a,2f9.2,a,f4.3,a)') '   transpose_to_z     ',tims_ave(6),tims_max(6)
   call print_message(message)
   write(message,'(a,2f9.2,a,f4.3,a)') '   transpose_from_z   ',tims_ave(7),tims_max(7)
   call print_message(message)
   write(message,'(a,2f9.2,a,f4.3,a)') '   transpose_to_x     ',tims_ave(8),tims_max(8)
   call print_message(message)
   write(message,'(a,2f9.2,a,f4.3,a)') '   transpose_from_x   ',tims_ave(9),tims_max(9)
   call print_message(message)
   write(message,'(a,2f9.2,a,f4.3,a)') '   transpose_to_y     ',tims_ave(10),tims_max(10)
   call print_message(message)
   write(message,'(a,2f9.2,a,f4.3,a)') '   transpose_from_y   ',tims_ave(11),tims_max(11)
   call print_message(message)
endif
if (tims_max(14)>0) then
   write(message,'(a,2f9.2,a,f4.3,a)') '   Biot-Savart        ',tims_ave(14),tims_max(14)
   call print_message(message)
endif
if (tims_max(15)>0) then
   write(message,'(a,2f9.2,a,f4.3,a)') '   compute_psi        ',tims_ave(15),tims_max(15)
   call print_message(message)
endif


call close_mpi


end program DNS









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  main time stepping loop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dns_solve(Q,Qhat,q1,work1,work2)
use params
use mpi

implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

!local variables
real*8,save :: q2(nx,ny,nz,n_var)
real*8,save :: q3(nx,ny,nz,n_var)


real*8  :: time=0
integer :: itime=0,ierr,n
integer :: itime_final=2**30
character(len=80) message
real*8 :: time_old=0
real*8 :: ea_new=0,ea_old
real*8 :: ke_new=0
real*8 :: ints_buf(nints)
real*8 :: tmx1,tmx2
integer,external :: lsf_time_remaining
integer :: lsftime



ints=0
maxs=0
delt=0

time=time_initial
if (time_final<0) then
   ! flag denoting run abs(time_final) time steps, instead of specifing 
   ! a final time
   itime_final=-time_final
   time_final=1e20
else
   time_final=time_final+time_initial
endif
write(message,'(a,f9.4)') "initial time: ",time_initial
call print_message(message)
write(message,'(a,f9.4)') "final time:   ",time_final
call print_message(message)
write(message,'(a,i10)') "max number timesteps: ",itime_final
call print_message(message)



do 
   time_old=time

   call rk4(time,Q,Qhat,q1,q2,q3,work1,work2)

   maxs(8)=-1  
   if (my_pe==io_pe) then
      if (lsf_time_remaining(lsftime)==0) then
         maxs(8)=lsftime/60.0
      endif
   endif

#ifdef USE_MPI
   ints_buf=ints
   call MPI_allreduce(ints_buf,ints,nints,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   ints_buf=maxs
   call MPI_allreduce(ints_buf,maxs,nints,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif
!  storage of some extra quantities:
   maxs(6)=time
   maxs(7)=time_old

   if (maxs(8)>=0 .and. maxs(8)<15 .and. enable_lsf_timelimit==1) then
      write(message,'(a,f20.10)') "LSF timelimit approaching. Stoping at time=",time
      call print_message("** ERROR ****************************************")
      call print_message(message)
      itime_final=itime
   endif
   
   if (maxval(maxs(1:3))> 1000) then
      write(message,'(a,f20.10)') "max U > 1000. Stoping at time=",time
      call print_message("** ERROR ****************************************")
      call print_message(message)
      itime_final=itime
   endif

   if (itime>=itime_final) time_final=time
   call time_control(itime,time,Q,Qhat,q1,q2,q3,work1,work2)

   ! start the main timer *AFTER* 1 complete time step:
   if (itime==0) call wallclock(tmx1)

   itime=itime+1
   if (time >= time_final) exit

enddo


call wallclock(tmx2)
tims(2)=tmx2-tmx1



end subroutine


















