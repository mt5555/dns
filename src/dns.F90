!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Fully resolved compressible Navier Stokes code
! copyright 2001 by Hal Marshal, Mark Taylor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program DNS
use params
use mpi
use structf
implicit none
real*8,save :: Q(nx,ny,nz,n_var)
real*8,save :: q1(nx,ny,nz,n_var)
real*8,save :: work1(nx,ny,nz)
real*8,save :: work2(nx,ny,nz)
character(len=80) message
integer ierr
real*8 tmx1,tmx2,tims_max(ntimers),tims_ave(ntimers)



call init_mpi       
call wallclock(tmx1)  ! wallclock may use MPI timers, call after init_mpi
call init_mpi_comm3d()
call init_grid      


write(message,'(a)') 'Running some tests'
!call print_message(message)
!call test           ! optional testing  routines go here



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Set initial data, make it divergence free
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(message,'(a)') 'Initial data'
call print_message(message)
if (init_cond==0) call init_data_khblob(Q)
if (init_cond==1) call init_data_kh(Q)
if (init_cond==2) call init_data_lwisotropic(Q,q1,work1,work2)

write(message,'(a)') 'Initial data projection'
call print_message(message)
call init_data_projection(Q,work1)  ! impose constrains on initial data
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
call wallclock(tmx1)
call dns_solve(Q,q1,work1,work2)
call wallclock(tmx2)
tims(2)=tmx2-tmx1







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
write(message,'(a,2f9.2,a,f4.3,a)') '   compute_pdf        ',tims_ave(12),tims_max(12)
call print_message(message)
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

call close_mpi


end program DNS









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  main time stepping loop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dns_solve(Q,Qhat,work1,work2)
use params
use mpi

implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

!local variables
real*8,save :: q1(nx,ny,nz,n_var)
real*8,save :: q2(nx,ny,nz,n_var)
real*8,save :: q3(nx,ny,nz,n_var)


real*8  :: time=0
integer :: itime=0,ierr,n
integer :: itime_final=2**30
character(len=80) message
real*8 :: time_old_old,time_old=0,delke_tot,delea_tot,dt
real*8 :: ea_new=0,ea_old
real*8 :: ke_new=0
real*8 :: ints_buf(nints)

ints_timeU=0
ints=0
maxs=0
delt=0

if (time_final<0) then
   ! flag denoting run abs(time_final) time steps, instead of specifing 
   ! a final time
   itime_final=-time_final
   time_final=-time_final
endif



do 
   time_old_old=time_old
   time_old=time
   ea_old=ea_new

   call rk4(time,Q,Qhat,q1,q2,q3,work1,work2)
#ifdef USE_MPI
   ints_buf=ints
   call MPI_allreduce(ints_buf,ints,nints,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   ints_buf=maxs
   call MPI_allreduce(ints_buf,maxs,nints,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif
   ! KE total dissapation 
   dt=time-time_old
   if (dt>0) delke_tot=(ints(1)-ints(6))/dt

   
   ea_new = ints(6) + .5*alpha_value**2 * ints(2) ! computed at time before rk4
   dt=time_old-time_old_old
   if (dt>0) delea_tot=(ea_new-ea_old) /dt



!  storage of some extra quantities:
   maxs(6)=ints_timeU
   maxs(7)=ints_timeDU
   maxs(8)=delea_tot
   maxs(9)=delke_tot
   
   if (maxval(maxs(1:3))> 1000) then
      print *,"max U > 1000. Stoping at time=",time
      itime_final=itime
   endif

   if (itime>=itime_final) time_final=time
   call time_control(itime,time,Q,q1,q2,q3,work1,work2)
   itime=itime+1
   if (time_final>0 .and. time >= time_final) exit

enddo
end subroutine


















