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
real*8 :: Q(nx,ny,nz,n_var)
character*80 message
integer ierr
real*8 tmx1,tmx2,temp(ntimers)


call wallclock(tmx1)
call init_mpi       
call init_mpi_comm3d()
call init_grid      


write(message,'(a)') 'Running some tests'
call print_message(message)
call test           ! optional testing  routines go here

write(message,'(a)') 'Initial data'
call print_message(message)
!call init_data_test(Q)             ! set up initial data 
!call init_data_khblob(Q)
call init_data_kh(Q)             ! set up initial data 

write(message,'(a)') 'Initial data projection'
call init_data_projection(Q)  ! impose constrains on initial data

call wallclock(tmx2)

write(message,'(a)') 'Time stepping...'
call print_message(message)

! set all counters to zero (so they dont include initialization)
tims=0
! tims(1) times the total initialization
tims(1)=tmx2-tmx1

call wallclock(tmx1)
call dns_solve(Q)
call wallclock(tmx2)
tims(2)=tmx2-tmx1

#ifdef USE_MPI
   temp=tims
   call MPI_allreduce(temp,tims,ntimers,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif

tims=tims/60
call print_message('CPU times:')
write(message,'(a,f9.2,a)') 'initialization: ',tims(1),' m'
call print_message(message)
write(message,'(a,f9.2,a)') 'dns_solve:      ',tims(2),' m'
call print_message(message)
write(message,'(a,f9.2,a,f4.3,a)') '   time_control + children:    ',tims(3),' m'
call print_message(message)
write(message,'(a,f9.2,a,f4.3,a)') '   RHS + children:             ',tims(5),' m'
call print_message(message)
write(message,'(a,f9.2,a,f4.3,a)') '   transpose + children:       ',tims(4),' m'
call print_message(message)

call close_mpi


end program DNS









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  main time stepping loop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dns_solve(Q)
use params
use mpi

implicit none
real*8 :: Q(nx,ny,nz,n_var)

!local variables
real*8  :: time=0
integer :: itime=0,ierr,n
character*80 message
real*8 :: ke_old,time_old
real*8 :: ints_buf(nints)

ints=0
maxs=0
time=0
itime=0
delt=0


do 
   time_old=ints_timeU
   ke_old=ints(1)

   call rk4(time,Q)

#ifdef USE_MPI
   ints_buf=ints
   call MPI_allreduce(ints_buf,ints,nints,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   ints_buf=maxs
   call MPI_allreduce(ints_buf,maxs,nints,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif
   ! KE total dissapation 
   delke_tot=ints_timeU-time_old
   if (delke_tot>0) delke_tot=(ints(1)-ke_old)/delke_tot

   
   if (maxval(maxs(1:3))> 1000) then
      print *,"max U > 1000. Stoping at time=",time
      time_final=time
   endif

   call time_control(itime,time,Q,ints,maxs)
   itime=itime+1
   if (time >= time_final) exit
enddo

end subroutine























