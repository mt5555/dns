!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Fully resolved compressible Navier Stokes code
! copyright 2001 by Hal Marshal, Mark Taylor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program DNS
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
character*80 message



call init_mpi       

write(message,'(a)') 'Initializing grid and b.c. and reading input file'
call print_message(message)
call init_grid      

write(message,'(a)') 'Setting initial data'
call print_message(message)
call init_data(Q)      ! set up initial data and impose constraints

write(message,'(a)') 'Running some tests'
call print_message(message)
call test           ! optional testing  routines go here

call dns_solve(Q)

write(message,'(a)') 'Cleaning up...'
call close_mpi


end program DNS









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  main time stepping loop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dns_solve(Q)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: time

!local variables
integer itime

call time_control(time,Q)  ! output initial data, choose delt

do 
   call rk4(time,Q)
   call time_control(time,Q)
#ifdef MPI
   call MPI_REDUCE(error_code,MPI_MAX...)
#endif

   if (error_code>0) then
      print *,"Error code = ",error_code
      print *,"Stoping at time=",time
      time_final=time
      call time_control(time,Q)
      exit
   endif
   if (time >= time_final) exit
enddo

end subroutine







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine to take one Runge-Kutta 4th order time step
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk4(time,Q)
use params
implicit none
real*8 :: time,Q(nx,ny,nz,n_var)

! local variables
real*8 :: Q_old(nx,ny,nz,n_var)
real*8 :: Q_tmp(nx,ny,nz,n_var)
real*8 :: rhs(nx,ny,nz,n_var)

Q_old=Q

! stage 1
call ns3D(rhs,Q_old,time)
Q=Q+delt*rhs/6.0


! stage 2
Q_tmp = Q_old + delt*rhs/2.0
call ns3D(rhs,Q_tmp,time+delt/2.0)
Q=Q+delt*rhs/3.0

! stage 3
Q_tmp = Q_old + delt*rhs/2.0
call ns3D(rhs,Q_tmp,time+delt/2.0)
Q=Q+delt*rhs/3.0

! stage 4
Q_tmp = Q_old + delt*rhs
call ns3D(rhs,Q_tmp,time+delt)
Q=Q+delt*rhs/6.0

time = time + delt

end subroutine rk4  

















