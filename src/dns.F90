!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Fully resolved compressible Navier Stokes code
! copyright 2001 by Hal Marshal, Mark Taylor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program DNS
use params
implicit none
character*80 message

write(message,'(a)') 'Initializing model '
call print_message(message)
call init_model     

write(message,'(a)') 'Initializing mpi '
call print_message(message)
call init_mpi       

write(message,'(a)') 'Initializing grid and b.c.'
call print_message(message)
call init_grid      

write(message,'(a)') 'Setting initial data'
call print_message(message)
call init_data      ! set up initial data and impose constraints

write(message,'(a)') 'Running some tests'
call print_message(message)
call test           ! optional testing  routines go here

stop
call dns_solve

end program DNS









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  main time stepping loop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dns_solve
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: time

!local variables
integer itime


! set delt based on CFL?
delt = .01

call out(time,Q)  ! output initial data

do itime=1,itime_max

   call rk4

   if (mod(itime,itime_output)==1 .or. itime==itime_max .or. error_code>0) then
      call out(time,Q)
   endif

   if (error_code>0) then
      print *,"Error code = ",error_code
      print *,"Stoping at time=",time
      exit
   endif
end do
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
call bc_loop(Q_tmp)
call ns3D(rhs,Q_tmp,time+delt/2.0)
Q=Q+delt*rhs/3.0

! stage 3
Q_tmp = Q_old + delt*rhs/2.0
call bc_loop(Q_tmp)
call ns3D(rhs,Q_tmp,time+delt/2.0)
Q=Q+delt*rhs/3.0

! stage 4
Q_tmp = Q_old + delt*rhs
call bc_loop(Q_tmp)
call ns3D(rhs,Q_tmp,time+delt)
Q=Q+delt*rhs/6.0

call bc_loop(Q)
time = time + delt

end subroutine rk4  

















