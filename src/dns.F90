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

write(message,'(a)') 'Running some tests'
call print_message(message)
call test           ! optional testing  routines go here

write(message,'(a)') 'Initial data'
call print_message(message)
call init_data(Q)             ! set up initial data 
call init_data_projection(Q)  ! impose constrains on initial data

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


!local variables
real*8  :: time=0
integer :: itime=0
real*8  :: ints(3)       ! ints(1) = ke
                         ! ints(2) = ke dissapation
                         ! ints(3) = ke dissapation from diffusion

ints=0
time=0
itime=0
call time_control(itime,time,Q,ints)  ! output initial data, choose delt

do 
   call rk4(time,Q,ints)
   itime=itime+1

#ifdef MPI
!   g_code=error_code
!   call MPI_allreduce(g_code,error_code,1,MPI_INTEGER,MPI_MAX,comm_3d,ierr)
#endif
   if (error_code>0) then
      print *,"Error code = ",error_code
      print *,"Stoping at time=",time
      time_final=time
   endif

   call time_control(itime,time,Q,ints)
   if (time >= time_final) exit
enddo

end subroutine







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine to take one Runge-Kutta 4th order time step
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk4(time,Q,ints)
use params
implicit none
real*8 :: time,Q(nx,ny,nz,n_var),ints(3)

! local variables
real*8 :: Q_old(nx,ny,nz,n_var)
real*8 :: Q_tmp(nx,ny,nz,n_var)
real*8 :: rhs(nx,ny,nz,n_var)
real*8 :: ke_old,ke_diff
integer i


ke_old=ints(1)
if (ke_old==0) then
   call compute_ke(Q,io_pe,ke_old)
endif


Q_old=Q

! stage 1
call ns3D(rhs,Q_old,time,ints(3))
call divfree(rhs,Q_tmp(1,1,1,1),Q_tmp(1,1,1,2),Q_tmp(1,1,1,3))
Q=Q+delt*rhs/6.0

! stage 2
Q_tmp = Q_old + delt*rhs/2.0
call ns3D(rhs,Q_tmp,time+delt/2.0,ke_diff)
call divfree(rhs,Q_tmp(1,1,1,1),Q_tmp(1,1,1,2),Q_tmp(1,1,1,3))
Q=Q+delt*rhs/3.0

! stage 3
Q_tmp = Q_old + delt*rhs/2.0
call ns3D(rhs,Q_tmp,time+delt/2.0,ke_diff)
call divfree(rhs,Q_tmp(1,1,1,1),Q_tmp(1,1,1,2),Q_tmp(1,1,1,3))
Q=Q+delt*rhs/3.0

! stage 4
Q_tmp = Q_old + delt*rhs
call ns3D(rhs,Q_tmp,time+delt,ke_diff)
Q=Q+delt*rhs/6.0
call divfree(Q,Q_tmp(1,1,1,1),Q_tmp(1,1,1,2),Q_tmp(1,1,1,3))
call compute_ke(Q,io_pe,ints(1))

! KE total dissapation 
ints(2)=(ints(1)-ke_old)/delt

time = time + delt



end subroutine rk4  

















