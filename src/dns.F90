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
integer ierr



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

write(message,'(a)') 'Initial data projection'
call init_data_projection(Q)  ! impose constrains on initial data


write(message,'(a)') 'Time stepping...'
call print_message(message)
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
integer :: itime=0,ierr
real*8 :: ke_old
real*8 ints(nints),maxs(nints)

ints=0
maxs=0
time=0
itime=0


delt=0
call rk4(time,Q,ints,maxs)
call time_control(itime,time,Q,ints,maxs)  ! output initial data, choose delt

do 
   ke_old=ints(1)
   call rk4(time,Q,ints,maxs)
   itime=itime+1

   ! KE total dissapation 
   ints(2)=(ints(1)-ke_old)/(1d-200+delt)


   if (error_code>0) then
      print *,"Error code = ",error_code
      print *,"Stoping at time=",time
      time_final=time
   endif

   call time_control(itime,time,Q,ints,maxs)
   if (time >= time_final) exit
enddo

end subroutine







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine to take one Runge-Kutta 4th order time step
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk4(time,Q,ints,maxs)
use params
use mpi
implicit none
real*8 :: time,Q(nx,ny,nz,n_var),ints(nints),maxs(nints)

! local variables
real*8 :: Q_old(nx,ny,nz,n_var)
real*8 :: Q_tmp(nx,ny,nz,n_var)
real*8 :: rhs(nx,ny,nz,n_var)
real*8 :: ke_diff
real*8 :: ints_buf(nints)
integer i,j,k,n,ierr


Q_old=Q


! stage 1
call ns3D(rhs,Q_old,time,1,ints)
call divfree(rhs,Q_tmp)
Q=Q+delt*rhs/6.0

! stage 2
Q_tmp = Q_old + delt*rhs/2.0
call ns3D(rhs,Q_tmp,time+delt/2.0,0,ints)
call divfree(rhs,Q_tmp)
Q=Q+delt*rhs/3.0

! stage 3
Q_tmp = Q_old + delt*rhs/2.0
call ns3D(rhs,Q_tmp,time+delt/2.0,0,ints)
call divfree(rhs,Q_tmp)
Q=Q+delt*rhs/3.0

! stage 4
Q_tmp = Q_old + delt*rhs
call ns3D(rhs,Q_tmp,time+delt,0,ints)
Q=Q+delt*rhs/6.0
call divfree(Q,Q_tmp)

ints(1)=0
maxs=0
do n=1,3
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   ints(1)=ints(1)+.5*Q(i,j,k,n)**2
   maxs(n)=max(maxs(n),Q(i,j,k,n))
enddo
enddo
enddo
enddo

time = time + delt

#ifdef USE_MPI
   ints_buf=ints
   call MPI_allreduce(ints_buf,ints,nints,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   ints_buf=maxs
   call MPI_allreduce(ints_buf,maxs,nints,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif


end subroutine rk4  

















