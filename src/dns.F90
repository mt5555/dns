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
call init_data(Q)             ! set up initial data 

write(message,'(a)') 'Initial data projection'
call init_data_projection(Q)  ! impose constrains on initial data

call wallclock(tmx2)
tims(1)=tmx2-tmx1

write(message,'(a)') 'Time stepping...'
call print_message(message)

call wallclock(tmx1)
call dns_solve(Q)
call wallclock(tmx2)
tims(2)=tmx2-tmx1

#ifdef USE_MPI
   temp=tims
   call MPI_allreduce(temp,tims,ntimers,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif


call print_message('CPU times:')
write(message,'(a,f10.5,a)') 'initialization: ',tims(1),' s'
call print_message(message)
write(message,'(a,f10.5,a)') 'dns_solve:      ',tims(2),' s'
call print_message(message)
write(message,'(a,f10.5,a,f4.3,a)') '   time_control:   ',tims(3),' s (',tims(3)/tims(2),')'
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
!call timestep(time,Q,ints,maxs)
!call time_control(itime,time,Q,ints,maxs)


do 
   ke_old=ints(1)
   call timestep(time,Q,ints,maxs)

   ! KE total dissapation 
   ints(2)=(ints(1)-ke_old)/(1d-200+delt)


   if (error_code>0) then
      print *,"Error code = ",error_code
      print *,"Stoping at time=",time
      time_final=time
   endif

   call time_control(itime,time,Q,ints,maxs)
   itime=itime+1
   if (time >= time_final) exit
enddo

end subroutine







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  time step
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine timestep(time,Q,ints,maxs)
use params
use mpi
implicit none
real*8 :: time,Q(nx,ny,nz,n_var),ints(nints),maxs(nints)

! local variables
real*8 :: ke_diff
real*8 :: ints_buf(nints)
integer i,j,k,n,ierr

call rk4(time,Q,ints,maxs)

   ! compute KE, max U, KE disspation
   ints(1)=0
   maxs(1:3)=0
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

#ifdef USE_MPI
   ints_buf=ints
   call MPI_allreduce(ints_buf,ints,nints,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   ints_buf=maxs
   call MPI_allreduce(ints_buf,maxs,nints,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif


end subroutine 








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine to take one Runge-Kutta 4th order time step
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk4(time,Q,ints,maxs)
use params
implicit none
real*8 :: time,Q(nx,ny,nz,n_var),ints(nints),maxs(nints)

! local variables
real*8 :: Q_old(nx,ny,nz,n_var)
real*8 :: Q_tmp(nx,ny,nz,n_var)
real*8 :: rhs(nx,ny,nz,n_var)
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


time = time + delt

end subroutine rk4  



#if 0

U = Un
G = F(U,tn)
U = U + 1/3 dt G
G = -5/9 G + F(U,tn+1/3 dt)
U = U + 15/16 dt G
G = - 153/128 G + F(U,tn+3/4dt)
U(n+1)=U + 8/15 G

#endif














