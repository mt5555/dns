!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Fully resolved compressible Navier Stokes code
! copyright 2001 by Hal Marshal, Mark Taylor, Rob Lowrie
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program DNS
implicit none

call init_model
call init_mpi
call init_data
call bc_preloop
call dns_solve

end program DNS









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  main time stepping loop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dns_solve
implicit none
use params
real*8 :: Q(nxd,nyd,nzd,nvard)
real*8 :: time

! set delt based on CFL?
delt = .01

call out(time,bigq)  ! output initial data

do itime=1,itime_max

   call rk4

   if (mod(i,itime_ouput)=1 .or. itime=itime_max .or. error_code>0) then
      call out(time,Q)
   endif

   if (error_code>0) then
      print *,"Error code = ",error_code
      print *,'Stoping at time=",time
      goto 100
   endif
end do
100 continue






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine to take one Runge-Kutta 4th order time step
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk4(time,bigq)
implicit none
use params
real*8 intent(in):: time,Q(nxd,nyd,nzd,nvard)
real*8 :: Q_old(nxd,nyd,nzd,nvard)
real*8 :: Q_tmp(nxd,nyd,nzd,nvard)
real*8 :: rhs(nxd,nyd,nzd,nvard)

biq_old=bigq

! stage 1
call NS(rhs,bigq_old,time)
bigq=bigq+delt*rhs/6.0

! stage 2
bigq_tmp = bigq_old + delt*rhs/2.0
call bc_loop(bigq_tmp)
call NS(rhs,bigq_tmp,time+delt/2.0)
bigq=bigq+delt*rhs/3.0

! stage 3
bigq_tmp = bigq_old + delt*rhs/2.0
call bc_loop(bigq_tmp)
call NS(rhs,bigq_tmp,time+delt/2.0)
bigq=bigq+delt*rhs/3.0

! stage 4
bigq_tmp = bigq_old + delt*rhs
call bc_loop(bigq_tmp)
call NS(rhs,bigq_tmp,time+delt)
bigq=bigq+delt*rhs/6.0

call bc_loop(bigq)
time = time + delt

end subroutine rk4  

















