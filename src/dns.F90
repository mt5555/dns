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
use fft_interface
implicit none
#ifdef MALLOC_Q
real*8,allocatable :: Q(:,:,:,:)        ! grid space state variable
real*8,allocatable :: Qhat(:,:,:,:)     ! Fourier space state variable
real*8,allocatable :: q1(:,:,:,:)       ! work arrays
real*8,allocatable :: q2(:,:,:,:)       ! work arrays
real*8,allocatable :: rhs(:,:,:,:)      ! right hand side of equations
real*8,allocatable :: work1(:,:,:)
real*8,allocatable :: work2(:,:,:)
#else
real*8,save :: Q(nx,ny,nz,n_var)        ! grid space state variable
real*8,save :: Qhat(nx,ny,nz,n_var)     ! Fourier space state variable
real*8,save :: q1(nx,ny,nz,n_var)       ! work arrays
real*8,save :: q2(nx,ny,nz,n_var)       ! work arrays
real*8,save :: rhs(nx,ny,nz,n_var)      ! right hand side of equations
real*8,save :: work1(nx,ny,nz)
real*8,save :: work2(nx,ny,nz)
#endif
character(len=80) message
integer :: ierr,i,j,k,n,itime=0
real*8 :: tmx1,tmx2,tims_max(ntimers),tims_ave(ntimers)
real*8 :: nf,flop,flops
logical :: power3


! initialize global constants:
ints=0
maxs=0

#ifdef MALLOC_Q
allocate(Q(nx,ny,nz,n_var))
allocate(Qhat(nx,ny,nz,n_var))
allocate(q1(nx,ny,nz,n_var))
allocate(q2(nx,ny,nz,n_var))
allocate(rhs(nx,ny,nz,n_var))
allocate(work1(nx,ny,nz))
allocate(work2(nx,ny,nz))
#endif


call set_sighandler()
call init_mpi       
call wallclock(tmx1)  ! wallclock may use MPI timers, call after init_mpi
call init_mpi_comm3d()
!if(do_mpi_io .and. io_pe==my_pe .and. g_nz>2000) &
!       call system("touch /users/taylorm/RUNNING")



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call init_model()
! fftw needs to do FFTs to initialize
call fft_interface_init(work1,work2,nx,ny,nz)  



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! optional testing  routines go here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!call iso_stats(Q,Qhat,work1,work2)
call test   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Set initial data, make it divergence free
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(message,'(a)') 'Initial data'
call print_message(message)
! set the random seed (used for some initial conditions) 
! otherwise it will be the same for all CPUs, producing a bad initial condition.  
call f90random_seed(my_pe)  
Q=0
call init_data(Q,Qhat,q1,work1,work2)
! reset the random seed (used for some forcings) based on initial_time, 
! so we dont have the same seed every time we restart. 
call f90random_seed(nint(my_pe+1000*time_initial))


#ifdef USE_MPI
call mpi_barrier(comm_3d,ierr)
#endif


call wallclock(tmx2)
tims(1)=tmx2-tmx1  ! total initialization


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Solve the equations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(message,'(a)') 'Time stepping...(calling dns_solve)'
call print_message(message)
call dns_solve(Q,Qhat,q1,q2,rhs,work1,work2,itime)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Print timing information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
tims_max=tims
tims_ave=tims
#ifdef USE_MPI
   call mpi_allreduce(tims,tims_max,ntimers,MPI_REAL8,MPI_MAX,comm_3d,ierr)
   call mpi_allreduce(tims,tims_ave,ntimers,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   tims_ave=tims_ave/ncpus
#endif
tims_max=tims_max/60
tims_ave=tims_ave/60

! sum the total of all transposes
! transpose_tot=sum(tims(6:11))




call print_message('CPU times (min):   (avg/max)')
write(message,'(a,2f9.2,a)') 'initialization: ',tims_ave(1),tims_max(1)
call print_message(message)
write(message,'(a,2f9.2,a)') 'initial output: ',tims_ave(4),tims_max(4)
call print_message(message)

if (itime==0) itime=1
write(message,'(a,2f9.2,a,2f10.5)') 'dns_solve:      ',tims_ave(2),tims_max(2),&
'  per timestep: ',tims_ave(2)/itime,tims_max(2)/itime
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
   write(message,'(a,2f9.4,a,f4.3,a)') '   transpose_to_z     ',tims_ave(6),tims_max(6)
   call print_message(message)
   write(message,'(a,2f9.4,a,f4.3,a)') '   transpose_from_z   ',tims_ave(7),tims_max(7)
   call print_message(message)
   write(message,'(a,2f9.4,a,f4.3,a)') '   transpose_to_x     ',tims_ave(8),tims_max(8)
   call print_message(message)
   write(message,'(a,2f9.4,a,f4.3,a)') '   transpose_from_x   ',tims_ave(9),tims_max(9)
   call print_message(message)
   write(message,'(a,2f9.4,a,f4.3,a)') '   transpose_to_y     ',tims_ave(10),tims_max(10)
   call print_message(message)
   write(message,'(a,2f9.4,a,f4.3,a)') '   transpose_from_y   ',tims_ave(11),tims_max(11)
   call print_message(message)
   write(message,'(a,2f9.4)')          '   traspose total     ',sum(tims_ave(6:11))
   call print_message(message)
endif
if (maxval(tims_max(18:19))>0) then
   write(message,'(a,2f9.4,a,f5.1,a)') '   FFT                ',tims_ave(18),tims_max(18),' count per timestep=',&
   ncalls(18)/real(itime),' N^2'
   call print_message(message)
   write(message,'(a,2f9.4,a,f5.1,a)') '   iFFT               ',tims_ave(19),tims_max(19),' count per timestep=',&
   ncalls(19)/real(itime),' N^2'
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
if (tims_max(16)>0) then
   write(message,'(a,2f9.2,a,f4.3,a)') '   tracer_advance     ',tims_ave(16),tims_max(16)
   call print_message(message)
endif
if (tims_max(17)>0) then
   write(message,'(a,2f9.2,a,f4.3,a)') '   elliptical contours',tims_ave(17),tims_max(17)
   call print_message(message)
endif

! if grid is a cube and only powers of 2 and at most one power of 3
! we can compute the flops:
! formulas accurate to < 1% on resolutions 16 ... 256
if (forcing_type==1 .and. equations==NS_UVW .and. dealias==2 .and. &
numerical_method==FOURIER) then 
if (g_nx==g_ny .and. g_nx == g_nz) then
   nf=g_nx
   power3=.false.
   if (mod(g_nx,3)==0) then
      nf=nf/3
      power3=.true.
   endif
   ! nf should only contain powers of 2:
   ! flop = floating point operations for 1 timestep
   if ( ( nf - 2**nint(log(nf)/log(2d0)) ) == 0 ) then
      if (power3) then
         ! 4(  340 N^3 + 2.24 27 N^3 log2(N) )
         nf=g_nx
         flop = 4 *(340*nf**3 + 2.24*27*(nf**3)*log(nf)/log(2d0) )
      else
         ! 4(  297 N^3 + 2.28 27 N^3 log2(N) )
         nf=g_nx
         flop = 4 *(297*nf**3 + 2.28*27*(nf**3)*log(nf)/log(2d0) )
      endif
      tmx1 = 60*tims_max(2)/itime  ! time in seconds
      flops = flop/tmx1/ncpu_x/ncpu_y/ncpu_z
      write(message,'(a,e12.4)')  'Estimated FLOP per timestep: ',flop
      call print_message(message)
      write(message,'(a,f10.2)')  'MFLOPS per cpu: ',flops/1e6
      call print_message(message)
   endif
endif
endif



call close_mpi


end program DNS









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  main time stepping loop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dns_solve(Q,Qhat,q1,q2,rhs,work1,work2,itime)
use params
use mpi
use tracers

implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: q2(nx,ny,nz,n_var)       ! work arrays
real*8 :: rhs(nx,ny,nz,n_var)      ! right hand side of equations
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
integer :: itime

!local variables

real*8  :: time=0
integer :: ierr,n
integer :: itime_final=2**30
character(len=80) message
real*8 :: time_old=0
real*8 :: ea_new=0,ea_old
real*8 :: ke_new=0
real*8 :: ints_buf(nints)
real*8 :: tmx1,tmx2
integer :: i,time_needed


delt=0
! quit if we are within 'time_needed' min of being killed:
time_needed = 15  
! for 2000^3 on 512 cpus, we need 10min per timestep:
! if (g_nx>=2048 .and. g_ny>=2048 .and. g_nz>=2048) time_needed=15

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

!   if(do_mpi_io .and. io_pe==my_pe .and. g_nz>2000) &
!       call system("touch /users/taylorm/RUNNING")
   call rk4(time,Q,Qhat,q1,q2,rhs,work1,work2)

   call caught_sig(i); maxs(9)=i;

#ifdef USE_MPI
   ints_buf=ints
   call mpi_allreduce(ints_buf,ints,nints,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   ints_buf=maxs
   call mpi_allreduce(ints_buf,maxs,nints,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif
   g_u2xave=ints(2)

!  storage of some extra quantities:
   maxs(6)=time
   maxs(7)=time_old

   if (maxs(9)>0) then
      write(message,'(a,f10.0)') "Stopping because of caught signal = ",maxs(9)
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
   call time_control(itime,time,Q,Qhat,q1,q2,rhs,work1,work2)

   if (itime==0) then
      ! set all timers to zero so they dont include initialization
      ! except for tims(1) and tims(4) which time initializaiton stuff.  
      tims(4) = tims(3)   ! tims(3) = total time in time_control(), tims(4)=time of first call to time_control()
      tims(2:3)=0
      tims(5:ntimers)=0
      call wallclock(tmx1)     ! start the main timer *AFTER* 1 complete time step:
   endif

   itime=itime+1
   if (time >= time_final-1e-15) exit

enddo
itime=itime -1  ! total number of time steps taken
                ! not couting the first time step with delt=0
                ! (that timestep is not included in cpu totals either)

call wallclock(tmx2)
tims(2)=tmx2-tmx1



end subroutine
















