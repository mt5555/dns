#include "macros.h"
subroutine time_control(itime,time,Q,Qhat,q1,q2,q3,work1,work2)
use params
use tracers
use transpose
use spectrum
implicit none
real*8 :: Q(nx,ny,nz,n_var)      
real*8 :: Qhat(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: q2(nx,ny,nz,n_var)
real*8 :: q3(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time
integer :: itime

! local variables
integer i,j,k,n
character(len=80) message
character(len=80) fname
real*8 remainder, time_target,mumax, umax,time_next,cfl_used_adv,cfl_used_vis,mx
real*8 tmx1,tmx2,del,lambda,H,ke_diss,epsilon,ens_diss
logical,external :: check_time
logical :: doit_output,doit_diag,doit_restart,doit_screen,doit_model
integer :: n1,n1d,n2,n2d,n3,n3d
real*8,save :: t0,t1=0,ke0,ke1=0,ea0,ea1=0,eta,ett
real*8,save :: ens0,ens1=0
real*8,save :: delea_tot=0,delke_tot=0,delens_tot=0


real*8,allocatable,save :: ints_save(:,:),maxs_save(:,:),ints_copy(:,:)
integer,save :: nscalars=0,nsize=0



call wallclock(tmx1)


time_target=time_final

!
! compute new time step
!
! CFL = delt*umax/delx     delt <= CFL*delx/umax
!
! viscous CFL =  delt*mu/delx^2  delt <= CFL*delx^2/mu
!  
!

umax=maxs(4)
if (ndim==3) then
   mumax = mu/(delx**2) + &
           mu/(dely**2) + &
           mu/(delz**2) 
else
   mumax = mu/(delx**2) + &
           mu/(dely**2) 
endif

if (equations==SHALLOW .and. grav>0) then
   umax=umax+fcor + sqrt(grav*H0)/min(delx,dely)   
endif


delt = cfl_adv/umax                         ! advective CFL
delt = min(delt,cfl_vis/mumax)              ! viscous CFL
delt = max(delt,delt_min)
delt = min(delt,delt_max)

! compute CFL used for next time step.
cfl_used_adv=umax*delt
cfl_used_vis=mumax*delt




!
! accumulate scalers into an array, output during 
! diagnositc output
!
if (nsize==0) then
   nsize=100
   ! scalar arrays need to be allocated
   allocate(ints_save(nints,nsize))
   allocate(maxs_save(nints,nsize))
endif

nscalars=nscalars+1
if (nscalars > nsize) then
   ! scalar arrays need to be enlarged - increase size by 100
   allocate(ints_copy(nints,nsize))  

   ints_copy=ints_save
   deallocate(ints_save)
   allocate(ints_save(nints,nsize+100))
   ints_save(1:nints,1:nsize)=ints_copy(1:nints,1:nsize)

   ints_copy=maxs_save
   deallocate(maxs_save)
   allocate(maxs_save(nints,nsize+100))
   maxs_save(1:nints,1:nsize)=ints_copy(1:nints,1:nsize)

   deallocate(ints_copy)

   nsize=nsize+100
endif
ints_save(1:nints,nscalars)=ints(1:nints)
maxs_save(1:nints,nscalars)=maxs(1:nints)


!
! KE total dissapation 
! t0,t1  three time levels
!
ens0=ens1
ens1=ints(7)   
ke0=ke1
ke1=ints(6)
ea0=ea1
ea1 = ints(6) + .5*alpha_value**2 * ints(2)

t1=max(t1,time_initial)
t0=t1                       ! previous previous time level 
t1=maxs(7)                  ! previous time level  (most quantities known at this time

if (t1-t0>0) then
   delke_tot=(ke1-ke0)/(t1-t0)
   delea_tot=(ea1-ea0)/(t1-t0)
   delens_tot=(ens1-ens0)/(t1-t0)
else
   delke_tot=0
   delea_tot=0
   delens_tot=0
endif



doit_restart=check_time(itime,time,restart_dt,0,0.0,time_next,1)
time_target=min(time_target,time_next)
! dont write a restart file if we just restarted:
if (time==time_initial .and. restart==1) doit_restart=.false.

! ouput of various fields.  
! right now, this is also our restart, so we also output these fields
! if the run is finished.
doit_output=check_time(itime,time,output_dt,ncustom,custom,time_next,1)
time_target=min(time_target,time_next)
! dont write output file if we just restarted:
if (time==time_initial .and. restart==1) doit_output=.false.

!
! diagnostic output (scalars, spectrum).  
!
doit_diag=check_time(itime,time,diag_dt,0,0.0,time_next,1)
time_target=min(time_target,time_next)

!
! model specific output.  
! dns:      structure functions and time averages
! disnvor:  elliptical contours
!
doit_model=check_time(itime,time,diag_dt,0,0.0,time_next,0)
time_target=min(time_target,time_next)

doit_screen=check_time(itime,time,screen_dt,0,0.0,time_next,1)
time_target=min(time_target,time_next)
! also output first 5 timesteps, unless screen_dt==0
if (itime<5 .and. screen_dt/=0) doit_screen=.true.



!
! display screen output
!
if (doit_screen) then
   call print_message("")
   write(message,'(a,f9.5,a,i6,a,f9.5,a,f6.0,i2)') 'time=',time,'(',itime,')  next output=',time_target, &
      '  LSF minutes left: ',maxs(8),enable_lsf_timelimit
   call print_message(message)	

   write(message,'(a,f9.7,a,f6.3,a,f6.3)') 'for next timestep: delt=',delt,' cfl_adv=',cfl_used_adv,' cfl_vis=',cfl_used_vis
   call print_message(message)	

   write(message,'(a,3f21.14)') 'max: (u,v,w) ',maxs(1),maxs(2),maxs(3)
   call print_message(message)	

   write(message,'(3(a,e12.5))') '<z-vor>=',ints(4),'   <hel>=',ints(5),&
           '   max(vor)',maxs(5)
   call print_message(message)	

   ke_diss = ints(10)  
   ! 
   ! using lambda**2 = <u1 u1>/<u1,1 u1,1>
   ! and   <u1,1 u1,1> = (1/15) || grad(u) ||^2
   !       < u    u  > = (1/3)  || u ||^2
   !
   if (mu>0 .and. ndim>2 .and. ke_diss<0) then
      lambda=sqrt(  5*(2*ints(6))/ints(2)  )
      epsilon=-ke_diss
      write(message,'(3(a,f12.5))') 'R_lambda=',lambda*sqrt(2*ints(6))/mu, &
           '  R=',1/mu
      call print_message(message)	
      
      
      ! K. microscale
      ! eta = (mu^3/epsilon)^.25
      ! epsilon = delke_tot
      eta = (mu**3 / epsilon)**.25
      write(message,'(a,3f13.4)') 'mesh spacing(eta): ',&
           delx/eta,dely/eta,delz/eta
      call print_message(message)	
      
      !eddy turn over time
      ! 
      ett=2*ke1/epsilon
      write(message,'(a,3f13.4)') 'eddy turnover time: ',ett
      call print_message(message)	
   endif
   
   !
   ! alpha model information
   !
   if (alpha_value>0) then
   write(message,'(a,f13.10,a,f13.4,a,f12.7)') 'Ea: ',&
        ea1,'   |gradU|: ',ints(2),'         total d/dt(Ea):',delea_tot
   call print_message(message)	

   write(message,'(3(a,f12.7))') 'd/dt(Ea) vis=',&
        ke_diss-mu*alpha_value**2*ints(1),&
        ' f=',ints(3) - alpha_value**2*ints(9),&
         '                      tot=',&
        ke_diss-mu*alpha_value**2*ints(1) + ints(3) - alpha_value**2*ints(9)
   call print_message(message)	
   endif

   if (ke1>99) then
   write(message,'(a,f13.8,a,f13.4,a,f12.7)') 'ke: ',ke1,'  enstropy: ',&
        ints(7),'         total d/dt(ke):',delke_tot
   else
   write(message,'(a,f13.10,a,f13.4,a,f12.7)') 'ke: ',ke1,'  enstropy: ',&
        ints(7),'         total d/dt(ke):',delke_tot
   endif
   call print_message(message)	

   if (equations==NS_PSIVOR .or. (equations==SHALLOW .and. alpha_value>0)) then
      ! ke dissapation not computed correctly in the above cases.
   else
      write(message,'(a,f12.7,a,f12.7,a,f12.7,a,f12.7)') &
           'd/dt(ke) vis=',ke_diss,' f=',ints(3),' alpha=',ints(8),&
           ' total=',(ke_diss+ints(3)+ints(8))
      call print_message(message)	
   endif
   if (equations==NS_PSIVOR) then
      !multiply by 2 because enstrophy is vor**2, not .5*vor**2
      ens_diss = 2*ints(5)  
      write(message,'(a,f15.7,a,f15.7)') &
           'total d/dt(w) vis=',delens_tot,'   mu <w,w_xx>= ',ens_diss
      call print_message(message)	
   endif
endif

!
!  restart dumps
!
if (doit_restart) then
   call print_message("writing restart file...")
   call multfile_io(time,Q)
endif





!
!  output dumps
!
if (doit_output) then
   write(message,'(a,f10.4)') "writing output files at t=",time
   call print_message(message)
   call tracers_save(io_pe,time)

   write(message,'(f10.4)') 10000.0000 + time

   if (equations==NS_UVW .and. rw_spec) then
      call transpose_from_z_3d(Qhat,q1)

      ! NS, primitive variables
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".us"
      call singlefile_io2(time,q1(1,1,1,1),fname,work1,work2,0,io_pe,.true.)

      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".vs"
      call singlefile_io2(time,q1(1,1,1,2),fname,work1,work2,0,io_pe,.true.)
      if (n_var==3) then
         fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".ws"
         call singlefile_io2(time,q1(1,1,1,n_var),fname,work1,work2,0,io_pe,.true.)
      endif

   else if (equations==NS_UVW) then

      ! NS, primitive variables
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".u"
      call singlefile_io(time,Q(1,1,1,1),fname,work1,work2,0,io_pe)
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".v"
      call singlefile_io(time,Q(1,1,1,2),fname,work1,work2,0,io_pe)
      if (n_var==3) then
         fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".w"
         call singlefile_io(time,Q(1,1,1,n_var),fname,work1,work2,0,io_pe)
      endif
      if (ndim==2) then
         call vorticity2d(q1,Q,work1,work2)
         fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".vor"
         call singlefile_io(time,q1,fname,work1,work2,0,io_pe)
      endif

   else if (equations==SHALLOW) then
      ! shallow water 2D
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".u"
      call singlefile_io(time,Q(1,1,1,1),fname,work1,work2,0,io_pe)
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".v"
      call singlefile_io(time,Q(1,1,1,2),fname,work1,work2,0,io_pe)
      if (n_var==3) then
         fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".h"
         call singlefile_io(time,Q(1,1,1,n_var),fname,work1,work2,0,io_pe)
      endif
      if (ndim==2) then
         call vorticity2d(q1,Q,work1,work2)
         fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".vor"
         call singlefile_io(time,q1,fname,work1,work2,0,io_pe)
      endif
   else if (equations==NS_PSIVOR) then
      ! 2D NS psi-vor formulation
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".vor"
      call singlefile_io(time,Qhat(1,1,1,1),fname,work1,work2,0,io_pe)
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".psi"
      call singlefile_io(time,Qhat(1,1,1,2),fname,work1,work2,0,io_pe)
   endif
endif


!
! diagnostic output
!
if (compute_transfer) then
   compute_transfer=.false.
   call compute_tran(time,Q,q1,work1,work2)
   call output_tran(time,Q,q1,q2,q3,work1,work2)
endif
if (doit_diag) then
   if ( g_bdy_x1==PERIODIC .and. &
        g_bdy_y1==PERIODIC .and. &
        g_bdy_z1==PERIODIC) then
      call compute_spec(time,Q,q1,work1,work2)
      call output_spec(time,Q,q1,q2,q3,work1,work2)

!     set this flag so that for next timestep, we will compute and save
!     spectral transfer functions:
      compute_transfer=.true.
   endif

   call output_scalars(time,ints_save,maxs_save,nints,nscalars)
   nscalars=0
else if (diag_dt==0) then
   ! if diagnostics are turned off, dont save the scalars!
   nscalars=0  
endif


!
! model specific output
!
if (doit_model) then
   ! model specific output:
   call output_model(time,Q,Qhat,q1,q2,q3,work1,work2)
endif





!
! restrict delt so we hit the next time_target
! (do this here so the courant numbers printed above represent the
! values used if the time step is not lowered for output control)
delt = min(delt,time_target-time)

call wallclock(tmx2)
tims(3)=tims(3) + (tmx2-tmx1)

end subroutine






logical function check_time(itime,time,dt,ncust,cust,time_next,include_final)
!
!  determine if current time is an "output" time, as set by:
!     dt:   output every 'dt' time
!   cust:   list of custome output times
!  
!  include_final=1   also output if time=time_final, even if time_final
!                    is not a multiple of dt or listed in 'cust'
!
use params
implicit none
integer ncust,itime
real*8 :: time,dt,cust(ncust)
real*8 :: time_next  ! output
integer :: include_final

!local variables
real*8 remainder
integer i
real*8 :: small=1e-7

check_time = .false.
time_next=time_final

! custom times always take precedence:
do i=1,ncust
   if (abs(time-cust(i))<small) then
      check_time=.true.
   else if (time<cust(i)) then
      time_next=min(time_next,cust(i))
      exit 
   endif
enddo

if (dt==0) return


if (include_final==1 .and. time>=time_final-small) check_time=.true.
if (time==0) check_time=.true.


if (dt<0) then
   if (mod(itime,nint(abs(dt)))==0) check_time=.true.
endif

if (dt>0) then
   remainder=time - int((small+time)/dt)*dt
   
   if (remainder < small) then
      check_time=.true.
   endif
   time_next = min(time_next,time-remainder+dt)

endif

end function
