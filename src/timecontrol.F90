#include "macros.h"
subroutine time_control(itime,time,Q,Qhat,q1,q2,q3,work1,work2)
use params
use transpose
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
real*8 remainder, time_target,psmax,mumax, umax,time_next,cfl_used_adv,cfl_used_vis,mx
real*8 :: cfl_used_psvis
real*8 tmx1,tmx2,del,lambda,H,ke_diss,epsilon,ens_diss,xtmp
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
! linear advection:
!      u_t + umax grad(u) + f cross u = 0
! FFT:   uk_t + i (2*pi*k*umax + f) uk = 0
! stabilityis based on:  2*pi*k*umax + f  =   pi*umax/delx + f
!                                         =   pi*(umax/delx + f/pi)
!  
! So if we base our CFL on umax/delx, we need to divide f by pi.
!
! CFL = delt*umax/delx       delt <= CFL*delx/umax
!
! with fcor:
!   CFL = delt*(umax/delx + f/pi)   delt <= CFL/(umax/delx + f/pi)
!
! viscous CFL =  delt*mu/delx^2  delt <= CFL*delx^2/mu
!  
!

umax=maxs(4)+fcor/pi

if (ndim==3) then
   mumax = mu/(delx**2) + &
           mu/(dely**2) + &
           mu/(delz**2) 
else
   mumax = mu/(delx**2) + &
           mu/(dely**2) 
endif
psmax=0
if (npassive>0) then
   do n=np1,np2
      if (schmidt(n)>0) psmax=max(mumax/schmidt(n),psmax)	
   enddo
endif


if (equations==SHALLOW .and. grav>0) then
   umax=umax + sqrt(grav*H0)/min(delx,dely)   
endif

!if u=0, take very small timesteps 
xtmp=1e-8
umax=max(umax,xtmp)
delt = cfl_adv/umax                         ! advective CFL
if (mu>0) delt = min(delt,cfl_vis/mumax)    ! viscous CFL
if (psmax>0) delt = min(delt,cfl_vis/psmax)    ! viscous CFL for passive scalars
delt = max(delt,delt_min)
delt = min(delt,delt_max)

! compute CFL used for next time step.
cfl_used_adv=umax*delt
cfl_used_vis=mumax*delt
cfl_used_psvis=psmax*delt



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
if (infinite_alpha==0) then
   ea1 = ints(6) + .5*alpha_value**2 * ints(2)
else
   ea1 = .5*ints(2)
endif

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



doit_restart=check_time(itime,time,restart_dt,0,0.0,time_next,1,0)
time_target=min(time_target,time_next)

! ouput of various fields.  
! right now, this is also our restart, so we also output these fields
! if the run is finished.
doit_output=check_time(itime,time,output_dt,ncustom,custom,time_next,1,0)
time_target=min(time_target,time_next)

!
! diagnostic output (output to .scalars file)
! every diag_dt, and final time
doit_diag=check_time(itime,time,diag_dt,0,0.0,time_next,1,0)
time_target=min(time_target,time_next)

!
! model specific output.  usually expensive, so dont
!                         compute on first or last timestep
!                         unless it happens to be a daig_dt interval
! dns:      spectrum, structure functions and time averages
! dnsvor:   tracers, contour info
!
doit_model=check_time(itime,time,model_dt,0,0.0,time_next,0,0)
time_target=min(time_target,time_next)

doit_screen=check_time(itime,time,screen_dt,0,0.0,time_next,1,0)
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

   write(message,'(a,f9.7,a,f6.3,a,f6.3,a,f6.3)') &
           'next delt=',delt,' cfl_adv=',cfl_used_adv,&
           ' cfl_vis=',cfl_used_vis,' cfl_passive_vis=',cfl_used_psvis
   call print_message(message)	

   write(message,'(a,3f21.14)') 'max: (u,v,w) ',maxs(1),maxs(2),maxs(3)
   call print_message(message)	

   if (npassive>0) then
   write(message,'(a,2f21.14)') 'min/max first passive scalar: ',-maxs(11),maxs(10)
   call print_message(message)	
   endif

   write(message,'(3(a,e12.5))') '<z-vor>=',ints(4),'   <hel>=',ints(5),&
           '   max(vor)',maxs(5)
   call print_message(message)	

   ke_diss = ints(10)  
   ! 
   ! epsilon = mu* || grad(u) ||^2
   ! using lambda**2 = <u1 u1>/<u1,1 u1,1>
   ! and   <u1,1 u1,1> = (1/15) || grad(u) ||^2   =   (1/15)*epsilon/mu
   !       < u1,u1  > =  (1/3)  || u ||^2         = (2/3) KE
   !
   ! lambda**2 = (2/3) KE / (1/15)epsilon/mu  = 10 KE mu / epsilon
   ! R_lambda = lambda*sqrt(2 KE/3)/mu  = sqrt(10 KE mu/epsilon) * sqrt(2 KE/3)/mu
   !                                    = sqrt(20/3 KE mu /epsilon)/mu
   !                                      sqrt(20/3 KE /epsilon) * sqrt(mu)
   !    
   ! units:  epsilon:  m^2/s^3
   !         KE        m^2/s^2
   !         mu        m^2/s
   !
   if (mu>0 .and. ndim>2 .and. ke_diss<0) then
!      lambda=sqrt(  5*(2*ints(6))/ints(2)  )

!      epsilon=-ke_diss
      if (infinite_alpha==0) then
         epsilon=-(  ke_diss-mu*alpha_value**2*ints(1) )
      else
         epsilon=-(  -mu*ints(1)  )
      endif
      lambda=sqrt( mu*(2*ea1/3) / (epsilon/15) )


      write(message,'(3(a,f12.5))') 'R_lambda=',lambda*sqrt(2*ea1/3)/mu, &
           '  R=',1/mu,' lambda=',lambda
      call print_message(message)	
      
      
      ! K. microscale
      ! eta = (mu^3/epsilon)^.25
      ! epsilon = delke_tot
      eta = (mu**3 / epsilon)**.25
      write(message,'(a,3f13.4)') 'mesh spacing/eta: ',&
           delx/eta,dely/eta,delz/eta
      call print_message(message)	
      if (dealias==1)  then
         write(message,'(a,3f13.4)') '2/3 dealiasing kmax*eta: ',&
           2*pi*g_nmin*eta/3
         call print_message(message)	
      endif
      if (dealias==2)  then
         write(message,'(a,3f13.4)') 'spherical dealiasing kmax*eta: ',&
           2*pi*sqrt(2.)*g_nmin*eta/3
         call print_message(message)	
      endif
      
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
      write(message,'(a,f13.8,a,f13.4,a,f12.7)') 'Ea: ',&
           ea1,'   |gradU|: ',ints(2),'         total d/dt(Ea):',delea_tot
      call print_message(message)	
      
      if (infinite_alpha==0) then
         write(message,'(3(a,f12.7))') 'd/dt(Ea) vis=',&
              ke_diss-mu*alpha_value**2*ints(1),&
              ' f=',ints(3) - alpha_value**2*ints(9),&
              '                      tot=',&
              ke_diss-mu*alpha_value**2*ints(1) + ints(3) - alpha_value**2*ints(9)
      else
         write(message,'(3(a,f12.7))') 'd/dt(Ea) vis=',&
              -mu*ints(1),' f=',-ints(9),&
              '                      tot=', -mu*ints(1) - ints(9)
      endif
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

   if (equations==NS_PSIVOR ) then
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
   call multfile_io(time,Q(1,1,1,1),0)  ! write headers
   do i=1,3	
      call multfile_io(time,Q(1,1,1,i),i) ! write u,v,w files
   enddo
endif





!
!  output dumps
!

if (doit_output) then
   write(message,'(a,f10.4)') "writing output files at t=",time
   call print_message(message)
   if (equations==NS_UVW .and. w_spec) then
      call transpose_from_z_3d(Qhat,q1)
      ! convert to nx,ny,nz dimensions for output:
      call output_uvw(runname,time,q1,q2,work1,work2,1)
   else if (equations==NS_PSIVOR) then
      ! 2D NS psi-vor formulation
      call output_uvw(runname,time,Qhat,q1,work1,work2,1)
   else
      call output_uvw(runname,time,Q,q1,work1,work2,1)
   endif
   call output_passive(runname,time,Q,work1,work2)
   call print_message("done with output")

endif


!
! diagnostic output
!
if (doit_diag) then
   call output_scalars(time,ints_save,maxs_save,nints,nscalars)
   nscalars=0
else if (diag_dt==0) then
   ! if diagnostics are turned off, dont save the scalars!
   nscalars=0  
endif


!
! model specific output
!
call output_model(doit_model,doit_diag,time,Q,Qhat,q1,q2,q3,work1,work2)





!
! restrict delt so we hit the next time_target
! (do this here so the courant numbers printed above represent the
! values used if the time step is not lowered for output control)
delt = min(delt,time_target-time)

call wallclock(tmx2)
tims(3)=tims(3) + (tmx2-tmx1)

end subroutine






logical function check_time(itime,time,dt,ncust,cust,time_next,include_final,&
   include_initial)
!
!  determine if current time is an "output" time, as set by:
!     dt:   output every 'dt' time
!   cust:   list of custome output times
!  
!  include_final=1   also output if time=time_final, even if time_final
!                    is not a multiple of dt or listed in 'cust'
!
!  include_initial=1   for restart runs, we skip the intial time
!                      unless this flag is set to 1.  
!
use params
implicit none
integer ncust,itime
real*8 :: time,dt,cust(ncust)
real*8 :: time_next  ! output
integer :: include_final,include_initial

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

! skip initial time on restart runs:
if (include_initial==0) then
   if (time==time_initial .and. restart==1) return
endif


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
