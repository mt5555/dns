#include "macros.h"
subroutine time_control(itime,time,Q)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: time
integer :: itime

! local variables
integer i,j,k,n
character(len=80) message
real*8 remainder, time_target,mumax, umax,time_next,cfl_used_adv,cfl_used_vis,mx
real*8 divx,divi,tmx1,tmx2,del,delke_tot
logical,external :: check_time
logical :: doit

real*8,allocatable,save :: ints_save(:,:),maxs_save(:,:)
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
if (g_nz>1) then
   mumax = mu*(1/delx**2 + 1/dely**2 + 1/delz**2)
else
   mumax = mu*(1/delx**2 + 1/dely**2)
endif

delt = cfl_adv/umax                         ! advective CFL
delt = min(delt,cfl_vis/mumax)              ! viscous CFL
delt = max(delt,delt_min)
delt = min(delt,delt_max)


!
!  restart dumps
!
doit=check_time(itime,time,restart_dt,0,0.0,time_next)
time_target=min(time_target,time_next)
if (doit) then
   call restart_write(time,Q)
endif




!
!  output dumps
!
doit=check_time(itime,time,output_dt,ncustom,custom,time_next)
time_target=min(time_target,time_next)
if (doit) then
   call output_write(time,Q)
endif


!
! accumulate scalers into an array, output during 
! diagnositc output
!
nscalars=nscalars+1
if (nscalars > nsize) then
   nsize=nscalars+100
   ! scalar arrays need to be enlarged
   if (allocated(ints_save)) deallocate(ints_save)
   allocate(ints_save(nints,nsize))
   if (allocated(maxs_save)) deallocate(maxs_save)
   allocate(maxs_save(nints,nsize))
endif
ints_save(1:nints,nscalars)=ints(1:nints)
maxs_save(1:nints,nscalars)=maxs(1:nints)




!
! diagnostic output
!
doit=check_time(itime,time,diag_dt,0,0.0,time_next)
time_target=min(time_target,time_next)
if (doit) then
   call output_diags(time,Q,ints_save,maxs_save,nints,nscalars)
   nscalars=0
endif






doit=check_time(itime,time,screen_dt,0,0.0,time_next)
time_target=min(time_target,time_next)
! also output first 5 timesteps, unless screen_dt==0
if (itime<5 .and. screen_dt/=0) doit=.true.



! compute CFL used for next time step.
cfl_used_adv=umax*delt
cfl_used_vis=mumax*delt


!
! display screen output
!
if (doit) then
   delke_tot=ints(6)

   write(message,'(a,f9.5,a,i5,a,f9.5)') 'time=',time,'(',itime,')  next output=',time_target
   call print_message(message)	

   write(message,'(a,f9.7,a,f6.3,a,f6.3)') 'for next timestep: delt=',delt,' cfl_adv=',cfl_used_adv,' cfl_vis=',cfl_used_vis
   call print_message(message)	

   write(message,'(a,3f22.15)') 'max: (u,v,w) ',maxs(1),maxs(2),maxs(3)
   call print_message(message)	

   call compute_div(Q,divx,divi)
   write(message,'(3(a,e12.5))') 'max(div)=',divx,'   max(vor)',maxs(5)
   call print_message(message)	

   write(message,'(3(a,e12.5))') '<z-vor>=',ints(4),'   <hel>=',ints(5)
   call print_message(message)	

   write(message,'(3(a,e12.5))') 'ke + .5*alpha*<vor**2>',&
      ints(1)+.5*alpha_model*ints(7)
   call print_message(message)	

   write(message,'(a,f13.10,a,f13.4,a,f12.7)') 'ke: ',ints(1),'  enstropy: ',&
        ints(7),'        total d/dt(ke): ',delke_tot
   call print_message(message)	
   write(message,'(a,f12.7,a,f12.7,a,f12.7)') &
     'd/dt(ke) from:  diffusion=',ints(3),' forcing=',ints(2),&
     ' total=',(ints(2)+ints(3))
   call print_message(message)	
   call print_message("")
endif


!
! restrict delt so we hit the next time_target
! (do this here so the courant numbers printed above represent the
! values used if the time step is not lowered for output control)
delt = min(delt,time_target-time)



call wallclock(tmx2)
tims(3)=tmx2-tmx1




end subroutine







subroutine restart_write(time,Q)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: time

! local variables
integer i,j,k,n
real*8 xnx,xny,xnz,xnv
character(len=80) message
character(len=20) tmp
CPOINTER :: fid
integer ierr

n=max(mpidims(1),mpidims(2),mpidims(3))
if (n<10) then
   n=5
else if (n<100) then
   n=4
else if (n<1000) then
   n=3
else if (n<10000) then
   n=2
else 
   call abort("opps, we assumed no more than 10000 cpus along one direction!")
endif
write(tmp,'(i5)') 10000+my_x
message="-" // tmp(n:5)
write(tmp,'(i5)') 10000+my_y
message=message(1:len_trim(message)) // "-" // tmp(n:5)
write(tmp,'(i5)') 10000+my_z
message=message(1:len_trim(message)) // "-" // tmp(n:5) // "-"

write(tmp,'(f10.4)') 10000.0000 + time
message=message(1:len_trim(message)) // tmp(2:10)

message = runname(1:len_trim(runname)) // message(1:len_trim(message)) // ".data"

!open(unit=10,file=message,form='binary')
call copen(message,"w",fid,ierr)
if (ierr/=0) then
   write(message,'(a,i5)') "restart_write(): Error opening file errno=",ierr
   call abort(message)
endif



call cwrite8(fid,time,1)
xnv=n_var
xnx=nslabx
xny=nslaby
xnz=nslabz
call cwrite8(fid,xnx,1)
call cwrite8(fid,xny,1)
call cwrite8(fid,xnz,1)
call cwrite8(fid,xnv,1)
do n=1,n_var
do k=nz1,nz2
do j=ny1,ny2
!do i=nx1,nx2
   call cwrite8(fid,Q(nx1,j,k,n),nx2-nx1+1)
!enddo
enddo
enddo
enddo
call cclose(fid)


end subroutine






subroutine output_diags(time,Q,ints_save,maxs_save,nv,nscalars)
use params
use structf
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: time
integer nv,nscalars
real*8 :: ints_save(nv,nscalars)
real*8 :: maxs_save(nv,nscalars)


! local variables
integer i,j,k,n
integer :: iwave,iwave_max,ierr
real*8 spec_x(0:g_nx/2)
real*8 spec_y(0:g_ny/2)
real*8 spec_z(0:g_nz/2)
real*8 x
real*8,allocatable  ::  spectrum(:),spectrum1(:)
character(len=80) :: message
character :: access
CPOINTER fid

! append to output files, unless time=0 create a new file 
access="a"
if (time==0) access="w"

iwave_max=max(g_nx,g_ny,g_nz)
allocate(spectrum(0:iwave_max))
allocate(spectrum1(0:iwave_max))
spectrum=0
spectrum1=0

do i=1,3
   iwave=iwave_max
   call compute_spectrum(Q(:,:,:,i),spectrum1,spec_x,spec_y,spec_z,iwave,io_pe)
   spectrum=spectrum+.5*spectrum1
enddo
write(message,'(a,f10.4)') " KE spectrum",time
call plotASCII(spectrum,iwave,message(1:25))
!call plotASCII(spec_x,g_nx/2,message)
!call plotASCII(spec_y,g_ny/2,message)
!call plotASCII(spec_z,g_nz/2,message)





if (my_pe==io_pe) then
!   write(message,'(f10.4)') 10000.0000 + time
!   message = runname(1:len_trim(runname)) // message(2:10) // ".spec"
   message = runname(1:len_trim(runname)) // ".spec"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "restart_write(): Error opening file errno=",ierr
      call abort(message)
   endif

   call cwrite8(fid,time,1)
   x=1+iwave; call cwrite8(fid,x,1)
   call cwrite8(fid,spectrum,1+iwave)
   x=1+g_nx/2; call cwrite8(fid,x,1)
   call cwrite8(fid,spec_x,1+g_nx/2)
   x=1+g_ny/2; call cwrite8(fid,x,1)
   call cwrite8(fid,spec_y,1+g_ny/2)
   x=1+g_nz/2; call cwrite8(fid,x,1)
   call cwrite8(fid,spec_z,1+g_nz/2)
   call cclose(fid)
endif



deallocate(spectrum)
deallocate(spectrum1)


!
! output structure functions
!
if (structf_init==1) then
if (my_pe==io_pe) then
!   write(message,'(f10.4)') 10000.0000 + time
!   message = runname(1:len_trim(runname)) // message(2:10) // ".sf"
   message = runname(1:len_trim(runname)) // ".sf"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "outputSF(): Error opening file errno=",ierr
      call abort(message)
   endif
endif
call outputSF(time,fid)
if (my_pe==io_pe) call cclose(fid)
endif



if (my_pe==io_pe) then
!   write(message,'(f10.4)') 10000.0000 + time
   message = runname(1:len_trim(runname)) // ".scalars"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "diag_output(): Error opening .scalars file errno=",ierr
      call abort(message)
   endif
   if (nscalars>0) then
      x=nv; call cwrite8(fid,x,1)
      x=nscalars; call cwrite8(fid,x,1)
      call cwrite8(fid,ints_save,nv*nscalars);
      call cwrite8(fid,maxs_save,nv*nscalars);
   endif
   call cclose(fid)
endif


end subroutine







subroutine output_write(time,Q)
use params
use transpose
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: time

! local variables
integer i,j,k,n
real*8 xnx,xny,xnz
real*8 :: vor(nx,ny,nz,n_var)
real*8 :: d1(nx,ny,nz),work(nx,ny,nz)
character(len=80) message
integer n_var_start,ierr
CPOINTER fid

call vorticity(vor,Q,d1,work)


if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time
   message = runname(1:len_trim(runname)) // message(2:10) // ".vor"
   !open(unit=11,file=message,form='binary')
   call copen(message,"w",fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "restart_write(): Error opening file errno=",ierr
      call abort(message)
   endif


   call cwrite8(fid,time,1)
   xnx=o_nx
   xny=o_ny
   xnz=o_nz
   call cwrite8(fid,xnx,1)
   call cwrite8(fid,xny,1)
   call cwrite8(fid,xnz,1)
   call cwrite8(fid,g_xcord(1),o_nx)
   call cwrite8(fid,g_ycord(1),o_ny)
   call cwrite8(fid,g_zcord(1),o_nz)
endif

call output1(vor(1,1,1,3),work,fid)
if (my_pe==io_pe) call cclose(fid)

end subroutine















logical function check_time(itime,time,dt,ncust,cust,time_next)
use params
implicit none
integer ncust,itime
real*8 :: time,dt,cust(ncust)
real*8 :: time_next  ! output

!local variables
real*8 remainder
integer i
real*8 :: small=1e-7

check_time = .false.
time_next=time_final

if (dt==0) return


if (time>=time_final-small) check_time=.true.
if (time==0) check_time=.true.


if (dt<0) then
   if (mod(itime,nint(abs(dt)))==0) check_time=.true.
endif

if (dt>0) then
   remainder=time - int((small+time)/dt)*dt
   
   if (remainder < small) then
      check_time=.true.
   endif
   time_next = time-remainder+dt

   do i=1,ncust
      if (abs(time-cust(i))<small) then
         check_time=.true.
      else if (time<cust(i)) then
         time_next=min(time_next,cust(i))
         exit 
      endif
   enddo
endif

end function


