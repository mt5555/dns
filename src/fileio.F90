#include "macros.h"
subroutine time_control(itime,time,Q)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: time
integer :: itime

! local variables
integer i,j,k,n
character*80 message
real*8 remainder, time_target,mumax, umax,time_next,cfl_used_adv,cfl_used_vis,mx
real*8 divx,divi,tmx1,tmx2,del
logical,external :: check_time
logical :: doit

call wallclock(tmx1)
time_target = time_final

!
! compute new time step
!
! CFL = delt*umax/delx     delt <= CFL*delx/umax
!
! viscous CFL =  delt*mu/delx^2  delt <= CFL*delx^2/mu
!  
!

if (g_nz>1) then
   umax=maxs(1)/delx+maxs(2)/dely+maxs(3)/delz
   mumax = mu*(1/delx**2 + 1/dely**2 + 1/delz**2)
else
   umax=maxs(1)/delx+maxs(2)/dely
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
! diagnostic output
!
doit=check_time(itime,time,diag_dt,0,0.0,time_next)
time_target=min(time_target,time_next)
if (doit) then
   call output_diags(time,Q)
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

   write(message,'(a,f9.5,a,i5,a,f8.5)') 'time=',time,'(',itime,')  next output=',time_target
   call print_message(message)	

   write(message,'(a,f9.7,a,f6.3,a,f6.3)') 'for next timestep: delt=',delt,' cfl_adv=',cfl_used_adv,' cfl_vis=',cfl_used_vis
   call print_message(message)	

   write(message,'(a,3f22.15)') 'max: (u,v,w) ',maxs(1),maxs(2),maxs(3)
   call print_message(message)	

   call compute_div(Q,divx,divi)
   write(message,'(3(a,e12.5))') 'max(div)=',divx,'   <z-vor>=',ints(4),'   <hel>=',ints(5)
   call print_message(message)	

   mx=ints(1)/g_nx/g_ny/g_nz  ! dont mult all the ints together, overflow!
   
   write(message,'(a,f13.10,a,f10.5,a,f10.5,a)') 'ke: ',mx,&
     ' d/dt log(ke) tot:',ints(2)/ints(1),&
     ' d/dt log(ke) diffusion:',ints(3)/ints(1)
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
character*80 message
character*20 tmp

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
open(unit=10,file=message,form='binary')



write(10) time
xnv=n_var
xnx=nslabx
xny=nslaby
xnz=nslabz
write(10) xnx,xny,xnz,xnv
do n=1,n_var
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   write(10) Q(i,j,k,n)	
enddo
enddo
enddo
enddo
close(10)


end subroutine






subroutine output_diags(time,Q)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: time

! local variables
integer i,j,k,n
integer :: iwave,iwave_max
real*8,allocatable  ::  spectrum(:),spectrum1(:)

integer,parameter :: numx=15,numy=10
character :: plot(0:numx,0:numy)
real*8 cspec(0:numx)
integer ix,iy

iwave_max=max(g_nx,g_ny,g_nz)
allocate(spectrum(0:iwave_max))
allocate(spectrum1(0:iwave_max))
spectrum=0
spectrum1=0

do i=1,3
   iwave=iwave_max
   call compute_spectrum(Q(:,:,:,i),spectrum1,iwave,io_pe)
   spectrum=spectrum+.5*spectrum1
enddo



if (my_pe==io_pe) then
! ASCII plot of KE spectrum
! number of waves:  1..iwave+1
! energy range:     1 .. 10e-10  log:  0 .. -10
!
cspec=0
plot=" "

do i=0,iwave
   if (i==0) then
      ix=0
   else
      ix=1 + ( (numx-1)*log10(real(i))/log10(real(iwave)) )     ! ranges from 1..numx
   endif
   if (ix<0) ix=0
   if (ix>numx) ix=numx
   cspec(ix)=cspec(ix)+spectrum(i)
enddo

do i=0,numx
   iy = -(numy/10.0) * log10(1d-200+cspec(i))  ! ranges from 0..numy
   if (iy<0) iy=0
   if (iy>numy) iy=numy
   cspec(i)=iy
enddo

do i=0,numx
    j=cspec(i)
    plot(i,j)="*"
enddo

print *
print *,"1E0   |",plot(:,0)," KE spectrum"
do i=1,numy-1
   print *,"      |",plot(:,i)
enddo
print *,"1E-10 |",plot(:,numy)
print *,"      +---------------+"
print *,"      k=0           k=",iwave

endif


deallocate(spectrum)
deallocate(spectrum1)

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
character*80 message
integer n_var_start

call vorticity(vor,Q,d1,work)


if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time
   message = runname(1:len_trim(runname)) // message(2:10) // ".vor"
   open(unit=11,file=message,form='binary')

   write(11) time
   xnx=o_nx
   xny=o_ny
   xnz=o_nz
   write(11) xnx,xny,xnz
   write(11) g_xcord(1:o_nx)
   write(11) g_ycord(1:o_ny)
   write(11) g_zcord(1:o_nz)
endif

call output1(vor(1,1,1,3),work,11)
if (my_pe==io_pe) close(11)

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


