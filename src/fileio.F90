#include "macros.h"
subroutine time_control(itime,time,Q,ints)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: time,ints(3)
integer :: itime

! local variables
integer i,j,k,n
character*80 message
real*8 remainder, time_target, maxv(3),umax,mumax,time_next,cfl_used_adv,cfl_used_vis,mx
real*8 divx,divi
logical,external :: check_time
logical :: doit

time_target = time_final

!
! compute new time step
!
! CFL = delt*umax/delx     delt <= CFL*delx/umax
!
! viscous CFL =  delt*mu/delx^2  delt <= CFL*delx^2/mu
!  
!
do i=1,3
   maxv(i) = maxval( abs(Q(nx1:nx2,ny1:ny2,nz1:nz2,i)) )
enddo
umax=pi*max(maxv(1)/delx,maxv(2)/dely,maxv(3)/delz)
#ifdef MPI
   call MPI_REDUCE(umax,MPI_MAX ... )
#endif
if (umax> pi*1000/min(delx,dely,delz)) error_code=1

! advective CFL
delt = cfl_adv/umax

! viscous CFL
mumax=min(delx,dely,delz)
mumax = pi*pi*mu/(mumax*mumax)
delt = min(delt,cfl_vis/mumax)


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




!
! restrict delt so we hit the next time_target
!
doit=check_time(itime,time,screen_dt,0,0.0,time_next)
time_target=min(time_target,time_next)

delt = min(delt,time_target-time)
cfl_used_adv=umax*delt
cfl_used_vis=mumax*delt


!
! display screen output
!
if (doit) then

   write(message,'(a,f8.5,a,i5,a,f8.5)') 'time=',time,'(',itime,')  next output=',time_target
   call print_message(message)	

   write(message,'(a,f9.7,a,f6.3,a,f6.3)') 'delt=',delt,' cfl_adv=',cfl_used_adv,' cfl_vis=',cfl_used_vis
   call print_message(message)	

   call compute_div(Q,io_pe,divx,divi)
   write(message,'(a,3f12.5,a,e12.5)') 'max: (u,v,w) ',maxv(1),maxv(2),maxv(3),' max div',divx
   call print_message(message)	

   mx = ints(3)/(1d-100+ints(2))
   if (mx>1) mx=1    ! round-off error if ke_viscous > ke_tot
   if (mx<-1) mx=-1  ! total ke is growing, but viscous term is negative
   
   write(message,'(a,f10.2,a,f6.2,a,f6.2,a)') 'ke: ',ints(1),&
     ' d/dt log(ke) total:',ints(2)/ints(1),&
     ' d/dt log(ke) diff:',ints(3)/ints(1)
   call print_message(message)	

   call print_message("")
endif

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

write(message,'(f9.4)') 1000.0000 + time
message = "data" // message(2:9) // ".out"
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
real*8 :: work(nx,ny,nz)
real*8 :: p(nx,ny,nz)
integer i,j,k,n
integer :: iwave
real*8,allocatable  ::  spectrum(:)

iwave=max(g_nx,g_ny,g_nz)
allocate(spectrum(0:iwave))

p=0
do i=1,3
   p=p+.5*Q(:,:,:,i)*Q(:,:,:,i)  
enddo
call compute_spectrum(p,spectrum,iwave)

#if 0
iwave=max(g_nx,g_ny,g_nz)/2
do i=1,iwave
   print *,i,spectrum(i)
enddo
rwave=0
do i=iwave+1,iwave_max
   rwave=rwave+spectrum(i)
enddo
print *,iwave+1,"..",iwave_max,rwave
#endif

deallocate(spectrum)

end subroutine






subroutine output_write(time,Q)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: time

! local variables
integer i,j,k,n
real*8 xnx,xny,xnz,xnv
real*8 :: vor(nx,ny,nz,n_var)
real*8 :: d1(nx,ny,nz),work(nx,ny,nz)
real*8 :: d2(1)
character*80 message
integer n_var_start

vor=0
! compute viscous terms (in rhs) and vorticity
do i=1,3
   ! compute u_x, u_xx
   call der(Q(1,1,1,i),d1,d2,work,DX_ONLY,1)
   if (i==3) vor(:,:,:,2) = vor(:,:,:,2) - d1
   if (i==2) vor(:,:,:,3) = vor(:,:,:,3) + d1

   ! compute u_y, u_yy
   call der(Q(1,1,1,i),d1,d2,work,DX_ONLY,2)
   if (i==3) vor(:,:,:,1) = vor(:,:,:,1) + d1
   if (i==1) vor(:,:,:,3) = vor(:,:,:,3) -d1

   ! compute u_z, u_zz
   call der(Q(1,1,1,i),d1,d2,work,DX_ONLY,3)
   if (i==2) vor(:,:,:,1) = vor(:,:,:,1) -d1
   if (i==1) vor(:,:,:,2) = vor(:,:,:,2) +d1
enddo



write(message,'(f9.4)') 1000.0000 + time
message = "vor" // message(2:9) // ".out"
open(unit=10,file=message,form='binary')




write(10) time
xnv=n_var
n_var_start=1
if (nz2==nz1) then
   xnv=1
   n_var_start=3
endif

xnx=nslabx
xny=nslaby
xnz=nslabz
write(10) xnx,xny,xnz,xnv
do n=n_var_start,n_var
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   write(10) vor(i,j,k,n)	
enddo
enddo
enddo
enddo
close(10)

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


