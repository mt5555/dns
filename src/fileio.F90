#include "macros.h"
subroutine time_control(itime,time,Q,ints,maxs)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: time,ints(*),maxs(*)
integer :: itime

! local variables
integer i,j,k,n
character*80 message
real*8 remainder, time_target, umax,mumax,time_next,cfl_used_adv,cfl_used_vis,mx
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
umax=pi*max(maxs(1)/delx,maxs(2)/dely,maxs(3)/delz)
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
   call compute_div(Q,divx,divi)
   write(message,'(a,3f12.5,a,e12.5)') 'max: (u,v,w) ',maxs(1),maxs(2),maxs(3),' max div',divx
   call print_message(message)	

   mx = ints(3)/(1d-100+ints(2))
   if (mx>1) mx=1    ! round-off error if ke_viscous > ke_tot
   if (mx<-1) mx=-1  ! total ke is growing, but viscous term is negative
   
   write(message,'(a,f10.2,a,f6.2,a,f6.2,a)') 'ke: ',ints(1),&
     ' d/dt log(ke) total:',ints(2)/ints(1),&
     ' d/dt log(ke) diffusion:',ints(3)/ints(1)
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
write(tmp,'(f9.4)') 1000.0000 + time
message=message(1:len_trim(message)) // tmp(2:9)

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
call compute_spectrum(p,spectrum,iwave,io_pe)

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






subroutine output1(p,pt,unit)
use params
use mpi
use transpose
implicit none
real*8 :: p(nx,ny,nz)
real*8 :: pt(g_nx2,nslaby,nz_2dx)

! local vars
real*8 buf(g_nx,nslaby)
integer sending_pe,ierr,tag,unit,z_pe,y_pe,x_pe
integer request,statuses(MPI_STATUS_SIZE),dest_pe3(3)
integer i,j,k,l
integer n1,n1d,n2,n2d,n3,n3d

call transpose_to_x(p,pt,n1,n1d,n2,n2d,n3,n3d)
do z_pe=0,ncpu_z-1
do x_pe=0,ncpu_x-1
do k=1,nz_2dx
do y_pe=0,ncpu_y-1

   ! output pt(1:g_nx,1:nslaby,k) from cpus: x_pe,y_pe,z_pe
   l=g_nx*nslaby

   dest_pe3(1)=my_x
   dest_pe3(2)=my_y
   dest_pe3(3)=my_z
   tag=1
   call mpi_cart_rank(comm_3d,dest_pe3,sending_pe,ierr)

   if (sending_pe==my_pe) then

      buf(:,:)=pt(1:g_nx,1:nslaby,k)
      if (my_pe == io_pe) then
         ! dont send message to self
      else
         tag=1
         call MPI_ISend(buf,l,MPI_REAL8,io_pe,tag,comm_3d,request,ierr)
         ASSERT("output1: MPI_ISend failure",ierr==0)
         call MPI_waitall(1,request,statuses,ierr) 	
         ASSERT("output1: MPI_waitalll failure",ierr==0)
      endif
   endif

   if (my_pe==io_pe) then
      if (sending_pe==my_pe) then
         ! dont recieve message from self
      else
         call MPI_IRecv(buf,l,MPI_REAL8,sending_pe,tag,comm_3d,request,ierr)
         ASSERT("output1: MPI_IRecv failure",ierr==0)
         call MPI_waitall(1,request,statuses,ierr) 	
         ASSERT("output1: MPI_waitalll failure",ierr==0)
      endif

      write(unit) buf
   endif

enddo
enddo
enddo
enddo


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

call vorticity(vor,Q,d1,d2)

if (my_pe==io_pe) then
   write(message,'(f9.4)') 1000.0000 + time
   message = runname(1:len_trim(runname)) // message(2:9) // ".vor"
   open(unit=11,file=message,form='binary')

   write(11) time
   xnv=n_var
   n_var_start=1
   if (nz2==nz1) then
      xnv=1
      n_var_start=3
  endif

   xnx=nslabx
   xny=nslaby
   xnz=nslabz
   write(11) xnx,xny,xnz,xnv
endif

call output1(vor,work,11)
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


