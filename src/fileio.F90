#include "macros.h"
subroutine time_control(time,Q)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: time

! local variables
integer i,j,k,n
character*80 message
real*8 remainder, time_target, maxv(3),umax,mumax,time_next,cfl_used,mx
logical,external :: check_time
real*8,external :: norm_divergence
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
umax=max(maxv(1)/delx,maxv(2)/dely,maxv(3)/delz)
#ifdef MPI
   call MPI_REDUCE(umax,MPI_MAX ... )
#endif
if (umax> 1000/min(delx,dely,delz)) error_code=1

! advective CFL
delt = cfl_adv/umax

! viscous CFL
mumax=min(delx,dely,delz)
mumax = mu/(mumax*mumax)
delt = min(delt,cfl_vis/mumax)


delt = max(delt,delt_min)
delt = min(delt,delt_max)

!
!  restart dumps
!
doit=check_time(time,restart_dt,0,0.0,time_next)
time_target=min(time_target,time_next)
if (doit) then
   call restart_write(time,Q)
endif


!
!  output dumps
!
doit=check_time(time,output_dt,0,0.0,time_next)
time_target=min(time_target,time_next)
if (doit) then
   call output_write(time,Q)
endif





!
! restrict delt so we hit the next time_target
!
doit=check_time(time,diag_dt,0,0.0,time_next)
time_target=min(time_target,time_next)

delt = min(delt,time_target-time)
cfl_used=umax*delt

!
! display diagnostic output
!
if (doit) then
   mx=norm_divergence(Q)

   write(message,'(a,f15.10,a,f12.7,a,f5.3,a,f10.5)') 'time=',time,' delt=',delt, &
   ' adv cfl=',cfl_used,' next output=',time_target
   call print_message(message)	
   write(message,'(a,4f12.5)') 'max: (u,v,w) div',maxv(1),maxv(2),maxv(3),mx
   call print_message(message)	

   mx=.5*sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,1:3)*Q(nx1:nx2,ny1:ny2,nz1:nz2,1:3))
   mx=mx/(nslabx*nslaby*nslabz)
   write(message,'(a,4f15.10)') 'ke:',mx
   call print_message(message)	

   call print_message("")
endif

end subroutine



real*8 function norm_divergence(Q)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: p(nx,ny,nz)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: mx

call divergence(p,Q,work1,work2)
mx=maxval(abs(p(nx1:nx2,ny1:ny2,nz1:nz2)))
#ifdef MPI
   call MPI_REDUCE(mx,MPI_MAX)
#endif
norm_divergence=mx
return
end function




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















logical function check_time(time,dt,ncustom,custom,time_next)
use params
implicit none
integer ncustom
real*8 :: time,dt,custom(ncustom)
real*8 :: time_next  ! output

!local variables
real*8 remainder
integer i
real*8 :: small=1e-5

check_time = .false.
time_next=time_final

if (dt<0) then
   check_time=.true.
endif

if (dt>0) then
   remainder=mod(time,dt)
   if (remainder < small) then
      check_time=.true.
   endif
   time_next = time-remainder+dt
   if (time>=time_final-small) check_time=.true.

   do i=1,ncustom
      if (abs(time-custom(i))<small) then
         check_time=.true.
      else if (time<custom(i)) then
         time_next=min(time_next,custom(i))
         exit 
      endif
   enddo
endif

end function


