#include "macros.h"



subroutine multfile_io(time,Q)
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

message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(1:len_trim(message)) // ".data"

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
call cclose(fid,ierr)


end subroutine









subroutine output_spec(time,Q,q1,q2,q3,work1,work2)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: q2(nx,ny,nz,n_var)
real*8 :: q3(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time

! local variables
integer i,j,k,n
integer :: iwave,iwave_max,ierr
real*8 spec_x(0:g_nx/2)
real*8 spec_y(0:g_ny/2)
real*8 spec_z(0:g_nz/2)
real*8 spec_x2(0:g_nx/2)
real*8 spec_y2(0:g_ny/2)
real*8 spec_z2(0:g_nz/2)
real*8 :: x,divx,divi
real*8 ::  spectrum(0:max(g_nx,g_ny,g_nz))
real*8 ::  spectrum1(0:max(g_nx,g_ny,g_nz))
character(len=80) :: message
character :: access
CPOINTER fid

! append to output files, unless time=0 create a new file 
access="a"
if (time==time_initial) access="w"

iwave_max=max(g_nx,g_ny,g_nz)
spectrum=0
spectrum1=0
spec_x=0
spec_y=0
spec_z=0

q1=Q


do i=1,ndim
   iwave=iwave_max
   call compute_spectrum(q1(1,1,1,i),work1,work2,spectrum1,spec_x2,spec_y2,spec_z2,iwave,io_pe)
   spectrum=spectrum+.5*spectrum1
   spec_x=spec_x + .5*spec_x2
   spec_y=spec_y + .5*spec_y2
   spec_z=spec_z + .5*spec_z2
enddo
write(message,'(a,f10.4)') " KE spectrum",time
call plotASCII(spectrum,iwave,message(1:25))
!call plotASCII(spec_x,g_nx/2,message)
!call plotASCII(spec_y,g_ny/2,message)
!call plotASCII(spec_z,g_nz/2,message)



! for incompressible equations, print divergence as diagnostic:
if (equations==NS_UVW) then
   call compute_div(Q,q1,work1,work2,divx,divi)
   write(message,'(3(a,e12.5))') 'max(div)=',divx
   call print_message(message)	
endif



if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time_initial
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".spec"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "spec_write(): Error opening file errno=",ierr
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
   call cclose(fid,ierr)
endif
end subroutine





subroutine output_scalars(time,ints_save,maxs_save,nv,nscalars)
use params
implicit none
real*8 :: time
integer nv,nscalars
real*8 :: ints_save(nv,nscalars)
real*8 :: maxs_save(nv,nscalars)

! local variables
real*8 :: x
integer i,j,k,n,ierr
character(len=80) :: message
character :: access
CPOINTER fid

access="a"
if (time==time_initial) access="w"


if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time_initial
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".scalars"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      call print_message(message)
      write(message,'(a,i5)') "diag_output(): Error opening .scalars file errno=",ierr
      call abort(message)
   endif
   x=nv; call cwrite8(fid,x,1)
   x=nscalars; call cwrite8(fid,x,1)
   call cwrite8(fid,mu,1)
   call cwrite8(fid,alpha_value,1)
   call cwrite8(fid,ints_save,nv*nscalars);
   call cwrite8(fid,maxs_save,nv*nscalars);

   x=0; call cwrite8(fid,x,1)  ! for historical file format reasons
   call cwrite8(fid,time,1)

   call cclose(fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "diag_output(): Error closing .scalars file errno=",ierr
      call abort(message)
   endif
endif

end subroutine











subroutine singlefile_io(time,p,fname,work,work2,read,fpe)
!
! I/O routines where all data goes through a single PE and is
! written to a single file
!
! read=0    write data to file fname
! read=1    read data from file fname
!
! fpe       processor to do the file I/O
!
use params
use mpi
use transpose
implicit none
integer :: read  ! =1 for read, 0 for write
integer :: fpe
real*8 :: time
real*8 :: p(nx,ny,nz)
real*8 :: work2(nx,ny,nz),work(nx,ny,nz)
character(len=*) :: fname

! local variables
integer i,j,k,n
real*8 xnx,xny,xnz
character(len=80) message
integer n_var_start,ierr
CPOINTER fid


if (my_pe==fpe) then

   if (read==1) then
      call copen(fname,"r",fid,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "singlefile_io(): Error opening file. Error no=",ierr
         call print_message(message)
         call print_message(fname)
         call abort("")
      endif
      call cread8e(fid,time,1,ierr)
      if (ierr/=1) then
         write(message,'(a,i5)') "singlefile_io(): Error reading file"
         call print_message(message)
         call print_message(fname)
         call abort("")
      endif
      xnx=o_nx
      xny=o_ny
      xnz=o_nz
      call cread8(fid,xnx,1)
      call cread8(fid,xny,1)
      call cread8(fid,xnz,1)
      if (int(xnx)/=o_nx) call abort("Error: data file nx <> nx set in params.h");
      if (int(xny)/=o_ny) call abort("Error: data file ny <> ny set in params.h");
      if (int(xnz)/=o_nz) call abort("Error: data file nz <> nz set in params.h");
      call cread8(fid,g_xcord(1),o_nx)
      call cread8(fid,g_ycord(1),o_ny)
      call cread8(fid,g_zcord(1),o_nz)
   else
      call copen(fname,"w",fid,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "singlefile_io(): Error opening file. Error no=",ierr
         call print_message(message)
         call print_message(fname)
         call abort("")
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
endif

#ifdef USE_MPI
call MPI_bcast(time,1,MPI_REAL8,io_pe,comm_3d ,ierr)
#endif

if (read==1) then
   call input1(p,work,work2,fid,fpe,.false.)
else
   call output1(p,work,work2,fid,fpe)
endif
if (my_pe==fpe) call cclose(fid,ierr)

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
   time_next = min(time_next,time-remainder+dt)

endif

end function


