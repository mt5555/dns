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
character,save :: access="0"
CPOINTER fid

! append to output files, unless this is first call.
if (access=="0") then
   access="w"
else
   access="a"
endif


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






subroutine singlefile_io(time,p,fname,work,work2,io_read,fpe)
!
! I/O routines where all data goes through a single PE and is
! written to a single file
!
! io_read=0    write data to file fname
! io_read=1    read data from file fname
!
! fpe       processor to do the file I/O
!
use params
use mpi
use transpose
implicit none
integer :: io_read  ! =1 for read, 0 for write
integer :: fpe
real*8 :: time
real*8 :: p(nx,ny,nz)
real*8 :: work2(nx,ny,nz),work(nx,ny,nz)
character(len=*) :: fname

call singlefile_io2(time,p,fname,work,work2,io_read,fpe,.false.)

end subroutine




subroutine singlefile_io2(time,p,fname,work,work2,io_read,fpe,output_spec)
!
! I/O routines where all data goes through a single PE and is
! written to a single file
!
! io_read=0    write data to file fname
! io_read=1    read data from file fname
!
! fpe       processor to do the file I/O
!
! output_spec=.true.     i/o on 2/3 dealiased spectral coefficients
! output_spec=.false.    i/o on full 3d array p 
!                     (grid point data or non-dealiased spec coeff.)
!
use params
use mpi
use transpose
implicit none
integer :: io_read  ! =1 for read, 0 for write
integer :: fpe
real*8 :: time
real*8 :: p(nx,ny,nz)
real*8 :: work2(nx,ny,nz),work(nx,ny,nz)
character(len=*) :: fname
logical :: output_spec
integer im_max,km_max,jm_max



! local variables
integer i,j,k,n
real*8 xnx,xny,xnz
character(len=80) message
integer n_var_start,ierr
CPOINTER fid

if (use_mpi_io .and. io_read==0) then
   call singlefile_mpi_io(time,p,fname,work,work2,io_read,fpe)
   return
endif


xnx=o_nx
xny=o_ny
xnz=o_nz


if (output_spec) then

   ! default is to write everything:
   xnx=g_nx
   xny=g_ny
   xnz=g_nz

   if (dealias==1) then
      xnx=2+2*(g_nx/3)
      xny=2+2*(g_ny/3)
      xnz=2+2*(g_nz/3)
   endif

   ! if some kind of spectral truncation:


endif

if (my_pe==fpe) then

   if (io_read==1) then
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
      call cread8(fid,xnx,1)
      call cread8(fid,xny,1)
      call cread8(fid,xnz,1)
      if (output_spec) then
         print *,'Spectrum input data'
         print *,'number of real coefficients: ',xnx,xny,xnz
         if (xnx>g_nx .or. xny>g_ny .or. xnz>g_nz) then
            ! we can only upsample low-res data.
            ! to run with high-res data, output a trucated form.  
            call print_message("Error: spectral input requires downsampling to lower resolution")
            call print_message("Input routines can only upsample.") 
            call print_message("Output routines can only downsample.") 
            call print_message("Run code at higher resolution, calling Output to downsample")
            call abort("error in singlefile_io2")
         endif
      else
         if (int(xnx)/=o_nx) call abort("Error: data file nx <> nx set in params.h");
         if (int(xny)/=o_ny) call abort("Error: data file ny <> ny set in params.h");
         if (int(xnz)/=o_nz) call abort("Error: data file nz <> nz set in params.h");
         call cread8(fid,g_xcord(1),o_nx)
         call cread8(fid,g_ycord(1),o_ny)
         call cread8(fid,g_zcord(1),o_nz)
      endif
   else
      call copen(fname,"w",fid,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "singlefile_io(): Error opening file. Error no=",ierr
         call print_message(message)
         call print_message(fname)
         call abort("")
      endif
      call cwrite8(fid,time,1)
      call cwrite8(fid,xnx,1)
      call cwrite8(fid,xny,1)
      call cwrite8(fid,xnz,1)
      if ( output_spec) then
         if (xnx>g_nx .or. xny>g_ny .or. xnz>g_nz) then
            ! we can only upsample low-res data.
            ! to run with high-res data, output a trucated form.  
            call print_message("Error: spectral output requires zero padding") 
            call print_message("Output routines can only downsample.") 
            call print_message("Input routines can input this data directly")
            call abort("error in singlefile_io2")
         endif
      else
         call cwrite8(fid,g_xcord(1),o_nx)
         call cwrite8(fid,g_ycord(1),o_ny)
         call cwrite8(fid,g_zcord(1),o_nz)
      endif
   endif
endif

#ifdef USE_MPI
call MPI_bcast(time,1,MPI_REAL8,io_pe,comm_3d ,ierr)
#endif

! max wave number to read/write for spectral I/O
im_max=xnx/2-1
jm_max=xny/2-1
km_max=xnz/2-1


if (io_read==1) then
   if (output_spec) then
      call input1_spec(p,work,work2,fid,fpe,im_max,jm_max,km_max)
   else
      call input1(p,work,work2,fid,fpe,.false.)
   endif
else
   if (output_spec) then
      call output1_spec(p,work,work2,fid,fpe,im_max,jm_max,km_max)
   else
      call output1(p,work,work2,fid,fpe,-1)
   endif

endif
if (my_pe==fpe) call cclose(fid,ierr)

end subroutine





subroutine singlefile_mpi_io(time,p,fname,work,work2,io_read,fpe)
!
! I/O routines where all data goes through a single PE and is
! written to a single file
!
! io_read=0    write data to file fname
! io_read=1    read data from file fname
!
! fpe       processor to do the file I/O
!
! output_spec=.true.     i/o on 2/3 dealiased spectral coefficients
! output_spec=.false.    i/o on full 3d array p 
!                     (grid point data or non-dealiased spec coeff.)
!
use params
use mpi
use transpose
implicit none
integer :: io_read  ! =1 for read, 0 for write
integer :: fpe
real*8 :: time
real*8 :: p(nx,ny,nz)
real*8 :: work2(nx,ny,nz),work(nx,ny,nz)
character(len=*) :: fname
integer im_max,km_max,jm_max



! local variables
integer i,j,k,n
real*8 xnx,xny,xnz
character(len=80) message
integer :: n_var_start,ierr,offset=0

integer statuses(MPI_STATUS_SIZE)
integer*8 infoin

CPOINTER fid

#ifndef USE_MPI_IO
call abort("singilefile_mpi_io: error: code not compiled with MPI-IO")   
#else


xnx=o_nx
xny=o_ny
xnz=o_nz

if (my_pe==fpe) then

   if (io_read==1) then
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
      call MPI_Info_create(infoin,ierr)
      call MPI_Info_set(infoin, "striping_factor", "64",ierr) 	
      call MPI_Info_set(infoin, "striping_unit",  "8388608",ierr)
      call MPI_File_open(comm_1,fname, &
           MPI_MODE_WRONLY + MPI_MODE_CREATE ,&
           infoin, fid,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "singlefile_io(): Error opening file. Error no=",ierr
         call print_message(message)
         call print_message(fname)
         call abort("")
      endif
      call MPI_File_write(fid,time,1,MPI_REAL8,statuses,ierr)
      call MPI_File_write(fid,xnx,1,MPI_REAL8,statuses,ierr)
      call MPI_File_write(fid,xny,1,MPI_REAL8,statuses,ierr)
      call MPI_File_write(fid,xnz,1,MPI_REAL8,statuses,ierr)
      call MPI_File_write(fid,g_xcord(1),o_nx,MPI_REAL8,statuses,ierr)
      call MPI_File_write(fid,g_ycord(1),o_ny,MPI_REAL8,statuses,ierr)
      call MPI_File_write(fid,g_zcord(1),o_nz,MPI_REAL8,statuses,ierr)
      offset=4 + o_nx + o_ny + o_nz
   endif
endif

#ifdef USE_MPI
call MPI_bcast(time,1,MPI_REAL8,io_pe,comm_3d ,ierr)
#endif

! max wave number to read/write for spectral I/O
im_max=xnx/2-1
jm_max=xny/2-1
km_max=xnz/2-1


if (io_read==1) then
   call input1(p,work,work2,fid,fpe,.false.)
else
   call output1(p,work,work2,fid,fpe,offset)
endif
if (my_pe==fpe) then
   call MPI_File_close(fid,ierr)
endif
#endif
end subroutine

















