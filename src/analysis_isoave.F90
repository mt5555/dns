#include "macros.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! compute structure functions for many different directions
! in the periodic cube.
! must be run on only 1 cpu.
!
!
! To run, set the base name of the file and the times of interest
! below.  For example:
!    tstart=0
!    tstop=1
!    tinc=.5
!    basename="temp"
!
! will result in looping over the files:   
!             temp0000.0000.[uvw]
!             temp0000.5000.[uvw]
!             temp0001.0000.[uvw]
!
!  to compile and run:   make analysis ; analysis
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program anal
use params
use mpi
use isoave
implicit none
real*8,save  :: Q(nx,ny,nz,n_var)
real*8,save  :: work1(nx,ny,nz)
real*8,save  :: work2(nx,ny,nz)
character(len=80) message,sdata
character(len=280) basename,fname
integer ierr,i,j,k,n,km,im,jm,icount
real*8 :: tstart,tstop,tinc,time,time2
real*8 :: u,v,w,x,y
real*8 :: kr,ke,ck,xfac
CPOINTER :: fid

! input file
basename="temp4_"
print *,basename
tstart=1.5
tstop=11.0
tinc=1.0
icount=0

!call set_byteswap_input(1);


call init_mpi
call init_mpi_comm3d()
call init_model

call writepoints()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if needed, initialize some constants.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Q=0


time=tstart
do
   icount=icount+1
   print *,'icount = ',icount

   write(sdata,'(f10.4)') 10000.0000 + time
   fname = basename(1:len_trim(basename)) // sdata(2:10) // ".u"
   print *,'filename: ',fname(1:len_trim(fname))
   call singlefile_io(time2,Q(1,1,1,1),fname,work1,work2,1,io_pe)

   fname = basename(1:len_trim(basename)) // sdata(2:10) // ".v"
   print *,'filename: ',fname(1:len_trim(fname))
   call singlefile_io(time2,Q(1,1,1,2),fname,work1,work2,1,io_pe)

   fname = basename(1:len_trim(basename)) // sdata(2:10) // ".w"
   print *,'filename: ',fname(1:len_trim(fname))
   call singlefile_io(time2,Q(1,1,1,3),fname,work1,work2,1,io_pe)


  
   call isoave1(Q)


   write(sdata,'(f10.4)') 10000.0000 + time
   fname = basename(1:len_trim(basename)) // sdata(2:10) // ".isostr"
   call copen(fname,"w",fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "output_model(): Error opening .sf file errno=",ierr
      call abort(message)
   endif
   call writeisoave(fid)
   call cclose(fid,ierr)



   time=time+tinc
   if (time>tstop) exit
enddo

call close_mpi
end program anal


