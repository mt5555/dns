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
real*8,allocatable  :: q1(:,:,:,:)
real*8,allocatable  :: q2(:,:,:,:)
real*8,allocatable  :: q3(:,:,:,:)
character(len=80) message,sdata
character(len=280) basename,fname
integer ierr,i,j,k,n,km,im,jm,icount
real*8 :: tstart,tstop,tinc,time,time2
real*8 :: u,v,w,x,y
real*8 :: kr,ke,ck,xfac
integer :: lx1,lx2,ly1,ly2,lz1,lz2,nxlen,nylen,nzlen
integer :: nxdecomp,nydecomp,nzdecomp
CPOINTER :: fid

tstart=3.25
tstop=3.25
tinc=1.0
icount=0

nxdecomp=1
nydecomp=1
nzdecomp=1

!call set_byteswap_input(1);


! these lines are modifed by some sed scripts for automatic running
! of this code by putting in new values of tstart, tstop, tinc,
! nxdecomp,nydecomp,nzdecom, etc.
!SEDbyteswap
!SEDtstart
!SEDdecomp



call init_mpi
call init_mpi_comm3d()
call init_model

!call writepoints(); stop


if (nxdecomp*nydecomp*nzdecomp>1) then
   ! in this case, we only need 2 work areas.  we are running in
   ! serial, so we have to save as much memory as possible
   allocate(q1(nx,ny,nz,1))
   allocate(q2(nx,ny,nz,1))
else
   allocate(q1(nx,ny,nz,ndim))
   allocate(q2(nx,ny,nz,ndim))
   allocate(q3(nx,ny,nz,ndim))
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if needed, initialize some constants.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Q=0



time=tstart
do
   icount=icount+1

   write(sdata,'(f10.4)') 10000.0000 + time
   fname = runname(1:len_trim(runname)) // sdata(2:10) // ".u"
   call print_message(fname(1:len_trim(fname)))
   call singlefile_io(time2,Q(1,1,1,1),fname,q1,q2,1,io_pe)

   fname = runname(1:len_trim(runname)) // sdata(2:10) // ".v"
   call print_message(fname(1:len_trim(fname)))
   call singlefile_io(time2,Q(1,1,1,2),fname,q1,q2,1,io_pe)

   fname = runname(1:len_trim(runname)) // sdata(2:10) // ".w"
   call print_message(fname(1:len_trim(fname)))
   call singlefile_io(time2,Q(1,1,1,3),fname,q1,q2,1,io_pe)
   time=time2

   do k=nz1,nz2
      do j=ny1,ny2
      do i=nx1,nx2
         Q(i,j,k,1)=i
         Q(i,j,k,2)=j
         Q(i,j,k,3)=k
      enddo
      enddo
      enddo



   do i=0,nxdecomp-1
   do j=0,nydecomp-1
   do k=0,nzdecomp-1  
      nxlen=nslabx/nxdecomp
      nylen=nslaby/nydecomp
      nzlen=nslabz/nzdecomp
      lx1=nx1 + i*nxlen     
      ly1=ny1 + j*nylen     
      lz1=nz1 + k*nzlen     

      lx2=lx1 + nxlen-1
      ly2=ly1 + nylen-1
      lz2=lz1 + nzlen-1


      if (nxdecomp*nydecomp*nzdecomp==1) then
         ! no subcubes:
         call isoavep(Q,q1,q2,q3)
         ! debug version:
         ! print *,'calling 1 cpu verison for debug'
         ! call isoave1(Q,q1,q2,nx1,nx2,ny1,ny2,nz1,nz2)
      else
         ! subcube version cannot run in parallel - it will abort
         call isoave1(Q,q1,q2,lx1,lx2,ly1,ly2,lz1,lz2)
      endif


      if (my_pe==io_pe) then
         write(sdata,'(f10.4)') 10000.0000 + time
         fname = runname(1:len_trim(runname)) // sdata(2:10) // ".isostr"
         if (nxdecomp*nydecomp*nzdecomp>1) then
            write(sdata,'(3i1)') i,j,k
            fname=fname(1:len_trim(fname)) // sdata(1:3)
         endif

         print *,fname
         call copen(fname,"w",fid,ierr)
         if (ierr/=0) then
            write(message,'(a,i5)') "output_model(): Error opening .sf file errno=",ierr
            call abort(message)
         endif
         call writeisoave(fid,time)
         call cclose(fid,ierr)
      endif

   enddo
   enddo
   enddo


   time=time+tinc
   if (time>tstop) exit
enddo

call close_mpi
end program anal


