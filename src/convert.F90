#include "macros.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read in a sequence of data files and perform analysis on them
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
program convert
use params
use mpi
use fft_interface
implicit none
real*8,save  :: Q(nx,ny,nz,n_var)
real*8,save  :: vor(nx,ny,nz,n_var)
real*8,save  :: work1(nx,ny,nz)
real*8,save  :: work2(nx,ny,nz)
character(len=80) message,sdata
character(len=280) basename,fname
integer ierr,i,j,k,n,km,im,jm,icount
real*8 :: tstart,tstop,tinc,time,time2
real*8 :: u,v,w,x,y
real*8 :: kr,ke,ck,xfac,dummy

! input file
tstart=.32
tstop=.32
tinc=1.0
icount=0
call init_mpi       
call init_mpi_comm3d()
call init_model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if needed, initialize some constants.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Q=0


time=tstart
do
icount=icount+1

   call input_uvw(time,Q,vor,work1,work2)
   print *,'max U: ',&
        maxval(Q(nx1:nx2,ny1:ny2,nz1:nz2,1)), &
        maxval(Q(nx1:nx2,ny1:ny2,nz1:nz2,2)), &
        maxval(Q(nx1:nx2,ny1:ny2,nz1:nz2,3))

   if (convert_opt==0) then  ! -cout uvw  
      ! just reoutput the variables:
      if (w_spec) then
         do n=1,3
            call fft3d(Q(1,1,1,n),work1)
         enddo
      endif
      basename=runname(1:len_trim(runname)) // "-new."
      call output_uvw(basename,time,Q,vor,work1,work2)
   endif

   if (convert_opt==1) then  ! -cout vor
      ! outputing vorticity
      basename=runname(1:len_trim(runname)) // "-vor."

      call print_message("computing vorticity...")
      call vorticity(vor,Q,work1,work2)
      call print_message("output vorticity...")
      call output_uvw(basename,time,vor,Q,work1,work2)
   endif

   if (convert_opt==2) then  ! -cout vorm
      call print_message("computing vorticity magnitude...")
      call vorticity(vor,Q,work1,work2)
      vor=Q
      do k=nz1,nz2
      do j=ny1,ny2
      do i=nx1,nx2
         work1(i,j,k)=sqrt(vor(i,j,k,1)**2+vor(i,j,k,2)**2+vor(i,j,k,3)**2)
      enddo
      enddo
      enddo
      ! output vorticity magnitude
      write(sdata,'(f10.4)') 10000.0000 + time
      basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))
      fname = basename(1:len_trim(basename)) // sdata(2:10) // ".vorm"
      call singlefile_io2(time,work1,fname,vor,work2,0,io_pe,.false.)
   endif

   time=time+tinc
   if (time > max(tstop,tstart)) exit
   if (time < min(tstop,tstart)) exit
enddo


call close_mpi
end program


