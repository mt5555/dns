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
real*8,allocatable  :: Q(:,:,:,:)
real*8,allocatable  :: vor(:,:,:,:)
real*8,save  :: work1(nx,ny,nz)
real*8,save  :: work2(nx,ny,nz)
character(len=80) message,sdata,tname
character(len=280) basename,fname
integer ierr,i,j,k,n,km,im,jm,icount
real*8 :: tstart,tstop,tinc,time,time2
real*8 :: u,v,w,x,y
real*8 :: kr,ke,ck,xfac,dummy
character(len=3) :: extension="uvwX"

! input file
tstart=.32
tstop=.32
tinc=1.0

! to read times from  file times.dat:
! tstart=-1; tinc=0; tname="times.dat"

! these lines are modifed by some sed scripts for automatic running
! of this code by putting in new values of tstart, tstop, tinc,
! nxdecomp,nydecomp,nzdecom, etc.
!SEDtstart


icount=0
call init_mpi       
call init_mpi_comm3d()
call init_model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if needed, initialize some constants.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (convert_opt == 3) then
   allocate(vor(1,1,1,1)) ! dummy variable -wont be used
   allocate(Q(nx,ny,nz,n_var))
else if (convert_opt == 4) then
   allocate(vor(1,1,1,1)) ! dummy variable -wont be used
   allocate(Q(nx,ny,nz,1))
else
   allocate(vor(nx,ny,nz,n_var))
   allocate(Q(nx,ny,nz,n_var))
endif


time=tstart
do
icount=icount+1
   if (tstart<0) then
      ! read times from unit 83
      fname= rundir(1:len_trim(rundir)) // tname(1:len_trim(tname))
      if (icount==1)  open(83,file=fname)
      read(83,*,err=100,end=100) time
   endif	

   Q=0

   if (convert_opt==0) then  ! -cout uvw  
      call input_uvw(time,Q,vor,work1,work2)
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
      call input_uvw(time,Q,vor,work1,work2)
      ! outputing vorticity
      basename=runname(1:len_trim(runname)) // "-vor."

      call print_message("computing vorticity...")
      call vorticity(vor,Q,work1,work2)
      call print_message("output vorticity...")
      call output_uvw(basename,time,vor,Q,work1,work2)
   endif

   if (convert_opt==2) then  ! -cout vorm
      call input_uvw(time,Q,vor,work1,work2)
      call print_message("computing vorticity magnitude...")
      call vorticity(vor,Q,work1,work2)
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
   if (convert_opt==3) then  ! -cout norm2
      ! 2048^3 needs 192GB * 1.66.  needs 256 cpus
      call input_uvw(time,Q,vor,work1,work2)
      print *,'max input: ',maxval(Q(nx1:nx2,ny1:ny2,nz1:nz2,1)), &
                            maxval(Q(nx1:nx2,ny1:ny2,nz1:nz2,2)), &
                            maxval(Q(nx1:nx2,ny1:ny2,nz1:nz2,3))
      call print_message("computing norm squared...")
      do k=nz1,nz2
      do j=ny1,ny2
      do i=nx1,nx2
         work1(i,j,k)=Q(i,j,k,1)**2+Q(i,j,k,2)**2+Q(i,j,k,3)**2
      enddo
      enddo
      enddo
      call print_message("outputting norm squared as REAL*4...")
      output_size=4
      write(sdata,'(f10.4)') 10000.0000 + time
      basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))
      fname = basename(1:len_trim(basename)) // sdata(2:10) // ".norm2"
      call singlefile_io3(time,work1,fname,Q,work2,0,io_pe,.false.,2)
      fname = basename(1:len_trim(basename)) // sdata(2:10) // ".norm2e"
      call singlefile_io3(time,work1,fname,Q,work2,0,io_pe,.false.,3)
   endif
   if (convert_opt==4) then  ! -cout 4uvw
      ! 2048^3 needs 192GB storage.  can run on 128 cpus.
      !  can run on 64 cpus, if running 2 jobs per node.
      write(sdata,'(f10.4)') 10000.0000 + time
      basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))

      do i=1,3  ! use this to loop over .u,.v,.w
         fname = basename(1:len_trim(basename)) // sdata(2:10) // &
                 "." // extension(i:i)
         call print_message("input file:")	
	 call print_message(fname(1:len_trim(fname)))
         call singlefile_io3(time,Q,fname,work1,work2,1,io_pe,.false.,1)
         print *,'max input: ',maxval(Q(nx1:nx2,ny1:ny2,nz1:nz2,1))

	 call print_message("outputting as REAL*4...")
         output_size=4
         fname = basename(1:len_trim(basename)) // sdata(2:10) // &
                 "." // extension(i:i) // "4"
	 call print_message(fname(1:len_trim(fname)))

         call singlefile_io3(time,Q,fname,work1,work2,0,io_pe,.false.,2)
      enddo	

   endif

   if (tstart>0) then
      time=time+tinc
      if (time > max(tstop,tstart)) exit
      if (time < min(tstop,tstart)) exit
   endif
enddo
if (io_pe==my_pe) then
   print *,'convert.F90 finished t=',time
   print *,tstart,tstop
endif

100 continue
call close_mpi
end program


