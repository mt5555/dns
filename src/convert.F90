#include "macros.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read in a sequence of data files and convert them
!
! type of conversion controled by -cout argument:
!
!  -cout uvw      (can be used with -si/-so  spectral input/output)
!                  and -smax spectral coefficient truncation to 
!                  perform downsampling or upsampling)
!  -count 4uvw     as above, but real*4 output
!  -cout vor
!  -cout vorm
!  -cout norm2
!  -cout passive    convert passive scalar file
!                   need to specify shmidt_in and type_in below
!  -cout gradu      output u_i,j
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
real*8 :: schmidt_in,mn,mx
integer :: type_in
character(len=4) :: extension="uvwX"
character(len=8) :: ext2,ext

! input file
tstart=.8
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

if (convert_opt==0 .or. convert_opt == 3 .or. convert_opt==5) then
   allocate(vor(1,1,1,1)) ! dummy variable -wont be used
   allocate(Q(nx,ny,nz,n_var))
else if (convert_opt == 4 .or. convert_opt==6) then
   allocate(vor(1,1,1,1)) ! dummy variable -wont be used
   allocate(Q(nx,ny,nz,1))
else if (convert_opt == 7) then
   allocate(vor(nx,ny,nz,1)) ! only first component used
   allocate(Q(nx,ny,nz,n_var))
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
      call input_uvw(time,Q,vor,work1,work2,1)  ! default headers
!      print *,'attempting to read headerless input data...'
!      call input_uvw(time,Q,vor,work1,work2,2)  ! no headers
      ! just reoutput the variables:
      if (w_spec) then
         do n=1,3
            call fft3d(Q(1,1,1,n),work1)
         enddo
      endif
      basename=runname(1:len_trim(runname)) // "-new."
!      call output_uvw(basename,time,Q,vor,work1,work2,1)  ! default headers
      basename=runname(1:len_trim(runname)) // "-raw."
      call output_uvw(basename,time,Q,vor,work1,work2,2)   ! no headers
   endif

   if (convert_opt==1) then  ! -cout vor
      call input_uvw(time,Q,vor,work1,work2,1)
      ! outputing vorticity
      basename=runname(1:len_trim(runname)) // "-vor."

      call print_message("computing vorticity...")
      call vorticity(vor,Q,work1,work2)
      call print_message("output vorticity...")
      call output_uvw(basename,time,vor,Q,work1,work2,1)
   endif

   if (convert_opt==2) then  ! -cout vorm
      call input_uvw(time,Q,vor,work1,work2,1)
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
      output_size=4
      write(sdata,'(f10.4)') 10000.0000 + time
      basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))

      fname = basename(1:len_trim(basename)) // sdata(2:10) // ".vorm"
      call singlefile_io3(time,work1,fname,vor,work2,0,io_pe,.false.,2)
!      fname = basename(1:len_trim(basename)) // sdata(2:10) // ".vorme"
!      call singlefile_io3(time,work1,fname,vor,work2,0,io_pe,.false.,3)

   endif
   if (convert_opt==3) then  ! -cout norm2
      ! 2048^3 needs 192GB * 1.66.  needs 256 cpus
      call input_uvw(time,Q,vor,work1,work2,1)
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
   if (convert_opt==5) then  ! -cout norm
      ! 2048^3 needs 192GB * 1.66.  needs 256 cpus
      call input_uvw(time,Q,vor,work1,work2,1)
      print *,'max input: ',maxval(Q(nx1:nx2,ny1:ny2,nz1:nz2,1)), &
                            maxval(Q(nx1:nx2,ny1:ny2,nz1:nz2,2)), &
                            maxval(Q(nx1:nx2,ny1:ny2,nz1:nz2,3))
      call print_message("computing norm squared...")
      do k=nz1,nz2
      do j=ny1,ny2
      do i=nx1,nx2
         work1(i,j,k)=sqrt(Q(i,j,k,1)**2+Q(i,j,k,2)**2+Q(i,j,k,3)**2)
      enddo
      enddo
      enddo
      call print_message("outputting norm squared as REAL*4...")
      output_size=4
      write(sdata,'(f10.4)') 10000.0000 + time
      basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))
      fname = basename(1:len_trim(basename)) // sdata(2:10) // ".norms"
      call singlefile_io3(time,work1,fname,Q,work2,0,io_pe,.false.,2)
      fname = basename(1:len_trim(basename)) // sdata(2:10) // ".normse"
      call singlefile_io3(time,work1,fname,Q,work2,0,io_pe,.false.,3)
   endif


   if (convert_opt==6) then  ! -cout passive
      schmidt_in=1.0
      type_in=0
      write(message,'(f10.4)') 10000.0000 + time
      write(ext,'(f8.3)') 1000 + schmidt_in
      write(ext2,'(i3)') 100+type_in
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
           // message(2:10) // '.t' // ext2(2:3) // '.s' // ext(2:8)
      call print_message(rundir)
      call print_message(runname)
      call print_message(fname)	
      call singlefile_io3(time,Q,fname,work1,work2,1,io_pe,.false.,1)
      call global_min(Q,mn)
      call global_max(Q,mx)
      write(message,'(a,2f17.5)') 'passive scalar min/max: ',mn,mx
      call print_message(message)	

      write(message,'(f10.4)') 10000.0000 + time
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
           // message(2:10) // '.t' // ext2(2:3) // '.s' // ext(2:8) &
          // '-raw'
      call print_message(fname)
      call singlefile_io3(time,Q,fname,work1,work2,0,io_pe,.false.,2)
   endif
   if (convert_opt==7) then  ! -cout gradu
      call input_uvw(time,Q,vor,work1,work2,1)
      call print_message("computing gradu ")
      call gradu_stats(time,Q,vor,work1,work2)
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







subroutine gradu_stats(time,Q,vor,work1,work2)
use params
use mpi
use fft_interface
use subcubes
implicit none
real*8 :: time
real*8 :: Q(nx,ny,nz,3)
real*8 :: vor(nx,ny,nz)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
integer i,j
character(len=280) basename,fname,sdata
character(len=80) :: message
character(len=4) :: extension="uvwX"
character(len=8) :: ext2,ext

integer :: sc
integer,parameter :: ssize=128
real*8 :: dx(1:2*ssize),dy(1:2*ssize),dz(1:2*ssize)
real*8 :: dslice(2*ssize,2*ssize,1)
real*8 :: mat(3,3,3)
real*8,allocatable :: data(:,:,:)
real*8,allocatable :: gradu(:,:,:)
real*8,allocatable :: gradu2(:,:,:)
integer :: i2,j2,k2,k

call setup_subcubes(ssize)
write(message,'(a,i5)') 'number of subcubes: ',nsubcube
call print_message(message)


allocate(gradu(3,3,nsubcube))
allocate(gradu2(3,3,nsubcube))
allocate(data(ssize,ssize,ssize))
gradu=0
gradu2=0


do i=1,3
do j=1,3    
   ! u_i,j
   write(message,'(a,i1,a,i1,a)') 'computing u_',i,',',j,' ...'
   call print_message(message)
   call der(Q(1,1,1,i),vor,work1,work2,DX_ONLY,j)
   do sc=1,nsubcube

      dx(1:ssize)=subcube_corner(1,sc)+subcube_cords(:)
      dy(1:ssize)=subcube_corner(2,sc)+subcube_cords(:)
      dz(1:ssize)=subcube_corner(3,sc)+subcube_cords(:)
      call interp_subcube(ssize,ssize,ssize,dx,dy,dz,data,vor,0,mat)

      if (io_pe==my_pe) then
         ! data() only computed on io_pe
         do k2=1,ssize
         do j2=1,ssize
         do i2=1,ssize
            gradu(i,j,sc)=gradu(i,j,sc)+data(i2,j2,k2)
            gradu2(i,j,sc)=gradu2(i,j,sc)+data(i2,j2,k2)**2
         enddo
         enddo
         enddo
         gradu(i,j,sc)=gradu(i,j,sc)/ssize/ssize/ssize
         gradu2(i,j,sc)=gradu2(i,j,sc)/ssize/ssize/ssize
      endif
   enddo

#if 0   
   write(sdata,'(f10.4)') 10000.0000 + time
   write(ext,'(2i1)') i,j
   basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))
   fname = basename(1:len_trim(basename)) // sdata(2:10) // '.' // ext(1:2)
   call singlefile_io3(time,vor,fname,work1,work2,0,io_pe,.false.,2)
   
   if (j==1) then
      basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
           // "-raw"
      fname = basename(1:len_trim(basename)) // sdata(2:10) // "." &
           // extension(i:i)
      call singlefile_io3(time,vor,fname,work1,work2,0,io_pe,.false.,2)
   endif
#endif
enddo
enddo

! output the gradu matricies:
if (io_pe==my_pe) then
write(sdata,'(f10.4)') 10000.0000 + time
write(ext,'(2i1)') i,j
basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))
fname = basename(1:len_trim(basename)) // sdata(2:10) // '.gradu'
open(15,file=fname,form='formatted')
do sc=1,nsubcube
   write(15,'(3f12.8,i5)') subcube_corner(1:3,sc),ssize
   do i=1,3
      write(15,'(3e18.10)') (gradu(i,j,sc),j=1,3)
   enddo
enddo
close(15)



write(sdata,'(f10.4)') 10000.0000 + time
write(ext,'(2i1)') i,j
basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))
fname = basename(1:len_trim(basename)) // sdata(2:10) // '.gradu2'
open(15,file=fname,form='formatted')
do sc=1,nsubcube
   write(15,'(3f12.8,i5)') subcube_corner(1:3,sc),ssize
   do i=1,3
      write(15,'(3e18.10)') (gradu2(i,j,sc),j=1,3)
   enddo
enddo
close(15)
endif


deallocate(data)
deallocate(gradu)
deallocate(gradu2)

return

call ghost_update(vor,1)

! loop thru gradu matrix, looking for selections
do sc=1,nsubcube
   write(sdata,'(f10.4)') 10000.0000 + time
   write(ext,'(i5)') sc  ! 10000
   basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))
   fname = basename(1:len_trim(basename)) // sdata(2:10) // '.sc' // ext(2:5)


   ! select a subcube 2x larger containg the sc'th subcube:
   dx(:)=wsubcube_corner(1,sc)+wsubcube_cords(:)
   dy(:)=wsubcube_corner(2,sc)+wsubcube_cords(:)
   dz(:)=wsubcube_corner(3,sc)+wsubcube_cords(:)
   do k=1,2*ssize
      call interp_subcube(2*ssize,2*ssize,1,dx,dy,dz(k),dslice,vor,0,mat)
      ! output this slab:
   enddo
enddo


end subroutine
