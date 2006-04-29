#include "macros.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read in a sequence of data files and convert them
!
! -zi       input file uses NS_UVW compression (only works with input_uvw)
! -zo       ouput data using NS_UVW compression (only works with ouput_uvw)
!
! -si       input file contains spectral coefficients instead of grid data
! -so       output spectral coefficieints instead of grid data
!
! -o4       output_size=4 (default is 8 bytes)
! -i4       input_size=4 (instead of default of 8) 
!
! type of conversion controled by -cout argument:
!
!  -cout uvw      (can be used with -si/-so  spectral input/output)
!                  and -smax spectral coefficient truncation to 
!                  perform downsampling or upsampling)
!  -cout 4uvw     as above, but loops over fields 1 by 1 
!                   (and only needs 3 scalar arrays, but cant do stats,compressed output) 
!                   use -o4 to get real*4 output
!  -cout vor
!  -cout vorm
!  -cout norm2
!  -cout passive    convert passive scalar file
!                   need to specify shmidt_in and type_in below
!  -cout gradu      output <u_i,j>  matrixes for subcubes
!  -cout stats      read data, print some stats
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
use spectrum
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
tstart=.0000
tstop=5.0
tinc=.5

! to read times from  file times.dat:
  tstart=-1; tinc=0; tname="times.dat"

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

if (convert_opt==0 .or. convert_opt == 3 .or. convert_opt==5 .or. &
    convert_opt==10 ) then
   allocate(vor(nx,ny,nz,n_var)) ! used for shallow water output routiens
   allocate(Q(nx,ny,nz,n_var))
else if (convert_opt == 4 .or. convert_opt==6) then
   allocate(vor(1,1,1,1)) ! dummy variable -wont be used
   allocate(Q(nx,ny,nz,1)) ! only need 1 slot
else if (convert_opt == 7 .or. convert_opt==8) then
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
   write(message,'(a,i4,a,f10.4)') 'iter=',icount,' attempting to read time=',time
   call print_message(message)

   Q=0

   if (convert_opt==0) then  ! -cout uvw  
      ! read data, header type =1, or specified in input file
      time2=time
      call input_uvw(time2,Q,vor,work1,work2,header_user)  
      if (g_nz < 2048) call print_stats(Q,vor,work1,work2)

      ! just reoutput the variables:
      if (w_spec) then
         do n=1,3
            write(message,'(a,i4)') 'w_spec fft3d: n=',n
            call print_message(message)
            call fft3d(Q(1,1,1,n),work1)
         enddo
      endif
      if ((w_spec /= r_spec) .or. (w_compressed /= r_compressed)) then
         ! converting from spec to grid (or vice versa) 
         basename=runname(1:len_trim(runname))
      else
         ! dont clobber input file!
         basename=runname(1:len_trim(runname)) // "-new."
      endif
      call output_uvw(basename,time2,Q,vor,work1,work2,header_user)  
      ! output headerless data:
      ! call output_uvw(basename,time,Q,vor,work1,work2,2)

   endif

   if (convert_opt==1) then  ! -cout vor
      call input_uvw(time,Q,vor,work1,work2,header_user)
      ! outputing vorticity
      basename=runname(1:len_trim(runname)) // "-vor."

      call print_message("computing vorticity...")
      call vorticity(vor,Q,work1,work2)
      call print_message("output vorticity...")
      call output_uvw(basename,time,vor,Q,work1,work2,header_user)
   endif

   if (convert_opt==2) then  ! -cout vorm
      call input_uvw(time,Q,vor,work1,work2,header_user)
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
      call input_uvw(time,Q,vor,work1,work2,header_user)
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
      ! this could will output 1 field at a time.  only needs 3 arrays.
      write(sdata,'(f10.4)') 10000.0000 + time
      basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))

      do i=1,3  ! use this to loop over .u,.v,.w
         fname = basename(1:len_trim(basename)) // sdata(2:10) // &
                 "." // extension(i:i)
         call print_message("input file:")	
	 call print_message(fname(1:len_trim(fname)))
         call singlefile_io3(time,Q,fname,work1,work2,1,io_pe,r_spec,header_user)
         print *,'max input: ',maxval(Q(nx1:nx2,ny1:ny2,nz1:nz2,1))

         if (w_spec) then
            do n=1,3
               call fft3d(Q(1,1,1,n),work1)
            enddo
            fname = basename(1:len_trim(basename)) // sdata(2:10) // &
                 "." // extension(i:i) // "s"
         else
            fname = basename(1:len_trim(basename)) // sdata(2:10) // &
                 "." // extension(i:i) // "-new"
         endif
         call print_message(fname(1:len_trim(fname)))

         call singlefile_io3(time,Q,fname,work1,work2,0,io_pe,w_spec,header_user)
         ! headerless:
         ! call singlefile_io3(time,Q,fname,work1,work2,0,io_pe,w_spec,2)
      enddo	
   endif
   if (convert_opt==5) then  ! -cout norm
      ! 2048^3 needs 192GB * 1.66.  needs 256 cpus
      call input_uvw(time,Q,vor,work1,work2,header_user)
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
      ! compute <gradu> for all subcubes, output in asci file
      call input_uvw(time,Q,vor,work1,work2,header_user)
      call print_message("computing gradu ")
      call gradu_stats(time,Q,vor,work1,work2)
   endif
   if (convert_opt==8) then  ! -cout extract_subcube
      ! read asci file for list of subcubes to extractg
      ! rotate (if necessary) to align gradient with x axis
      ! output raw brick-of-floats for each rotated subcube
      call input_uvw(time,Q,vor,work1,work2,header_user)
      call print_message("extracting subcubes ")
      call gradu_rotate_subcubes(time,Q,vor,work1,work2)
   endif
   if (convert_opt==9) then  ! -cout spec_window  
      ! read input data, detrend, window, output spectrum and cospectrum
      !call input_uvw(time,Q,vor,work1,work2,header_user) ! DNS default headers
      call input_uvw(time,Q,vor,work1,work2,2) ! no headers
      call print_message("computing spectrum ")

      ! compute U using random vorticity: (for testing)
      ! call ranvor(Q,vor,work1,work2,1)

      do n=1,3
         call detrend_data(Q(1,1,1,n),0)
      enddo
      call compute_spec(time,Q,vor,work1,work2)
      call output_spec(time,time)
      call output_helicity_spec(time,time) 
   endif

   if (convert_opt==10) then ! -cout iotest
      ! i/o performance
      write(message,*) 'computing random data'
      call print_message(message)
      do k=nz1,nz2
         do j=ny1,ny2
            call gaussian( Q(nx1,j,k,1), nslabx  )
         enddo
      enddo
      call print_message('performing io test')
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
           // ".u"

      ! stripe, num_cpus, 

      if (ncpu_z<16) then  ! only run these tests when debugging:
      if (ncpu_z>=4) then
         call iotest(1,2,fname,Q,work1,work2)
         call iotest(2,2,fname,Q,work1,work2)
         call iotest(4,2,fname,Q,work1,work2)
      endif
      if (ncpu_z>=8) then
         call iotest(1,4,fname,Q,work1,work2)
         call iotest(2,4,fname,Q,work1,work2)
         call iotest(4,4,fname,Q,work1,work2)
      endif
      endif

      if (ncpu_z>=16) then
         call iotest(2,8,fname,Q,work1,work2)
         call iotest(4,8,fname,Q,work1,work2)
         call iotest(8,8,fname,Q,work1,work2)
      endif
      if (ncpu_z>=32) then
         call iotest(4,16,fname,Q,work1,work2)
         call iotest(8,16,fname,Q,work1,work2)
         call iotest(16,16,fname,Q,work1,work2)
      endif
      if (ncpu_z>=64) then
         call iotest(8,32,fname,Q,work1,work2)
         call iotest(16,32,fname,Q,work1,work2)
         call iotest(32,32,fname,Q,work1,work2)
      endif
      if (ncpu_z>=128) then
         call iotest(16,64,fname,Q,work1,work2)
         call iotest(32,64,fname,Q,work1,work2)
         call iotest(64,64,fname,Q,work1,work2)
      endif
      if (ncpu_z>=256) then
         call iotest(16,128,fname,Q,work1,work2)
         call iotest(32,128,fname,Q,work1,work2)
         call iotest(64,128,fname,Q,work1,work2)
      endif


   endif


   if (convert_opt==11) then ! -cout stats    test compression data
      time2=time
      call input_uvw(time2,Q,vor,work1,work2,header_user)  
      call print_stats(Q,vor,work1,work2)
   endif


   if (tstart>=0) then
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




subroutine iotest(stripe,num_io_cpu,fname,p,work1,work2)
use transpose
implicit none
real*8 :: p(nx,ny,nz)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
character(len=*) :: fname
integer :: num_io_cpu,stripe,ierr
real*8 :: tmx1,tmx2,time=1111.1111    

#ifdef USE_MPI
mpi_maxio=num_io_cpu
write(mpi_stripe,'(i3)') stripe
call mpi_io_init(0)

call mpi_barrier(comm_3d,ierr)
call wallclock(tmx1)
call singlefile_io3(time,p,fname,work1,work2,0,io_pe,.false.,2)
call mpi_barrier(comm_3d,ierr)
call wallclock(tmx2)
tmx2=tmx2-tmx1

if (io_pe==my_pe) then
   print *,'stripe=',mpi_stripe
   write(*,'(a,f12.5,a,f12.5,a)') 'cpu time for output: ',tmx2,'s (',tmx2/60,'m)'
   print *,'data rate MB/s: ',8.*3.*g_nx*g_ny*g_nz/1024./1024./tmx2
endif
#endif
end subroutine iotest






subroutine detrend_data(p,wtype)
use params
use mpi
use ghost
implicit none
real*8 :: p(nx,ny,nz)
integer :: wtype
integer :: i,j,k,ierr
real*8 :: gradp(3),meanp,gradp2(3),px,py,pz
real*8 :: xwindow(nx)
real*8 :: ywindow(ny)
real*8 :: zwindow(nz)
integer :: xinside(nx),yinside(ny),zinside(nz),grad_count

if (nx1<3 .or. nx2+2>nx)&
     call abort('Error: insufficient ghost cells in x direction')
if (ny1<3 .or. ny2+2>ny)&
     call abort('Error: insufficient ghost cells in y direction')
if (nz1<3 .or. nz2+2>nz) then
     call abort('Error: insufficient ghost cells in z direction')
endif

! update ghost cell data
call ghost_update_x(p,1)
call ghost_update_y(p,1)
call ghost_update_z(p,1)


! compute window function
if (wtype==0) then
   ! domain: 0..1   
   ! 0-0.25    .5-.5*cos(x*pi/.25)
   ! .25-.75   1
   ! .75-1.0   .5+.5*cos((x-.75)*pi/.25)
   !
   xinside=1
   yinside=1
   zinside=1
   xwindow=1
   ywindow=1
   zwindow=1
   do i=nx1,nx2
      if (xcord(i)<=.25) then
         xwindow(i)=.5-.5*cos(xcord(i)*pi*4)
         xinside(i)=0
      endif
      if (xcord(i)>=.75) then
         xwindow(i)=.5+.5*cos((xcord(i)-.75)*pi*4)
         xinside(i)=0
      endif
   enddo
   do j=ny1,ny2
      if (ycord(j)<=.25) then
         ywindow(j)=.5-.5*cos(ycord(j)*pi*4)
         yinside(j)=0
      endif
      if (ycord(j)>=.75) then
         ywindow(j)=.5+.5*cos((ycord(j)-.75)*pi*4)
         yinside(j)=0
      endif
   enddo
   do k=nz1,nz2
      if (zcord(k)<=.25) then
         zwindow(k)=.5-.5*cos(zcord(k)*pi*4)
         zinside(k)=0
      endif
      if (zcord(k)>=.75) then
         zwindow(k)=.5+.5*cos((zcord(k)-.75)*pi*4)
         zinside(k)=0
      endif
   enddo

endif

#ifdef DETREND_DEBUG
if (io_pe==my_pe) then
!print *,'xwindow: ',xwindow
!print *,'ywindow: ',ywindow
!print *,'zwindow: ',zwindow
endif
#endif


gradp=0
grad_count=0
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   if (xinside(i)==1 .and. yinside(j)==1 .and. zinside(k)==1) then 
      px=(2*(p(i+1,j,k)-p(i-1,j,k))/3 - (p(i+2,j,k)-p(i-2,j,k))/12)/delx
      py=(2*(p(i,j+1,k)-p(i,j-1,k))/3 - (p(i,j+2,k)-p(i,j-2,k))/12)/dely
      pz=(2*(p(i,j,k+1)-p(i,j,k-1))/3 - (p(i,j,k+2)-p(i,j,k-2))/12)/delz
      grad_count=grad_count+1
      gradp(1)=gradp(1)+px
      gradp(2)=gradp(2)+py
      gradp(3)=gradp(3)+pz
   endif
enddo
enddo
enddo

#ifdef USE_MPI
gradp2=gradp
call mpi_allreduce(gradp2,gradp,3,MPI_REAL8,MPI_SUM,comm_3d,ierr)
i=grad_count
call mpi_allreduce(i,grad_count,1,MPI_INTEGER,MPI_SUM,comm_3d,ierr)
#endif
gradp=gradp/grad_count


! remove gradiant
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   px = gradp(1)*(xcord(i)-.5) + gradp(2)*(ycord(j)-.5) + gradp(3)*(zcord(k)-.5)
   p(i,j,k)=p(i,j,k) -  px 
enddo
enddo
enddo


meanp=0
grad_count=0
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   meanp=meanp+p(i,j,k)
enddo
enddo
enddo
#ifdef USE_MPI
px=meanp
call mpi_allreduce(px,meanp,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
meanp=meanp/g_nx/g_ny/g_nz

! remove mean
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   p(i,j,k)=p(i,j,k) - meanp
enddo
enddo
enddo


if (io_pe==my_pe) then
!   print *,'mean: ',meanp
   write(*,'(a,3f15.5)') 'grap*4         ',4*gradp(1:3)
endif




! compute means again and see if we really detrended the data:
#define DETREND_DEBUG
#ifdef DETREND_DEBUG
meanp=0
gradp=0
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   meanp=meanp+p(i,j,k)
   px=(2*(p(i+1,j,k)-p(i-1,j,k))/3 - (p(i+2,j,k)-p(i-2,j,k))/12)/delx
   py=(2*(p(i,j+1,k)-p(i,j-1,k))/3 - (p(i,j+2,k)-p(i,j-2,k))/12)/dely
   pz=(2*(p(i,j,k+1)-p(i,j,k-1))/3 - (p(i,j,k+2)-p(i,j,k-2))/12)/delz
   gradp(1)=gradp(1)+px
   gradp(2)=gradp(2)+py
   gradp(3)=gradp(3)+pz
enddo
enddo
enddo

#ifdef USE_MPI
px=meanp
call mpi_allreduce(px,meanp,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
gradp2=gradp
call mpi_allreduce(gradp2,gradp,3,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif

meanp=meanp/g_nx/g_ny/g_nz
gradp=gradp/g_nx/g_ny/g_nz
if (io_pe==my_pe) then
!   print *,'detrended mean :   ',meanp
!   write(*,'(a,3f15.5)') 'detrended 4*gradp:   ',4*gradp(1:3)
endif
#endif


!apply window
do k=nz1,nz2
do j=ny2,ny2
do i=nx1,nz2
   p(i,j,k)=p(i,j,k)*xwindow(i)*ywindow(j)*zwindow(k)
enddo
enddo
enddo


end subroutine


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
real*8,allocatable :: data(:,:,:)
real*8,allocatable :: gradu(:,:,:,:)
integer :: i2,j2,k2,k,ip

call setup_subcubes(ssize)
write(message,'(a,i5)') 'number of subcubes: ',nsubcube
call print_message(message)


allocate(gradu(3,3,nsubcube,2))
allocate(data(ssize,ssize,ssize))
gradu=0


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
      call extract_subcube(ssize,ssize,ssize,dx,dy,dz,data,vor)

      if (io_pe==my_pe) then
         ! data() only computed on io_pe
         do k2=1,ssize
         do j2=1,ssize
         do i2=1,ssize
            gradu(i,j,sc,1)=gradu(i,j,sc,1)+data(i2,j2,k2)
            gradu(i,j,sc,2)=gradu(i,j,sc,2)+data(i2,j2,k2)**2
         enddo
         enddo
         enddo
         gradu(i,j,sc,:)=gradu(i,j,sc,:)/ssize/ssize/ssize
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
   do ip=1,2
      write(sdata,'(f10.4)') 10000.0000 + time
      write(ext,'(i1)') ip
      basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))
      fname = basename(1:len_trim(basename)) // sdata(2:10) // '.gradu' //  ext(1:1)
      open(15,file=fname,form='formatted')
      do sc=1,nsubcube
         write(15,'(3f12.8,i5)') subcube_corner(1:3,sc),ssize
         do i=1,3
            write(15,'(3e18.10)') (gradu(i,j,sc,ip),j=1,3)
         enddo
      enddo
      close(15)
   enddo
endif


deallocate(data)
deallocate(gradu)

end subroutine







subroutine gradu_rotate_subcubes(time,Q,vor,work1,work2)
use params
use mpi
use fft_interface
use subcubes
use ghost
implicit none
real*8 :: time
real*8 :: Q(nx,ny,nz,3)
real*8 :: vor(nx,ny,nz)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
integer :: i,j,k,ierr
character(len=280) basename,fname,sdata
character(len=80) :: message
character(len=4) :: extension="uvwX"
character(len=8) :: ext2,ext

integer :: sc
integer :: ssize,irot,il,jl
real*8,allocatable :: dslice(:,:,:)
real*8 :: mat(3,3,3),gradu(3,3)
real*8 :: corner(3),x0(3),x1(3),x2(3),x3(3),e01(3),e02(3),e03(3),S12,c0(3)
CPOINTER :: fid

! check to make sure we have 2 ghost cells in all directions:
if (nx1<3 .or. nx2+2>nx)&
     call abort('Error: insufficient ghost cells in x direction')
if (ny1<3 .or. ny2+2>ny)&
     call abort('Error: insufficient ghost cells in y direction')
if (nz1<3 .or. nz2+2>nz)&
     call abort('Error: insufficient ghost cells in z direction')

! update ghost cell data
call ghost_update_x(Q,3)
call ghost_update_y(Q,3)
call ghost_update_z(Q,3)


write(sdata,'(f10.4)') 10000.0000 + time
basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))
fname = basename(1:len_trim(basename)) // sdata(2:10) // '.select'
call print_message(fname)
open(84,err=100,file=fname(1:len_trim(fname)),form='formatted')

sc=0  ! subcube counter
do
   sc=sc+1
   ! read coordinates of subcube corner and size (gridpoints)
   ! corner(1:3), ssize, gradu
   read(84,*,err=100,end=100) corner(1),corner(2),corner(3),ssize
   do i=1,3
      read(84,*,err=100,end=100) (gradu(i,j),j=1,3)
   enddo

   if (io_pe==my_pe) then
      write(*,'(a,3f12.5)') 'corner: ',corner(1:3)
   endif

   ! compute corner of larger cube containing subcube (x0(:))
   x0(1) = corner(1) - ssize*delx/2
   x0(2) = corner(2) - ssize*dely/2
   x0(3) = corner(3) - ssize*delz/2

   ! compute 3 corner points connected to x0(:)
   x1 = x0; x1(1)=x0(1) + 2*ssize*delx
   x2 = x0; x2(2)=x0(2) + 2*ssize*dely
   x3 = x0; x3(3)=x0(3) + 2*ssize*delz

   ! compute center of cube:
   c0(1)=x0(1) + ssize*delx
   c0(2)=x0(2) + ssize*dely
   c0(3)=x0(3) + ssize*delz


   ! compute rotation matrix from gradu(3,3):
   ! apply rotation to x0,x1,x2,x3
  
   ! find largest element of gradu()
   S12=gradu(1,2)
   do i=1,3
   do j=1,3
      if (abs(gradu(i,j))>abs(S12)) then
         il=i
         jl=j
         S12=gradu(i,j)
      endif
   enddo
   enddo



   ! rotate, if needed
   if (il==1 .and. jl==2) then
      ! no rotation needed
   else if (il==1 .and. jl==3) then
      ! swap 2 & 3
      call rotcube(2,3,x0,x1,x2,x3,c0)
   else if (il==2 .and. jl==3) then
      ! swap 1 & 2, swap 3 & 2
      call rotcube(1,2,x0,x1,x2,x3,c0)
      call rotcube(2,3,x0,x1,x2,x3,c0)
   else if (il==2 .and. jl==1) then
      ! swap 1 & 2
      call rotcube(1,2,x0,x1,x2,x3,c0)
   else if (il==3 .and. jl==1) then
      ! swap 1 & 2, swap 3 & 1
      call rotcube(1,2,x0,x1,x2,x3,c0)
      call rotcube(1,3,x0,x1,x2,x3,c0)
   else if (il==3 .and. jl==2) then
      ! swap 3 & 1
      call rotcube(1,3,x0,x1,x2,x3,c0)
   else
      call abort("error: largest <gradu> element is on diagional")
   endif


   irot=0  ! irot=0: uses closest gridpoint.  for 0 or 90 degree rotations
           ! irot=1: 4th order interpolation 

   ! compute unit vectors
   ! x0 + unit vectors will be used to generate grid of points
   ! where we need to interpolate
   e01=(x1-x0)/(2*ssize*delx)
   e02=(x2-x0)/(2*ssize*dely)
   e03=(x3-x0)/(2*ssize*delz)


   ! allocate memory:
   if (sc==1) then
      allocate(dslice(2*ssize,2*ssize,1) )
   endif


   do j=1,3
      ! compute filename.  use runname+"-sc#_"+time+".uvw"
      write(sdata,'(f10.4)') 10000.0000 + time
      write(ext,'(i5)') 10000+sc  ! 10000
      basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))
      fname = basename(1:len_trim(basename)) // "-sc" // ext(2:5) // &
              "_" // sdata(2:10) // '.' // extension(j:j)
      if (io_pe==my_pe) &
          call copen(fname,"w",fid,ierr)
      
      do k=1,2*ssize ! loop over all z-slices
         ! interpolate 1 slice of data
         call interp_subcube(2*ssize,k,k,x0,e01,e02,e03,dslice,Q(1,1,1,j),irot)
         ! output this slab:
         if (io_pe==my_pe) &
             call cwrite8(fid,dslice,2*ssize*2*ssize)
      enddo
      
      if (io_pe==my_pe) &
          call cclose(fid,ierr)
   enddo

    

enddo
return

100 continue
call print_message("Error reading subcube list file...")
end subroutine



subroutine print_stats(Q,div,work1,work2)
use params
use mpi
implicit none
real*8 :: Q(nx,ny,nz,3)
real*8 :: div(nx,ny,nz,3)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: mx(3),mx2(3),divx,divi,ens,ke
integer :: n,ierr,i,j,k
character(len=280) :: message

call vorticity(div,Q,work1,work2)
ens=0
ke=0
do n=1,3
do k=nz1,nz2
   do j=ny1,ny2
      do i=nx1,nx2
         ens=ens + div(i,j,k,n)**2
         ke=ke+.5*Q(i,j,k,n)**2
      enddo
   enddo
enddo
enddo
ens=ens/g_nx/g_ny/g_nz
ke=ke/g_nx/g_ny/g_nz
#ifdef USE_MPI
   divi=ke
   call mpi_allreduce(divi,ke,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   divi=ens
   call mpi_allreduce(divi,ens,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif



write(message,'(a,3f18.14)') 'STATS:  KE = ',ke
call print_message(message)


do n=1,3
   mx(n)=maxval(abs(Q(nx1:nx2,ny1:ny2,nz1:nz2,n)))
enddo 
#ifdef USE_MPI
   mx2=mx
   call mpi_allreduce(mx2,mx,3,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif

write(message,'(a,3f18.14)') 'STATS:  maxU = ',mx
call print_message(message)

call compute_div(Q,div,work1,work2,divx,divi)
write(message,'(3(a,e12.5))') 'STATS:  max(div)=',divx
call print_message(message)	


write(message,'(3(a,f18.12))') 'STATS:  enstrophy=',ens
call print_message(message)	




end subroutine
