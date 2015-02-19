!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Copyright 2007.  Los Alamos National Security, LLC. This material was
!produced under U.S. Government contract DE-AC52-06NA25396 for Los
!Alamos National Laboratory (LANL), which is operated by Los Alamos
!National Security, LLC for the U.S. Department of Energy. The
!U.S. Government has rights to use, reproduce, and distribute this
!software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
!LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
!FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
!derivative works, such modified software should be clearly marked, so
!as not to confuse it with the version available from LANL.
!
!Additionally, this program is free software; you can redistribute it
!and/or modify it under the terms of the GNU General Public License as
!published by the Free Software Foundation; either version 2 of the
!License, or (at your option) any later version. Accordingly, this
!program is distributed in the hope that it will be useful, but WITHOUT
!ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
!for more details.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!  -cout 4uvw (4)    as above, but loops over fields 1 by 1 
!                   (and only needs 3 scalar arrays, but cant do stats,compressed output) 
!                   use -o4 to get real*4 output
!  -cout vor  (1)
!  -cout vorm (2)
!  -cout norm  (5)
!  -cout norm2 (3)
!  -cout passive (6) convert passive scalar file
!                    need to specify shmidt_in and type_in below
!  -cout gradu (7)  output <u_i,j>  matrixes for subcubes
!  -cout extract_subcube (8)  ?
!  -cout spec_window  (9)
!  -cout iotest(10) ?
!  -cout stats (11) read data, print some stats
!  -cout uwbar (12) ?
!  -cout trunc (13) read data, fft, truncate, fft back. (Only does uvw truncation)

!
! option used to check aliasing error:   
!  -cout nlout (14) compute random U,V,UV.  output U,V,UV
!  -cout nlin  (15) read in U,V,UV.  recompute UV (on a higher res grid)
!                   compare UV2,UV in spectral space
!
!  -cout dfilter  (16)   read in U,V,W, output delta-filtered U
!  -cout dpdf     (17)   read in U,V,W, output delta-filtered U pdfs
!  -cout hpass (18) read data, fft, high-pass filter, fft back.
!
!  -cout coarsegrain (19)  read data, coarse grain, output
!  -cout dudx  (20)
!  -cout ehor  (21)    horizontal energy field
!  -cout potens (22)   potential enstrophy field
!  -cout pe (23) potential energy field
!  -cout CH_decomp (24) write wave modes and vortical modes field
!  -cout ttrunc (25) read scalar, fft, truncate fft back
!  -cout zerocr (26) compute zero-crossing freqency PDFs for u,v,w in x,yand z

! To run, set the base name of the file and the times of interest
! below.  For example:
!    tstart=4.000000000 edit below
!    tstop=4.5
!    tinc=.5
!    basename="temp"
!
! will result in looping over the files:   
!             temp0000.0000.[uvw]
!             temp0000.5000.[uvw]
!             temp0001.0000.[uvw]
!
!  to compile and run:   make convert ; convert
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program convert
use params
use mpi
use fft_interface
use spectrum
use pdf
use transpose , only : input1
implicit none
real*8,allocatable  :: Q(:,:,:,:)
real*8,allocatable  :: Q2(:,:,:,:)
real*8,allocatable  :: vor(:,:,:,:)
real*8,save  :: work1(nx,ny,nz)
real*8,save  :: work2(nx,ny,nz)
character(len=80) message,sdata,tname,sdata2
character(len=280) basename,fname
integer ierr,i,j,k,n,km,im,jm,icount
real*8 :: tstart,tstop,tinc,time,time2,time_tmp
real*8 :: u,v,w,x,y
real*8 :: kr,ke,ck,xfac,dummy,xtmp
real*8 :: schmidt_in,mn,mx,a0,a1
real*8 :: xwerr(1000),xw,binsize
CPOINTER :: null=0,fidu,fid
integer :: type_in,ntot,nzero,nerr, kshell_max,k2
character(len=4) :: extension="uvwX"
character(len=8) :: ext2,ext
integer :: Mval(4) = (/4,16,32,64/)   ! values used for coarse graining


! input file
tstart=0.0387
tstop=0.0387
tinc=0.0001


! to read times from  file times.dat:
!tstart=-1; tinc=0; tname="times.dat"

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

if (convert_opt == 4 .or. convert_opt==6) then
   ! special option to save storage, for processing very large data
   allocate(vor(1,1,1,1)) ! dummy variable -wont be used
   allocate(Q(nx,ny,nz,1)) ! only need 1 slot
else if (convert_opt == 7 .or. convert_opt==8 ) then
   ! special option to save storage, for processing very large data
   allocate(vor(nx,ny,nz,1)) ! only first component used
   allocate(Q(nx,ny,nz,n_var))
else
   ! default: allocate 2 3D arrays
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
      time_tmp=time
      call input_uvw(time_tmp,Q,vor,work1,work2,header_user)  
      if (.not. r_spec) then  ! r_spec reader will print stats, so we can skip this:
         call print_stats(Q,vor,work1,work2)
      endif
      
      ! just reoutput the variables:
      if (w_spec) then
         do n=1,3
            write(message,'(a,i4)') 'w_spec fft3d: n=',n
            call print_message(message)
            call fft3d(Q(1,1,1,n),work1)
         enddo
      endif
      if ((w_spec .neqv. r_spec) .or. (w_compressed .neqv. r_compressed)) then
         ! converting from spec to grid (or vice versa) 
         basename=runname(1:len_trim(runname))
      else
         ! dont clobber input file!
         basename=runname(1:len_trim(runname)) // "raw."
      endif
!      call output_uvw(basename,time2,Q,vor,work1,work2,header_user)  
      ! output headerless data:
       call output_uvw(basename,time,Q,vor,work1,work2,2)
       
    endif

    if (convert_opt==1) then  ! -cout vor
       time_tmp=time
       call input_uvw(time_tmp,Q,vor,work1,work2,header_user)
       ! outputing vorticity
       if (ndim==3) then
          basename=runname(1:len_trim(runname)) // "-vor."
          call print_message("computing vorticity...")
          call vorticity(vor,Q,work1,work2)
          call print_message("output vorticity...")
          ! call output_uvw(basename,time,vor,Q,work1,work2,header_user)
          ! output headerless
          call output_uvw(basename,time,vor,Q,work1,work2,2)
          
       endif
       if (ndim==2) then
          call print_message("computing 2D vorticity...")
          call vorticity2d(vor,Q,work1,work2)
          write(sdata,'(f10.4)') 10000.0000 + time
          fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) // ".vor"
          call singlefile_io3(time,vor,fname,work1,work2,0,io_pe,.false.,header_user)
          
          ! compute the v vorticity, in vor(:,:,:,2)
          if (alpha_value>0) then
             call print_message("computing 2D v vorticity...")
             vor(:,:,:,2)=vor(:,:,:,1)
             call v_vorticity2d(vor(1,1,1,2),work1)
             fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) // ".vvor"
             call singlefile_io3(time,vor(1,1,1,2),fname,work1,work2,0,io_pe,.false.,header_user)
          endif
          
          ! compute the stream function, stored in vor(:,:,:,2)
          call print_message("computing 2D stream function...")
          vor(:,:,:,2)=vor(:,:,:,1)
          call fft3d(vor(1,1,1,2),work1)
          a0=0; a1=-1
          call fft_laplace_inverse(vor(1,1,1,2),a0,a1)
          call ifft3d(vor(1,1,1,2),work1)
          
          fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) // ".psi"
          call singlefile_io3(time,vor(1,1,1,2),fname,work1,work2,0,io_pe,.false.,header_user)
       endif

       
       
    endif
    
    if (convert_opt==2) then  ! -cout vorm
       time_tmp=time
       call input_uvw(time_tmp,Q,vor,work1,work2,header_user)
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
       time_tmp=time
       call input_uvw(time_tmp,Q,vor,work1,work2,header_user)
       do i=1,ndim
          call print_max(Q(1,1,1,i),i)
       enddo
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
          call print_max(Q,i)
          
          if (w_spec) then
             write(message,'(a,i4)') 'w_spec fft3d: n=',i
             call print_message(message)
             call fft3d(Q(1,1,1,1),work1)
             fname = basename(1:len_trim(basename)) // sdata(2:10) // &
                  "." // extension(i:i) // "s"
          else
             fname = basename(1:len_trim(basename)) // "-new" // sdata(2:10) // &
                  "." // extension(i:i)
          endif
          call print_message(fname(1:len_trim(fname)))
          
          if (w_spec .or. r_spec) then
             ! converting to/from spectral, keep same headers
             call singlefile_io3(time,Q,fname,work1,work2,0,io_pe,w_spec,header_user)
          else
             ! assume we are just converting to headerless brick-of-floats
             call singlefile_io3(time,Q,fname,work1,work2,0,io_pe,w_spec,2)
          endif
       enddo
    endif
    
    if (convert_opt==5) then  ! -cout norm
       ! 2048^3 needs 192GB * 1.66.  needs 256 cpus
       call input_uvw(time,Q,vor,work1,work2,header_user)
       do i=1,ndim
          call print_max(Q(1,1,1,i),i)
       enddo
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
      type_in=4
      write(message,'(f10.4)') 10000.0000 + time
      write(ext,'(f8.3)') 1000 + schmidt_in
      write(ext2,'(i3)') 100+type_in
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
           // message(2:10) // '.t' // ext2(2:3) // '.s' // ext(2:8)
      call print_message(rundir)
      call print_message(runname)
      call print_message(fname)
      ! time_tmp gets replaced by value -1 for headerless data and so we don't want to use time
      time_tmp=time	
      call singlefile_io3(time_tmp,Q,fname,work1,work2,1,io_pe,.false.,header_user)
      call global_min(Q,mn)
      call global_max(Q,mx)
      write(message,'(a,2f17.5)') 'passive scalar min/max: ',mn,mx
      call print_message(message)	
      
      write(message,'(f10.4)') 10000.0000 + time
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
           // 'raw_'// message(2:10) // '.t' // ext2(2:3) // '.s' // ext(2:8)
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
      
      !      do n=1,3
      !         call detrend_data(Q(1,1,1,n),0)
      !      enddo
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
      if (r_compressed) then
         !         skip this.  reader will print stats, and div(u)=0
         !         call print_stats(Q,vor,work1,work2)
      else
         call print_stats(Q,vor,work1,work2)
      endif
   endif
   
   
   if (convert_opt==12) then ! -cout uwbar  
      time_tmp=time
      call input_uvw(time_tmp,Q,vor,work1,work2,header_user)  
      call print_stats(Q,vor,work1,work2)
      ! compute the vorticity, store in vor(:,:,:,1)
      call vorticity2d(vor,Q,work1,work2)
      
      call print_message("computing uwbar...")
      ! overwrite Q with bar(u vor) - bar(u) bar(vor)
      call uw_filter(Q,vor(1,1,1,1),vor(1,1,1,2),work1,work2)
      
      call print_message("outputing uwbar...")
      write(sdata,'(f10.4)') 10000.0000 + time
      basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
           // sdata(2:10) // ".kct" 
      call singlefile_io3(time,Q(1,1,1,1),basename,work1,work2,0,io_pe,.false.,header_user)
      basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
           // sdata(2:10) // ".ct" 
      call singlefile_io3(time,Q(1,1,1,2),basename,work1,work2,0,io_pe,.false.,header_user)
   endif
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Low-pass filter or truncation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   if (convert_opt==13) then  ! -cout trunc
      if (w_spec) then
         ! -so (spectral output) uses spec_max to specify coefficients to write
         ! -cout trunc uses spec_max as a truncation parameter, so we cant
         ! use both of these together 
         call abortdns("-cout trunc option cant be used with spectral output")
      endif
      
      ! read data, header type =1, or specified in input file
      time_tmp=time
      call input_uvw(time_tmp,Q,vor,work1,work2,header_user)  
      if (.not. r_spec) then  ! r_spec reader will print stats, so we can skip this:
         call print_stats(Q,vor,work1,work2)
      endif
      
      ! compute the FFT
      do n=1,3
         write(message,'(a,i4)') 'w_spec fft3d: n=',n
         call print_message(message)
         call fft3d(Q(1,1,1,n),work1)
      enddo
      
      ! apply Filter
      do n=1,3
         ! call the truncation filter subroutine in fftops.F90
         call fft_filter_trunc(Q(1,1,1,n)) 
      enddo
      
      ! compute iFFT
      do n=1,3
         write(message,'(a,i4)') 'w_spec ifft3d: n=',n
         call print_message(message)
         call ifft3d(Q(1,1,1,n),work1)
      enddo
      
      ! give the output file a new name
      write(message, '(i5)') 10000 + spec_max
      basename=runname(1:len_trim(runname)) // "trunc"//message(2:5)//"_"
      
      call output_uvw(basename,time2,Q,vor,work1,work2,header_user)  
      ! output headerless data:
      ! call output_uvw(basename,time,Q,vor,work1,work2,2)
   endif
   
   
   if (convert_opt==14) then  ! -cout nlin
      time2=0
      w_spec = .true.
      basename=runname(1:len_trim(runname))
      
      ! compute two random fields
      call input1(Q(1,1,1,1),work2,work1,null,io_pe,.true.,-1)  
      call input1(Q(1,1,1,2),work2,work1,null,io_pe,.true.,-1)  
      call compute_nonlinear(Q,vor,work1,work2)
      do n=1,3
         call fft3d(Q(1,1,1,n),work1) ! take FFT so we can output spectral coefficients
      enddo
      ntot=0
      do k=nz1,nz2
         do j=ny1,ny2
            do i=nx1,nx2
               if (dealias_remove( abs(imcord(i)),abs(jmcord(j)), abs(kmcord(k))  )) then
               else
                  ntot=ntot+1
               endif
            enddo
         enddo
      enddo
      spec_max = g_nx  ! entire field, even if user specified spec_max on command line
      call output_uvw(basename,time2,Q,vor,work1,work2,header_user)  
      print *,'number of retained modes: ',ntot
   endif

   if (convert_opt==15) then  ! -cout nlout
      time2=0
      r_spec = .true.
      call input_uvw(time2,Q,vor,work1,work2,header_user)   ! read spec, but returns grid point values.
      ! compute nonlinear term 
      Q(:,:,:,3)=Q(:,:,:,1)*Q(:,:,:,2)
      call fft3d(Q(1,1,1,3),work1)
      vor(:,:,:,3)=Q(:,:,:,3)  ! save "exact" solution
      
      ! read in dealiased solution (again)
      call input_uvw(time2,Q,vor,work1,work2,header_user)   
      call fft3d(Q(1,1,1,3),work1)
      ! compare, mode by mode,  Q3 with saved Q3
      ntot=0
      nzero=0
      nerr=0
      mx=0
      xwerr=0
      do k=nz1,nz2
         do j=ny1,ny2
            do i=nx1,nx2
               ntot=ntot+1
               if (abs(Q(i,j,k,3))<1e-12) then
                  nzero=nzero+1
               else
                  a0 = abs(Q(i,j,k,3)-vor(i,j,k,3))
                  mx=max(mx,a0)
                  if (a0>1e-12) then
                     nerr=nerr+1
                     xw=sqrt(real(imcord(i)**2+jmcord(j)**2+kmcord(k)**2))
                     xwerr(int(xw))=max(xwerr(int(xw)),a0)
                     write(*,'(f6.2,3i4,3e16.7)') xw,&
                          imcord(i),jmcord(j),kmcord(k),&
                          Q(i,j,k,3),vor(i,j,k,3),a0
                  endif
               endif
            enddo
         enddo
      enddo
      do i=1,1000
         if (xwerr(i)/=0) then
            print *,'wave number=',i,'  max error=',xwerr(i)
         endif
      enddo
      print *,'number of modes:               ',ntot
      print *,'number of non-zero modes:      ',ntot-nzero
      print *,'max error over non-zero modes: ',mx
   endif
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  output delta-filtered u velocity component
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (convert_opt==16) then  ! -cout dfilter
      ! read data, header type =1, or specified in input file
      time_tmp=time
      call input_uvw(time_tmp,Q,vor,work1,work2,header_user)  
      if (.not. r_spec) then  ! r_spec reader will print stats, so we can skip this:
         call print_stats(Q,vor,work1,work2)
      endif
      
      do n=2,2  !  n=1,2,3 loops over (u,v,w) velocity components
         write(message,'(a,i4)') 'computing fft3d of input: n=',n
         call print_message(message)
         call fft3d(Q(1,1,1,n),work1)
         
         ! compute delta-filtered component, store in "vor()" array
         do k=1,g_nmin/3
            vor(:,:,:,n)=Q(:,:,:,n)
            call fft_filter_shell(vor(1,1,1,n),k) 
            call ifft3d(vor(1,1,1,n),work1)
            
            write(message,'(a,i4)') 'outputting delta filtered u for kshell=',k
            call print_message(message)
            
            write(sdata,'(f10.4)') 10000.0000 + time
            write(sdata2,'(i5)') 10000 + k
            basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
                 // sdata(2:10) // ".k" // sdata2(2:5)
            ! call singlefile_io3(time,vor(1,1,1,n),basename,work1,work2,0,io_pe,.false.,header_user)
            call singlefile_io3(time,vor(1,1,1,n),basename,work1,work2,0,io_pe,.false.,2)
            
         enddo
      enddo
   endif
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  output PDFs of delta-filtered u velocity component
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (convert_opt==17) then  ! -cout dpdf
      ! number of kshells:  
      kshell_max = nint(dealias_sphere_kmax)
      
      
      ! tell PDF module what we will be computing: 
      number_of_cpdf = 2*kshell_max  
      compute_uvw_pdfs = .true.    ! velocity increment PDFs
      compute_uvw_jpdfs = .false.    ! velocity increment joint PDFs
      compute_passive_pdfs = .false.  ! passive scalar PDFs
      
      ! read data, header type =1, or specified in input file
      time_tmp=time
      call input_uvw(time_tmp,Q,vor,work1,work2,header_user)  
      if (.not. r_spec) then  ! r_spec reader will print stats, so we can skip this:
         call print_stats(Q,vor,work1,work2)
      endif
      
      call global_max_abs(Q(1,1,1,1),mx)  ! max of abs(U) component
      uscale = mx/100     ! binsize for velocity increment PDFs
      epsscale = uscale   ! binsize for epsilon
      call init_pdf_module()
      
      
      call print_message("computing velocity increment PDFs")
      call compute_all_pdfs(Q,vor,work1)
      
      do n=2,2  !  n=1,2,3 loops over (u,v,w) velocity components
         write(message,'(a,i4)') 'computing fft3d of input: n=',n
         call print_message(message)
         call fft3d(Q(1,1,1,n),work1)
         
         ! compute delta-filtered component, store in "vor()" array
         do k=1,kshell_max
            k2=kshell_max+k
            
            work2 = Q(:,:,:,n)
            call fft_filter_shell(work2,k)
            call ifft3d(work2,work1)
            
            ! compute PDFs
            call global_max_abs(work2,mx)
            binsize = mx/100   ! should produce about 200 bins
            
            write(message,'(a,i4,a,e10.3,a,e10.3)') 'PDF delta filtered k=',k,' max|u|=',mx,' binsize=',binsize
            call print_message(message)
            call compute_pdf_scalar(work2,cpdf(k),binsize)
            
            
            work2 = Q(:,:,:,n)
            call fft_filter_shell1(work2,k)
            call ifft3d(work2,work1)
            
            ! compute PDFs
            call global_max_abs(work2,mx)
            binsize = mx/100   ! should produce about 200 bins
            
            write(message,'(a,i4,a,e10.3,a,e10.3)') 'PDF simplified delta filtered k=',k,' max|u|=',mx,' binsize=',binsize
            call print_message(message)
            call compute_pdf_scalar(work2,cpdf(k2),binsize)

            
         enddo
      enddo
      !
      ! now output the PDF data
      !
      write(sdata,'(f10.4)') 10000.0000 + time
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) // ".cpdf"
      call print_message(fname)
      if (my_pe==io_pe) then
         call copen(fname,"w",fid,ierr)
         if (ierr/=0) then
            write(message,'(a,i5)') "output_model(): Error opening .cpdf file errno=",ierr
            call abortdns(message)
         endif
      endif
      
      write(sdata,'(f10.4)') 10000.0000 + time
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) // ".sf"
      call print_message(fname)
      if (my_pe==io_pe) then
         call copen(fname,"w",fidu,ierr)
         if (ierr/=0) then
            write(message,'(a,i5)') "output_model(): Error opening .sf file errno=",ierr
            call abortdns(message)
         endif
      endif
      call output_pdf(time,fidu,NULL,NULL,fid,NULL)
      
      if (my_pe==io_pe) call cclose(fid,ierr)
      if (my_pe==io_pe) call cclose(fidu,ierr)
      
   endif
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! high-pass filter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (convert_opt==18) then  ! -cout hpass
      if (w_spec) then
         ! -so (spectral output) uses spec_max to specify coefficients to write
         ! -cout trunc uses spec_max as a truncation parameter, so we cant
         ! use both of these together 
         call abortdns("-cout hpass option cant be used with spectral output")
      endif
      
      ! read data, header type =1, or specified in input file
      time_tmp=time
      !     DNS-headered data input
      !     call input_uvw(time2,Q,vor,work1,work2,header_user)  
      !     headerless data input
      call input_uvw(time_tmp,Q,vor,work1,work2,2)  
      
      if (.not. r_spec) then  ! r_spec reader will print stats, so we can skip this:
         call print_stats(Q,vor,work1,work2)
      endif
      
      ! compute the FFT
      do n=1,3
         write(message,'(a,i4)') 'w_spec fft3d: n=',n
         call print_message(message)
         call fft3d(Q(1,1,1,n),work1)
      enddo
      
      ! apply Filter
      do n=1,3
         ! call the truncation filter subroutine in fftops.F90
         call fft_filter_hpass(Q(1,1,1,n)) 
      enddo
      
      ! compute iFFT
      do n=1,3
         write(message,'(a,i4)') 'w_spec ifft3d: n=',n
         call print_message(message)
         call ifft3d(Q(1,1,1,n),work1)
      enddo
      
      ! give the output file a new name
      write(message, '(i5)') 10000 + spec_max
      basename=runname(1:len_trim(runname)) // "hpass"//message(2:5)//"."
      
      !output headered data:
      !call output_uvw(basename,time2,Q,vor,work1,work2,header_user)  
      !output headerless data:
      call output_uvw(basename,time,Q,vor,work1,work2,2)
   endif
   
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! coarse grain  coarsegrain
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (convert_opt==19) then  ! -cout coarsegrain
      ! read data, header type =1, or specified in input file
      if (.not.allocated(Q2)) allocate(Q2(nx,ny,nz,n_var))
      time_tmp=time
      call input_uvw(time_tmp,Q,vor,work1,work2,header_user)  
      if (.not. r_spec) then  ! r_spec reader will print stats, so we can skip this:
         call print_stats(Q,vor,work1,work2)
      endif
      
      do j=1,size(Mval)
         do n=1,3
            ! call the truncation filter subroutine in fftops.F90
            call coarse_grain(Q(1,1,1,n),Q2(1,1,1,n),work1,Mval(j) ) 
         enddo
         
         ! give the output file a new name
         write(message, '(i5)') 10000 + Mval(j)
         basename=runname(1:len_trim(runname)) // "-cg"//message(2:5)//"."
         
         call output_uvw(basename,time2,Q2,vor,work1,work2,header_user)  
         ! output headerless data:
         ! call output_uvw(basename,time,Q,vor,work1,work2,2)
      enddo
   endif
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! output du/dx and PDF of du/dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (convert_opt==20) then  ! -cout dudx
      call input_uvw(time,Q,vor,work1,work2,header_user)
      call print_message("computing du/dx...")
      
      
      ! output vorticity magnitude
      write(sdata,'(f10.4)') 10000.0000 + time
      basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))
      
      call der(Q(1,1,1,3),vor,work1,work2,DX_ONLY,3)
      fname = basename(1:len_trim(basename)) // sdata(2:10) // ".dwdz"
      call singlefile_io3(time,vor,fname,work1,work2,0,io_pe,.false.,2)
      
      
      call der(Q(1,1,1,2),vor,work1,work2,DX_ONLY,2)
      fname = basename(1:len_trim(basename)) // sdata(2:10) // ".dvdy"
      call singlefile_io3(time,vor,fname,work1,work2,0,io_pe,.false.,2)
      
      call der(Q(1,1,1,1),vor,work1,work2,DX_ONLY,1)
      fname = basename(1:len_trim(basename)) // sdata(2:10) // ".dudx"
      call singlefile_io3(time,vor,fname,work1,work2,0,io_pe,.false.,2)
      
      
      number_of_cpdf = 2
      compute_uvw_pdfs = .false.    ! velocity increment PDFs
      compute_uvw_jpdfs = .false.    ! velocity increment joint PDFs
      compute_passive_pdfs = .false.  ! passive scalar PDFs
      call init_pdf_module()
      
      
      ! compute PDF of du/dx
      call global_max_abs(work1,mx)  ! compute max, used to determine binsize
      binsize = mx/100   ! should produce about 200 bins
      call compute_pdf_scalar(work1,cpdf(1),binsize)
      
      
      ! compute epsilon, store in work1.  
      !      call hyperder(Q,work1,work2)
      call coshder(Q(1,1,1,1),vor,work2)	
      work1=vor(:,:,:,1)*Q(:,:,:,1)
      call coshder(Q(1,1,1,2),vor,work2)	
      work1=work1+vor(:,:,:,1)*Q(:,:,:,2)
      call coshder(Q(1,1,1,3),vor,work2)	
      vor(:,:,:,1)=work1+vor(:,:,:,1)*Q(:,:,:,3)
      call global_max_abs(vor,mx)
      binsize = mx/100   ! should produce about 200 bins
      call compute_pdf_scalar(vor,cpdf(2),binsize)
      
      ! output coshder field
      fname = basename(1:len_trim(basename)) // sdata(2:10) // ".coshder"
      call singlefile_io3(time,vor,fname,work1,work2,0,io_pe,.false.,2)
      
      
      write(sdata,'(f10.4)') 10000.0000 + time
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) // ".dudxpdf"
      call print_message(fname)
      if (my_pe==io_pe) then
         call copen(fname,"w",fid,ierr)
         if (ierr/=0) then
            write(message,'(a,i5)') "output_model(): Error opening .cpdf file errno=",ierr
            call abortdns(message)
         endif
      endif
      
      call output_pdf(time,NULL,NULL,NULL,fid,NULL)
   endif
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write out headerless horizontal velocity field in physical space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   if (convert_opt==21) then  ! -cout ehor
      ! 2048^3 needs 192GB * 1.66.  needs 256 cpus
      time_tmp=time
      call input_uvw(time_tmp,Q,vor,work1,work2,header_user)
      do i=1,ndim
         call print_max(Q(1,1,1,i),i)
      enddo
      call print_message("computing horizontal energy...")
      do k=nz1,nz2
         do j=ny1,ny2
            do i=nx1,nx2
               work1(i,j,k)=0.5*(Q(i,j,k,1)**2+Q(i,j,k,2)**2)
            enddo
         enddo
      enddo
      write(sdata,'(f10.4)') 10000.0000 + time
      basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))
      fname = basename(1:len_trim(basename)) // sdata(2:10) // ".ehor"
      call singlefile_io3(time,work1,fname,Q,work2,0,io_pe,.false.,2)
   endif
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Write out headerless potential vorticity and potential enstrophy field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (convert_opt==22) then  ! -cout potens
      write(message,'(f10.4)') 10000.0000 + time
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
           // message(2:10) // '.pv'
      call singlefile_io3(time,Q,fname,work1,work2,1,io_pe,.false.,1)
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
           // message(2:10) // '.potvor'
      call print_message("writing headerless potential vorticity...")
      call singlefile_io3(time,Q,fname,work1,work2,0,io_pe,.false.,2)
      call print_message("computing potential enstrophy...")	
      do k=nz1,nz2
         do j=ny1,ny2
            do i=nx1,nx2
               work1(i,j,k)=0.5*(Q(i,j,k,1)**2)
            enddo
         enddo
      enddo
      
      write(message,'(f10.4)') 10000.0000 + time
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
           // message(2:10) // '.pv2'
      call print_message(fname)
      call singlefile_io3(time,work1,fname,Q,work2,0,io_pe,.false.,2)
   endif
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Write out headerless potential energy field in physical space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   if (convert_opt==23) then  ! -cout pe
      ! 2048^3 needs 192GB * 1.66.  needs 256 cpus
      time_tmp=time
      call input_uvw(time_tmp,Q,vor,work1,work2,header_user)
      time_tmp=time
      call input_passive(runname, time_tmp, Q, work1, work2,header_user)
      !      do i=1,ndim
      !         call print_max(Q(1,1,1,i),i)
      !      enddo
      call print_message("computing potential energy...")
      do k=nz1,nz2
         do j=ny1,ny2
            do i=nx1,nx2
               work1(i,j,k)=0.5*(Q(i,j,k,np1)**2)
            enddo
         enddo
      enddo
      write(sdata,'(f10.4)') 10000.0000 + time
      basename=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))
      fname = basename(1:len_trim(basename)) // sdata(2:10) // ".pe"
      call singlefile_io3(time,work1,fname,Q(1,1,1,np1),work2,0,io_pe,.false.,2)
   endif
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Write out wave modes field and vortical modes field in physical space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   if (convert_opt==24) then  ! -cout CH_decomp
      schmidt_in=1.0
      type_in=4
      time_tmp=time
      call input_uvw(time_tmp,Q,vor,work1,work2,header_user)

      time_tmp=time 
      call input_passive(runname,time_tmp,Q,work1,work2,header_user)

      if (.not.allocated(Q2)) allocate(Q2(nx,ny,nz,n_var))


      !      do i=1,ndim
      !         call print_max(Q(1,1,1,i),i)
      !      enddo
      call print_message("computing wave and vortical projected fields ...")
      call compute_project_CHfield(Q,Q2,work1,work2)	
      
      ! check energy of wave field
      dummy=0
      do n = 1,n_var
         do k=nz1,nz2
            do j=ny1,ny2
               do i=nx1,nx2
                  dummy= dummy+0.5*(Q(i,j,k,n)**2)
               enddo
            enddo
         enddo
      enddo
#ifdef USE_MPI
      xtmp = dummy
      call mpi_allreduce(xtmp,dummy,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
      dummy=dummy/g_nx/g_ny/g_nz
      write(message,'(a,e15.8)') "grid space energy in wave field = ",dummy
      call print_message(message)
      
      
      ! check energy of vortical field
      dummy=0
      do n = 1,n_var
         do k=nz1,nz2
            do j=ny1,ny2
               do i=nx1,nx2
                  dummy= dummy+0.5*(Q2(i,j,k,n)**2)
               enddo
            enddo
         enddo
      enddo
#ifdef USE_MPI
      xtmp = dummy
      call mpi_allreduce(xtmp,dummy,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
      dummy=dummy/g_nx/g_ny/g_nz
      write(message,'(a,e15.8)') "grid space energy in vortical field = ",dummy
      call print_message(message)
      
      write(sdata,'(f10.4)') 10000.0000 + time
      write(ext,'(f8.3)') 1000 + schmidt_in
      write(ext2,'(i3)') 100+type_in
      
      !     write out headerless wave component u,v,w,t fields
      basename=runname(1:len_trim(runname)) // 'wave_'
      call output_uvw(basename,time,Q,vor,work1,work2,2)
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
           //'wave_'// sdata(2:10) // '.t' // ext2(2:3) // '.s' // ext(2:8)
      call singlefile_io3(time,Q(1,1,1,np1),fname,work1,work2,0,io_pe,.false.,2)
      write(sdata,'(f10.4)') 10000.0000 + time
      
      !     write out headerless vortical component u,v,w,t fields
      basename=runname(1:len_trim(runname)) // 'vort_'
      call output_uvw(basename,time,Q2,vor,work1,work2,2) 
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
           // 'vort_' // sdata(2:10) // '.t' // ext2(2:3) // '.s' // ext(2:8) 
      call singlefile_io3(time,Q2(1,1,1,np1),fname,work1,work2,0,io_pe,.false.,2)
      
   endif
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Low-pass filter or truncation for scalar field 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   if (convert_opt==25) then  ! -cout ttrunc
      if (w_spec) then
         ! -so (spectral output) uses spec_max to specify coefficients to write
         ! -cout trunc uses spec_max as a truncation parameter, so we cant
         ! use both of these together 
         call abortdns("-cout trunct option cant be used with spectral output")
      endif
      
      ! read data, header type =1, or specified in input file
      time2=time
      call input_passive(runname, time2, Q, work1, work2,header_user)
      if (.not. r_spec) then  ! r_spec reader will print stats, so we can skip this:
         call print_stats(Q,vor,work1,work2)
      endif
      
      ! compute the FFT
      write(message,'(a,i4)') 'w_spec fft3d: n=',np1
      call print_message(message)
      call fft3d(Q(1,1,1,np1),work1)
      
      ! apply Filter
      ! call the truncation filter subroutine in fftops.F90
      call fft_filter_trunc(Q(1,1,1,np1)) 
      
      ! compute iFFT
      write(message,'(a,i4)') 'w_spec ifft3d: n=',np1
      call print_message(message)
      call ifft3d(Q(1,1,1,np1),work1)
      
      ! give the output file a new name
      write(message, '(i5)') 10000 + spec_max
      write(sdata,'(f10.4)') 10000.0000 + time
      write(ext,'(f8.3)') 1000 + schmidt_in
      write(ext2,'(i3)') 100+type_in
      
      
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
           //'trunc'//message(2:5)//"_"//sdata(2:10) // '.t' // ext2(2:3) // '.s' // ext(2:8)
      call singlefile_io3(time,Q(1,1,1,np1),fname,work1,work2,0,io_pe,.false.,1)
   endif
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  output PDFs of zero-crossings of velocity components in all cartesian
!! directions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (convert_opt==26) then  ! -cout zerocr
      
      ! tell PDF module what we will be computing: 
      number_of_cpdf = 0           
      compute_uvw_pdfs = .false.    ! velocity increment PDFs
      compute_uvw_jpdfs = .false.    ! velocity increment joint PDFs
      compute_passive_pdfs = .false.  ! passive scalar PDFs
      compute_zerocrossing_pdfs = .true. ! zero crossing PDFs      
      
      
      ! read data, header type =2 for headerless data or specified in input file
      time_tmp=time
      call input_uvw(time_tmp,Q,vor,work1,work2,2)  
      if (.not. r_spec) then  ! r_spec reader will print stats, so we can skip this:
         call print_stats(Q,vor,work1,work2)
      endif
      
      mx = max(1.0d0*g_nx,1.0d0*g_ny)  ! max grid points is the maximum number of crossings
      mx = max(mx,g_nz*1.0d0)
      if (my_pe==io_pe) print *, 'max num zero-crosssing',mx
      zcscale = 2     ! binsize for number of zerocrossings
      call init_pdf_module()
      
      
      call print_message("computing zero-crossing PDFs")
      call compute_all_pdfs(Q,vor,work1)
      

      !
      ! now output the PDF data
      !
      write(sdata,'(f10.4)') 10000.0000 + time
      fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) // ".zcl_pdf"
      call print_message(fname)
      if (my_pe==io_pe) then
         call copen(fname,"w",fid,ierr)
         if (ierr/=0) then
            write(message,'(a,i5)') "output_model(): Error opening .zcpdf file errno=",ierr
            call abortdns(message)
         endif
      endif
      
      call output_pdf(time,fid,NULL,NULL,NULL,NULL)
      
      if (my_pe==io_pe) call cclose(fid,ierr)
      
   endif
   
   
   
   if (tstart>=0) then
      time=time+tinc
      if (time > max(tstop,tstart)) exit
      if (time < min(tstop,tstart)) exit
   endif
enddo
if (io_pe==my_pe) then
   print *,'convert() finished'
   print *,'tstart = ',tstart
   print *,'tstop  = ',tstop
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
     call abortdns('Error: insufficient ghost cells in x direction')
if (ny1<3 .or. ny2+2>ny)&
     call abortdns('Error: insufficient ghost cells in y direction')
if (nz1<3 .or. nz2+2>nz) then
     call abortdns('Error: insufficient ghost cells in z direction')
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
     call abortdns('Error: insufficient ghost cells in x direction')
if (ny1<3 .or. ny2+2>ny)&
     call abortdns('Error: insufficient ghost cells in y direction')
if (nz1<3 .or. nz2+2>nz)&
     call abortdns('Error: insufficient ghost cells in z direction')

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
      call abortdns("error: largest <gradu> element is on diagional")
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



subroutine print_max(p,n)
use params
use mpi
implicit none
integer n,ierr
real*8 :: p(nx,ny,nz)
real*8 :: mx

call global_max_abs(p,mx)
if (my_pe==io_pe) then
   write(*,'(a,i2,a,e20.15)') 'n=',n,' max abs val =',mx
endif
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


!call vorticity(div,Q,work1,work2)
ens=0
ke=0
do n=1,3
do k=nz1,nz2
   do j=ny1,ny2
      do i=nx1,nx2
!         ens=ens + div(i,j,k,n)**2
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
   call global_max_abs(Q(1,1,1,n),mx(n))
enddo 

write(message,'(a,3f18.14)') 'STATS:  maxU = ',mx
call print_message(message)

call compute_div(Q,div,work1,work2,divx,divi)
write(message,'(3(a,e12.5))') 'STATS:  max(div)=',divx
call print_message(message)	


!write(message,'(3(a,f18.12))') 'STATS:  enstrophy=',ens
!call print_message(message)	




end subroutine









subroutine compute_nonlinear(Q,Qhat,work,work_hat)
use params
implicit none
real*8 Q(nx,ny,nz,n_var)
real*8 Qhat(g_nz2,nx_2dz,ny_2dz,n_var)           ! Fourier data at time t
real*8 work_hat(g_nz2,nx_2dz,ny_2dz)
real*8 work(nx,ny,nz)
! local
integer :: n,i,j,k,im,jm,km



! compute fourier coefficieints
do n=1,2
   call z_fft3d_trashinput(Q(1,1,1,n),Qhat(1,1,1,n),work)
enddo

! apply dealias filter, and additional filter if user specifed -smax X
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nx_2dz
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

         ! dealias           
         if ( dealias_remove(abs(im),abs(jm),abs(km))) then
            Qhat(k,i,j,1)=0
            Qhat(k,i,j,2)=0
         endif
         if ( spec_max>0 .and. im**2+jm**2+km**2 > spec_max**2 ) then
            Qhat(k,i,j,1)=0
            Qhat(k,i,j,2)=0
         endif


      enddo
   enddo
enddo
! compute nonlinear term in grid space
do n=1,2
   call z_ifft3d(Qhat(1,1,1,n),Q(1,1,1,n),work)
enddo


Q(:,:,:,3)=Q(:,:,:,1)*Q(:,:,:,2)
! back to spectral space
call z_fft3d_trashinput(Q(1,1,1,3),Qhat(1,1,1,3),work)


if (use_phaseshift) then
   ! compute nonlinear product with phase shifted values, then
   ! take the average of this calculation and result above

   ! phase shift and FFT Q2 into Q3
   call z_phaseshift(Qhat(1,1,1,2),1,work)    ! phaseshift Qhat
   call z_ifft3d(Qhat(1,1,1,2),Q(1,1,1,3),work)

   ! phase shift and FFT Q1 into work
   call z_phaseshift(Qhat(1,1,1,1),1,work)    ! phaseshift Qhat
   call z_ifft3d(Qhat(1,1,1,1),work,work_hat)

   Q(:,:,:,3)=Q(:,:,:,3)*work(:,:,:)
   call z_fft3d_trashinput(Q(1,1,1,3),work_hat,work)
   call z_phaseshift(work_hat,-1,work)  ! un-phaseshift result
   Qhat(:,:,:,3) = .5*Qhat(:,:,:,3) + .5*work_hat(:,:,:)
endif



! apply filter to nonlinear product
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nx_2dz
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)
         ! dealias           
!         if (  abs(Qhat(k,i,j,3))>1e-10 .or. abs(work_hat(k,i,j))>1e-10 ) then
!            print *,'uv: ',im,jm,km,Qhat(k,i,j,3),work_hat(k,i,j)
!         endif


         if ( dealias_remove(abs(im),abs(jm),abs(km))) then
            Qhat(k,i,j,3)=0
         endif
      enddo
   enddo
enddo

! return filtered initial condition and final answer in grid space
call z_ifft3d(Qhat(1,1,1,3),Q(1,1,1,3),work)


end subroutine




