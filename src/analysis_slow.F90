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
program anal
use params
use mpi
use fft_interface
implicit none
real*8,save  :: Q(nx,ny,nz,n_var)
real*8,save  :: Qslow(nx,ny,nz,n_var)
real*8,save  :: vor2d(nx,ny)
real*8,save  :: dx(nx,ny,nz),dy(nx,ny,nz)
real*8,save  :: div(nx,ny,nz)
real*8,save  :: potvor(nx,ny,nz)
real*8,save  :: div2d(nx,ny)
real*8,save  :: work1(nx,ny,nz)
real*8,save  :: work2(nx,ny,nz)
real*8,save  :: dxx(1) ! dummy varray 
!real*8,save  :: vorfilt(nx,ny,nz),pvfilt(nx,ny,nz)
character(len=80) message,sdata
character(len=280) basename,fname,fnamewrite
integer ierr,i,j,k,n,km,im,jm,icount,ictof,ivorsave,ntrunc
integer icentroid,iscalars,i3dspec,i2dspec,ireaduvh
real*8 :: tstart,tstop,tinc,time,time2
real*8 :: u,v,w,x,y, n_r, spec_r
real*8 :: kr,ke,ck,xfac,vor_slow,pe_slow,kev_slow,keh_slow,kew_slow,potens_slow
integer ni,ns,nb


! input file
basename="n100f1"
if (my_pe==io_pe) print *,basename
tstart=0.0
tstop=0.0
tinc=.1
ireaduvh = 0 ! =0 Do not read in the u,v,h fields
icentroid=0 ! Compute the KE centroid? 0 = no, 1 = yes
icount=0  ! A counter for averaging to compute time averages
ictof = 1 ! If 1 then flip endian sign
ivorsave=0 ! If 1 then compute an averaged vorticity
ntrunc = 0 ! If not zero find Q with wave numbers ntrunc and above
iscalars = 0 ! If not zero read in scalars (for energy as a function&
             ! of time, etc)
i3dspec = 1 ! = 1 Read in and write out a plot for the 3d spectrum
            ! = 0 Do nothing
            ! = 2 Compute the 2d spectrum from the 3d spectrum
i2dspec = 0 ! = 1 Read in the 2d spectrum and plot it out
            ! = 0 Do nothing
call init_mpi       
call init_mpi_comm3d()
call init_model

! byte swap the input data:
call set_byteswap_input(1);


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Further Initialization such as arrays and opening files
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Q=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! File for the slow conserved quantities 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     open(55,file='keh_slow.txt')
     open(56,file='kev_slow.txt')
     open(57,file='vor_slow.txt')
     open(58,file='pe_slow.txt')
     open(59,file='potens_slow.txt')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Read in time dependent files
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   time=tstart
   do
      icount=icount+1 
      if(my_pe ==io_pe) print *,'icount = ',icount

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         
!         read in the Q array
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(sdata,'(f10.4)') 10000.0000 + time
      if(ireaduvh /= 0) then
         fname = basename(1:len_trim(basename)) // sdata(2:10) // ".u"
         if (my_pe==io_pe) then
            print *,'filename: ',fname(1:len_trim(fname))
         endif
         call singlefile_io(time2,Q(1,1,1,1),fname,work1,work2,1,io_pe)
         if (ictof == 1) then
            fnamewrite = basename(1:len_trim(basename)) // sdata(2:10) // ".Fu"
            open(42,file=fnamewrite,status='unknown',action='write',form='unformatted')
            write(42) Q(nx1:nx2,ny1:ny2,nz1:nz2,1)
         endif
         
         fname = basename(1:len_trim(basename)) // sdata(2:10) // ".v"
         call singlefile_io(time2,Q(1,1,1,2),fname,work1,work2,1,io_pe)
         if (ictof == 1) then
            fnamewrite = basename(1:len_trim(basename)) // sdata(2:10) // ".Fv"
            open(42,file=fnamewrite,status='unknown',action='write',form='unformatted')
            write(42) Q(nx1:nx2,ny1:ny2,nz1:nz2,1)
         endif
         
         fname = basename(1:len_trim(basename)) // sdata(2:10) // ".w"
         call singlefile_io(time2,Q(1,1,1,3),fname,work1,work2,1,io_pe)
         if (ictof == 1) then
            fnamewrite = basename(1:len_trim(basename)) // sdata(2:10) // ".Fw"
            open(42,file=fnamewrite,status='unknown',action='write',form='unformatted')
            write(42) Q(nx1:nx2,ny1:ny2,nz1:nz2,1)
         endif

         fname = basename(1:len_trim(basename)) // sdata(2:10) // ".t04.s001.000"
         call singlefile_io(time2,Q(1,1,1,4),fname,work1,work2,1,io_pe)
         if (ictof == 1) then
            fnamewrite = basename(1:len_trim(basename)) // sdata(2:10) // ".Ft04.s001.000"
            open(42,file=fnamewrite,status='unknown',action='write',form='unformatted')
            write(42) Q(nx1:nx2,ny1:ny2,nz1:nz2,1)
         endif
      else
         if(my_pe == io_pe) write(6,*) "You need to read in u v w theta."
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         Compute the slow variables for Ro -> 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        First compute the horizontally divergence free part to
!        subtract from the velocity
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call der(Q(1,1,1,1),dx,dxx,work1,DX_ONLY,1)
      call der(Q(1,1,1,2),dy,dxx,work1,DX_ONLY,2)
      div = dx+dy
      call fft3d(div,work1)
      !
      !        Compute the inverse horizontal laplacian in 2 dimension
      !
      do k=nz1,nz2
         do i=nx1,nx2
            do j=ny1,ny2
               jm=(jmcord(j))
               im=(imcord(i))
               div(i,j,k) = -div(i,j,k)/((im**2+jm**2)*pi2_squared)
            enddo
         enddo
      enddo
      call ifft3d(div,work1)
      !
      !        Now take derivatives of this to subtract off of each component
      !
      call der(div,dx,dxx,work1,DX_ONLY,1)
      call der(div,dy,dxx,work1,DX_ONLY,2)
      !
      Qslow = 0.
      do k=nz1,nz2
         do j=ny1,ny2
            do i=nx1,nx2
               Qslow(i,j,k,1) = Q(i,j,k,1) - dx(i,j,k)
               Qslow(i,j,k,2) = Q(i,j,k,2) - dy(i,j,k)
               Qslow(i,j,k,3) = Q(i,j,k,3)
               Qslow(i,j,k,4) = Q(i,j,k,4)
            enddo
         enddo
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Compute the slow variables by integrating u, v, and w in the vertical
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call zaverage(Qslow(1,1,1,1),dx)
      Qslow(:,:,:,1)=dx
      call zaverage(Qslow(1,1,1,2),dx)
      Qslow(:,:,:,2)=dx
      call zaverage(Qslow(1,1,1,3),dx)
      Qslow(:,:,:,3)=dx


      ! compute potential vorticity, (pvtype=1) 
      call potential_vorticity(potvor,div,Qslow,work1,work2,1)


      ! Compute the horizontal divergence of the slow variables
      call der(Qslow(1,1,1,1),dx,dxx,work1,DX_ONLY,1)
      call der(Qslow(1,1,1,2),dy,dxx,work1,DX_ONLY,2)
      do i=nx1,nx2
         do j=ny1,ny2
            div2d(i,j) = dx(i,j,nz1) + dy(i,j,nz1)
         enddo
      enddo

!     compute horizontal vorticity of slow variablles
      call der(Qslow(1,1,1,1),dy,dxx,work1,DX_ONLY,2)
      call der(Qslow(1,1,1,2),dx,dxx,work1,DX_ONLY,1)
      do i=nx1,nx2
         do j=ny1,ny2
            vor2d(i,j) = dx(i,j,nz1)-dy(i,j,nz1)
         enddo
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        Write out files for the vertical vorticity and the
!        horizontal divergence
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         if (ictof == 1) then
!            fnamewrite = basename(1:len_trim(basename)) // sdata(2:10) // ".Fvor"
!            open(42,file=fnamewrite,status='unknown',action='write',form='unformatted')
!            write(42) vor(nx1:nx2,ny1:ny2,nz1:nz2,3)
!            fnamewrite = basename(1:len_trim(basename)) // sdata(2:10) // ".Fpvor"
!            open(42,file=fnamewrite,status='unknown',action='write',form='unformatted')
!            write(42) vor2d(nx1:nx2,ny1:ny2,nz1:nz2)
!         endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         
!         compute the conserved quantities
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         First the 2D quantities
!             
      keh_slow = 0.
      kev_slow = 0.
      vor_slow = 0.
      do j=ny1,ny2
         do i=nx1,nx2
            vor_slow = vor_slow + vor2d(i,j)**2
            keh_slow = keh_slow + (Qslow(i,j,nz1,1)**2+Qslow(i,j,nz1,2)**2)
            kev_slow = kev_slow + Qslow(i,j,nz1,3)**2
         enddo
      enddo
      keh_slow=.5*keh_slow/g_nx/g_ny
      kew_slow=.5*kev_slow/g_nx/g_ny
      vor_slow = .5*vor_slow/g_nx/g_ny
      !
      !        Then the 3d quantities
      !
      potens_slow = 0.
      pe_slow=0.
      do k=nz1,nz2
         do j=ny1,ny2
            do i=nx1,nx2
               pe_slow = pe_slow + (1. + bous*zcord(k)*Qslow(i,j,k,4))**2
               potens_slow = potens_slow + .5*potvor(i,j,k)**2
            enddo
         enddo
      enddo

      pe_slow = .5*pe_slow/g_nz/g_nx/g_ny
      potens_slow = .5*potens_slow/g_nz/g_nx/g_ny


#ifdef USE_MPI
      xfac=keh_slow
      call mpi_allreduce(xfac,keh_slow,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
      xfac=kew_slow
      call mpi_allreduce(xfac,kew_slow,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
      xfac=vor_slow
      call mpi_allreduce(xfac,vor_slow,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
      xfac=pe_slow
      call mpi_allreduce(xfac,pe_slow,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
      xfac=potens_slow
      call mpi_allreduce(xfac,potens_slow,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif


      
#if 0
      fnamewrite = basename(1:len_trim(basename)) // sdata(2:10) // ".Tvor"
      call singlefile_io(time2,vorfilt,fnamewrite,work1,work2,0,io_pe)
      fnamewrite = basename(1:len_trim(basename)) // sdata(2:10) // ".FTvor"
      open(42,file=fnamewrite,status='unknown',action='write',form='unformatted')
      write(42) vorfilt(nx1:nx2,ny1:ny2,nz1:nz2)
      
      fnamewrite = basename(1:len_trim(basename)) // sdata(2:10) // ".Tpvor"
      call singlefile_io(time2,pvfilt,fnamewrite,work1,work2,0,io_pe)
      fnamewrite = basename(1:len_trim(basename)) // sdata(2:10) // ".FTpvor"
      open(42,file=fnamewrite,status='unknown',action='write',form='unformatted')
      write(42) pvfilt(nx1:nx2,ny1:ny2,nz1:nz2)
#endif
         
      if (my_pe==io_pe) then
         print *,'ke = ',ke
         write(55,'(e14.6,2x,e14.6)') real(time),keh_slow
         write(56,'(e14.6,2x,e14.6)') real(time),kev_slow
         write(57,'(e14.6,2x,e14.6)') real(time),vor_slow
         write(58,'(e14.6,2x,e14.6)') real(time),pe_slow
         write(59,'(e14.6,2x,e14.6)') real(time),potens_slow
      endif
      
      
      time=time+tinc
      if (time>tstop) exit
   enddo ! DO
   
   
   close(55)
   if (ivorsave /= 0) then
   endif
   
   call close_mpi
end program
   

