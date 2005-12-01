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
program anal_spec2d
use params
use mpi
use fft_interface
implicit none
real*8 :: spec_2d(nx+1,nz/2+1)
character(len=80) message,sdata
character(len=280) basename,fname,fnamewrite
integer ierr,i,j,k,n,km,im,jm,icount,ictof,ivorsave,ntrunc
integer icentroid,iscalars,i3dspec,i2dspec,ireaduvh
real*8 :: tstart,tstop,tinc,time,time2
real*8 :: u,v,w,x,y
real*8 :: kr,ke,ck,xfac
CPOINTER :: fspec2d
real*8, allocatable :: spec2d(:,:)
integer nkh, nkz

allocate(spec2d(130,130))

! input file
basename="r16"
if (my_pe==io_pe) print *,basename
tstart=10.5
tstop=10.5
tinc=.1
ireaduvh = 0 ! =0 Do not read in the u,v,h fields
icentroid=0 ! Compute the KE centroid? 0 = no, 1 = yes
icount=0  ! A counter for averaging to compute time averages
ictof = 0 ! If 1 then flip endian sign
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
!call set_byteswap_input(1);

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!Read in and plot the 2d spectrum
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(i3dspec == 1) then
   if (my_pe==io_pe)  then
      call copen("/scratch/wingate/r16/r160003.5000.spec2d","r",fspec2d,ierr)
      print *,'opened the file ierr=',ierr
      do 
         call cread8e(fspec2d,time,1,ierr)
         print *,'read first number  time=',time,' ierr=',ierr
         if (ierr /= 1) exit
         call cread8e(fspec2d,x,1,ierr)
         nkh = x
         call cread8e(fspec2d,x,1,ierr)
         nkz = x
         write(6,*) "time = ",time, "nkh = ",nkh, "nkz = ",nkz
         call cread8e(fspec2d,spec2d,nkh*nkz,ierr)
      enddo

   endif
endif ! if(i3dspec == 1) then



call close_mpi
end program anal_spec2d


