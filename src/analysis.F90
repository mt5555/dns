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
use structf
use fft_interface
implicit none
real*8,save  :: Q(nx,ny,nz,n_var)
real*8,save  :: Qhat(nx,ny,nz,n_var)
real*8,save  :: vor(nx,ny,nz,n_var)
real*8,save  :: vor2d(nx,ny,nz)
real*8,save  :: dx(nx,ny,nz)
real*8,save  :: dxx(nx,ny,nz)
real*8,save  :: div(nx,ny,nz)
real*8,save  :: work1(nx,ny,nz)
real*8,save  :: work2(nx,ny,nz)
character(len=80) message,sdata
character(len=280) basename,fname
integer ierr,i,j,k,n,km,im,jm
real*8 :: tstart,tstop,tinc,time,time2
real*8 :: u,v,w

! input file
basename="kh/khP"
tstart=0
tstop=.75
tinc=.25

call init_mpi       
call init_mpi_comm3d()
call init_grid      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if needed, initialize some constants.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Q=0
if (init_cond==3) call init_data_sht(Q,vor,work1,work2)

time=tstart
do
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

   ! compute vorticity 
   call vorticity(vor,Q,work1,work2)
   ! for 2D problems, 1st and 2nd components of vor() will be zero
   vor2d(:,:,:)=vor(:,:,:,3)


   ! compute divergence
   !call divergence(div,Q,work1,work2)


   ! sample loop:
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      u=Q(i,j,k,1)
      v=Q(i,j,k,2)
      w=Q(i,j,k,3)
   enddo
   enddo
   enddo

   n=2 

   ! compute the x derivative of n'th component of Q:   
   call der(Q(1,1,1,n),dx,dxx,work1,DX_AND_DXX,1)   

   ! compute the y derivative of n'th component of Q:   
   call der(Q(1,1,1,n),dx,dxx,work1,DX_AND_DXX,2)   

   ! compute the z derivative of n'th component of Q:   
   call der(Q(1,1,1,n),dx,dxx,work1,DX_AND_DXX,3)   


   ! compute an FFT:
  
   Qhat=Q
   do n=1,3
      call fft3d(Qhat(1,1,1,n),work1)
   enddo

   ! loop over FFT over component 'n'
   n=1
   do k=nz1,nz2
      km=(kmcord(k))
      do j=ny1,ny2
         jm=(jmcord(j))
         do i=nx1,nx2
            im=(imcord(i))
            ! print *,im,jm,km,Qhat(i,j,k,n)
         enddo
      enddo
   enddo
   


   time=time+tinc
   if (time>tstop) exit
enddo


call close_mpi
end program anal

