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
real*8,save  :: vorsave(nx,ny,nz)
real*8,save  :: dx(nx,ny,nz)
real*8,save  :: dxx(nx,ny,nz)
real*8,save  :: div(nx,ny,nz)
real*8,save  :: work1(nx,ny,nz)
real*8,save  :: work2(nx,ny,nz)
character(len=80) message,sdata
character(len=280) basename,fname
integer ierr,i,j,k,n,km,im,jm,icount
real*8 :: tstart,tstop,tinc,time,time2
real*8 :: u,v,w,x,y
real*8 :: kr,ke,ck,xfac

! input file
basename="I_a0"
print *,basename
tstart=0
tstop=11.3
tinc=.10
icount=0
vorsave=0
call init_mpi       
call init_mpi_comm3d()
call init_model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if needed, initialize some constants.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Q=0
!  Open this file to store the Kinetic Energy Centroid!
   open(55,file='KeCentroid.gp')
!
if (init_cond==3) call init_data_sht(Q,vor,work1,work2)

time=tstart
do
icount=icount+1
print *,'icount = ',icount
   write(sdata,'(f10.4)') 10000.0000 + time
   fname = basename(1:len_trim(basename)) // sdata(2:10) // ".u"
   if (my_pe==io_pe) then
      print *,'filename: ',fname(1:len_trim(fname))
   endif
   call singlefile_io(time2,Q(1,1,1,1),fname,work1,work2,1,io_pe)

   fname = basename(1:len_trim(basename)) // sdata(2:10) // ".v"
   call singlefile_io(time2,Q(1,1,1,2),fname,work1,work2,1,io_pe)

   fname = basename(1:len_trim(basename)) // sdata(2:10) // ".w"
   call singlefile_io(time2,Q(1,1,1,3),fname,work1,work2,1,io_pe)


#if 0
   ! compute vorticity 
   call vorticity(vor,Q,work1,work2)
   ! for 2D problems, 1st and 2nd components of vor() will be zero
   vor2d(:,:,:)=(vor(:,:,:,3) + fcor)/(1.d0+Q(:,:,:,3))
   print *,'Fo1 -> ',time, vor2d(10,10,1),vorsave(10,10,1)
   vorsave(:,:,:) = vorsave(:,:,:) + vor2d(:,:,:)
   print *,'Fo2 -> ',time, vorsave(10,10,1)

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
!   call der(Q(1,1,1,n),dx,dxx,work1,DX_AND_DXX,1)   

   ! compute the y derivative of n'th component of Q:   
!   call der(Q(1,1,1,n),dx,dxx,work1,DX_AND_DXX,2)   

   ! compute the z derivative of n'th component of Q:   
!   call der(Q(1,1,1,n),dx,dxx,work1,DX_AND_DXX,3)   
#endif




   ! compute an FFT:
  
   Qhat=Q
   do n=1,3
      call fft3d(Qhat(1,1,1,n),work1)
   enddo


   ! loop over FFT over component 'n'
   ck=0
   ke=0
   do k=nz1,nz2
      km=(kmcord(k))
      do j=ny1,ny2
         jm=(jmcord(j))
         do i=nx1,nx2
            im=(imcord(i))
            ! print *,im,jm,km,Qhat(i,j,k,n)

            xfac = 2*2*2
            if (km==0) xfac=xfac/2
            if (jm==0) xfac=xfac/2
            if (im==0) xfac=xfac/2


!            kr = sqrt(real(im**2) + real(jm**2) + real(km**2))
            kr =  sqrt(real(im**2) + real(jm**2))
            ck=ck + xfac*kr*(Qhat(i,j,k,1)**2 + Qhat(i,j,k,2)**2)
            ke=ke + xfac*(Qhat(i,j,k,1)**2 + Qhat(i,j,k,2)**2)
         enddo
      enddo
   enddo
#ifdef USE_MPI
!   xfac=ke
!   call MPI_allreduce(xfac,ke,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
!   xfac=ck
!  call MPI_allreduce(xfac,ck,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
   if (my_pe==io_pe) then
      print *,'ck = ',ck
      print *,'ke = ',ke
      print *, 'centroid = ',ck/ke
      write(55,'(e14.6,2x,e14.6)') real(time),ck/ke
   endif
  
   ke=0
   k=1
   do j=ny1,ny2
   do i=nx1,nx2
      u=Q(i,j,k,1)
      v=Q(i,j,k,2)
      ke=ke+u**2+v**2
   enddo
   enddo
	print *,'grid ke',ke/g_nx/g_ny




   time=time+tinc
   if (time>tstop) exit
enddo
!  Close the ke centroid file
   close(55)

!  Write out the averaged PV
!   delx = real(1.d0/double(nx2-nx1))
!   dely = real(1.d0/double(ny2-ny1))
   open(55,file='AveragePV.gp')
   do k=nz1,nz2
      km=(kmcord(k))
      y = 0 - dely
      do j=ny1,ny2
         jm=(jmcord(j))
         y = y + dely
         x = 0 - delx
         do i=nx1,nx2
            im=(imcord(i))
            x = x + delx
            write(55,'(e14.6,2x,e14.6,2x,e14.6)') real(x),real(y),real(vorsave(i,j,k))/real(icount)
         enddo
         write(55,*)
      enddo
   enddo
!   do i=1,int(n)
!      write(55,'(i4,e14.6)') i-1,spec(i)
!   enddo
   close(55)

call close_mpi
end program anal


