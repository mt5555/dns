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
real*8,save  :: Qhat(nx,ny,nz,n_var)
real*8,save  :: vor(nx,ny,nz,n_var)
real*8,save  :: vor2d(nx,ny,nz),temp(nx,ny,nz)
real*8,save  :: vorsave(nx,ny,nz)
real*8,save  :: dx(nx,ny,nz)
real*8,save  :: dxx(nx,ny,nz)
real*8,save  :: div(nx,ny,nz)
real*8,save  :: work1(nx,ny,nz)
real*8,save  :: work2(nx,ny,nz)
real*8,save  :: Qhatfilt(nx,ny,nz,n_var),Qfilt(nx,ny,nz,n_var)
real*8,save  :: vorfilt(nx,ny,nz),pvfilt(nx,ny,nz)
character(len=80) message,sdata
character(len=280) basename,fname,fnamewrite
integer ierr,i,j,k,n,km,im,jm,icount,ictof,ivorsave,ntrunc
integer icentroid
real*8 :: tstart,tstop,tinc,time,time2
real*8 :: u,v,w,x,y
real*8 :: kr,ke,ck,xfac

! input file
basename="Ia0"
print *,basename
tstart=10.5
tstop=10.5
tinc=.1
icentroid=0 ! Compute the KE centroid? 0 = no, 1 = yes
icount=0  ! A counter for averaging to compute time averages
ictof = 0 ! If 1 then flip endian sign
ivorsave=0 ! If 1 then compute an averaged vorticity
ntrunc = 32 ! If not zero find Q with wave numbers ntrunc and above
call init_mpi       
call init_mpi_comm3d()
call init_model

! byte swap the input data:
call set_byteswap_input(1);




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

   fname = basename(1:len_trim(basename)) // sdata(2:10) // ".h"
   call singlefile_io(time2,Q(1,1,1,3),fname,work1,work2,1,io_pe)
   if (ictof == 1) then
     fnamewrite = basename(1:len_trim(basename)) // sdata(2:10) // ".Fh"
     open(42,file=fnamewrite,status='unknown',action='write',form='unformatted')
     write(42) Q(nx1:nx2,ny1:ny2,nz1:nz2,1)
   endif


   ! compute vorticity 
   call vorticity(vor,Q,work1,work2)
   ! for 2D problems, 1st and 2nd components of vor() will be zero
   vor2d(:,:,:)=(vor(:,:,:,3) + fcor)/(1.d0+Q(:,:,:,3))
   if (ivorsave /=0) then
     print *,'Fo1 -> ',time, vor2d(10,10,1),vorsave(10,10,1)
     vorsave(:,:,:) = vorsave(:,:,:) + vor2d(:,:,:)
     print *,'Fo2 -> ',time, vorsave(10,10,1)
   endif

   if (ictof == 1) then
     fnamewrite = basename(1:len_trim(basename)) // sdata(2:10) // ".Fvor"
     open(42,file=fnamewrite,status='unknown',action='write',form='unformatted')
     write(42) vor(nx1:nx2,ny1:ny2,nz1:nz2,3)
     fnamewrite = basename(1:len_trim(basename)) // sdata(2:10) // ".Fpvor"
     open(42,file=fnamewrite,status='unknown',action='write',form='unformatted')
     write(42) vor2d(nx1:nx2,ny1:ny2,nz1:nz2)
   endif
!#if 0
   ! compute divergence
   !call divergence(div,Q,work1,work2)


   ! sample loop:
!   do k=nz1,nz2
!   do j=ny1,ny2
!   do i=nx1,nx2
!      u=Q(i,j,k,1)
!      v=Q(i,j,k,2)
!      w=Q(i,j,k,3)
!   enddo
!   enddo
!   enddo

!   n=2 

   ! compute the x derivative of n'th component of Q:   
!   call der(Q(1,1,1,n),dx,dxx,work1,DX_AND_DXX,1)   

   ! compute the y derivative of n'th component of Q:   
!   call der(Q(1,1,1,n),dx,dxx,work1,DX_AND_DXX,2)   

   ! compute the z derivative of n'th component of Q:   
!   call der(Q(1,1,1,n),dx,dxx,work1,DX_AND_DXX,3)   
!#endif




   if (ntrunc > 0) then
   ! compute an FFT:
  
   Qhat=Q
   do n=1,3
      call fft3d(Qhat(1,1,1,n),work1)
   enddo

   ! compute the grid points values, but with the truncated series 'n'
   ! First store the truncated spectrual coeffs in Qhatfilt
   Qhatfilt = 0
   do k=nz1,nz2
      km=(kmcord(k))
      do j=ny1,ny2
         jm=(jmcord(j))
         do i=nx1,nx2
            im=(imcord(i))

            xfac = 2*2*2
            if (km==0) xfac=xfac/2
            if (jm==0) xfac=xfac/2
            if (im==0) xfac=xfac/2

!            print *,'im jm for filter -> ',i,j,im,jm
            if ((abs(im) <= ntrunc) .and. (abs(jm) <= ntrunc)) then
!               print *,i,j,im,jm,Qhat(i,j,k,3)
               Qhatfilt(i,j,k,1) = Qhat(i,j,k,1)           
               Qhatfilt(i,j,k,2) = Qhat(i,j,k,2)           
               Qhatfilt(i,j,k,3) = Qhat(i,j,k,3)
             else
!               print *,i,j,im,jm,'zero'
               Qhatfilt(i,j,k,1) = 0.d0           
               Qhatfilt(i,j,k,2) = 0.d0           
               Qhatfilt(i,j,k,3) = 0.d0           
            endif
         enddo
      enddo
   enddo
   ! Compute the filtered values of the fields 
   Qfilt=Qhatfilt
   do n=1,3
      call ifft3d(Qfilt(1,1,1,n),work1)
   enddo
   ! Now compute the vorticity
   call vorticity(vor,Qfilt,work1,work2)
   ! for 2D problems, 1st and 2nd components of vor() will be zero
   vorfilt(:,:,:)=vor(:,:,:,3)
   pvfilt(:,:,:)=(vor(:,:,:,3) + fcor)/(1.d0+Qfilt(:,:,:,3))
   print *,'Filtered Vorticity -> ',time, vorfilt(10,10,1)
 
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

   endif

   if (icentroid > 0) then
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
close(55)
if (ivorsave /= 0) then
!  Close the ke centroid file
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
endif

call close_mpi
end program anal


