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
use sforcing
implicit none
real*8,save  :: Q(nx,ny,nz,n_var)
real*8,save  :: f(nx,ny,nz,n_var)
real*8,save  :: Q2ave(nx,ny,nz,n_var)
real*8,save  :: ufave(nx,ny,nz)
real*8,save  :: uxuxave(nx,ny,nz)
real*8,save  :: dxx(nx,ny,nz)
real*8,save  :: work1(nx,ny,nz)
real*8,save  :: work2(nx,ny,nz)
character(len=80) message,sdata
character(len=280) basename,fname
integer ierr,i,j,k,n,km,im,jm,icount,nave,n1,n2
real*8 :: tstart,tstop,tinc,time,time2
real*8 :: u,v,w,x,y
real*8 :: kr,ke,ck,xfac,f_diss

! input file
basename="iso12_250G_"
print *,basename
tstart=5
tstop=24
tinc=1.0
icount=0
call init_mpi       
call init_mpi_comm3d()
call init_model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if needed, initialize some constants.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Q=0
!
if (init_cond==3) call init_data_sht(Q,f,work1,work2)

nave=0
Q2ave=0
ufave=0
uxuxave=0
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
   nave=nave+1


   

   do n=1,3
   do k=nz1,nz2
      do j=ny1,ny2
         do i=nx1,nx2
            Q2ave(i,j,k,n)=Q2ave(i,j,k,n) + Q(i,j,k,n)**2
         enddo
      enddo
   enddo
   enddo


   do n1=1,3
   do n2=1,3
       ! Q(:,:,:,n1)*Q_n2(:,:,:,n1)
      call der(Q(1,1,1,n1),work2,dxx,work1,DX_AND_DXX,n2)   
      uxuxave(:,:,:)=uxuxave(:,:,:) + mu*Q(:,:,:,n1)*dxx(:,:,:)
   enddo
   enddo

   do i=1,n
      call fft3d(Q(1,1,1,n),work1)
   enddo
   f=0
   call sforce(f,Q,f_diss)
   do i=1,n
      call ifft3d(f(1,1,1,n),work1)
   enddo
   do i=1,n
      ufave=ufave+f(:,:,:,n)*Q(:,:,:,n)
   enddo


   
enddo
uxuxave=uxuxave/nave
Q2ave=Q2ave/nave
ufave=ufave/nave



fname = basename(1:len_trim(basename)) // ".urms"
call singlefile_io(time2,Q2ave(1,1,1,1),fname,work1,work2,0,io_pe)

fname = basename(1:len_trim(basename)) // ".vrms"
call singlefile_io(time2,Q2ave(1,1,1,2),fname,work1,work2,0,io_pe)

fname = basename(1:len_trim(basename)) // ".wrms"
call singlefile_io(time2,Q2ave(1,1,1,3),fname,work1,work2,0,io_pe)


fname = basename(1:len_trim(basename)) // ".udiss"
call singlefile_io(time2,uxuxave(1,1,1),fname,work1,work2,1,io_pe)

fname = basename(1:len_trim(basename)) // ".uf"
call singlefile_io(time2,ufave(1,1,1),fname,work1,work2,1,io_pe)


call close_mpi
end program anal


