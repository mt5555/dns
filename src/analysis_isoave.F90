#include "macros.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! compute structure functions for many different directions
! in the periodic cube.
! must be run on only 1 cpu.
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
program anal
use params
use mpi
use isoave
implicit none

! 1 cpu version uses less memory then parallel version
! but to test parallel version in serial, set this to 1.
integer :: always_use_parallel_code = 0;
integer :: use_serial

real*8,allocatable  :: Q(:,:,:,:)
real*8,allocatable  :: q1(:,:,:,:)
real*8,allocatable  :: q2(:,:,:,:)
real*8,allocatable  :: q3(:,:,:,:)
real*8,allocatable   :: work1(:,:,:)
real*8,allocatable   :: work2(:,:,:)
real*8,allocatable  :: work3(:,:,:)
real*8,allocatable  :: work4(:,:,:)

character(len=80) message,sdata
character(len=280) basename,fname
integer ierr,i,j,k,n,km,im,jm,icount
real*8 :: tstart,tstop,tinc,time,time2
real*8 :: u,v,w,x,y
real*8 :: kr,ke,ck,xfac,range(3,2)
integer :: lx1,lx2,ly1,ly2,lz1,lz2,nxlen,nylen,nzlen
integer :: nxdecomp,nydecomp,nzdecomp
CPOINTER :: fid




tstart=3.5
tstop=1.0
tinc=-.50
icount=0

nxdecomp=1
nydecomp=1
nzdecomp=1

call set_byteswap_input(1);
comp_sk_helical=.true.;  print *,'Also computing H_ltt...'



! these lines are modifed by some sed scripts for automatic running
! of this code by putting in new values of tstart, tstop, tinc,
! nxdecomp,nydecomp,nzdecom, etc.
!SEDbyteswap
!SEDtstart
!SEDdecomp



call init_mpi
call init_mpi_comm3d()
call init_model

if (ncpus==1) then
   use_serial=1;
else
   use_serial=0;
endif
if (always_use_parallel_code==1) use_serial=0;
if (use_serial==1) then
   call print_message("using serial version of isoave code")
else
   call print_message("using parallel version of isoave code")
endif


!call writepoints(); stop

allocate(Q(nx,ny,nz,ndim))
allocate(work1(nx,ny,nz))
allocate(work2(nx,ny,nz))

if (use_serial==0) then
   ! parallel version requires extra data:
   allocate(q1(nx,ny,nz,ndim))
   allocate(q2(nx,ny,nz,ndim))
   allocate(q3(nx,ny,nz,ndim))
   allocate(work3(nx,ny,nz))
   allocate(work4(nx,ny,nz))
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if needed, initialize some constants.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Q=0

!  read in SK's data, and compute helical structure function
!call convert_sk(Q,work1,work2);  stop;


time=tstart
do
   icount=icount+1

   call dataio(time,Q,work1,work2,1)

   do i=0,nxdecomp-1
   do j=0,nydecomp-1
   do k=0,nzdecomp-1  

      if (my_pe==io_pe) then
         write(sdata,'(f10.4)') 10000.0000 + time
         fname = runname(1:len_trim(runname)) // sdata(2:10) // ".isostr"
         if (nxdecomp*nydecomp*nzdecomp>1) then
            write(sdata,'(3i1)') i,j,k
            fname=fname(1:len_trim(fname)) // "_" // sdata(1:3)
         endif

         print *,fname
      endif

      if (use_serial==1) then
         nxlen=nslabx/nxdecomp
         nylen=nslaby/nydecomp
         nzlen=nslabz/nzdecomp
         
         lx1=nx1 + i*nxlen     
         ly1=ny1 + j*nylen     
         lz1=nz1 + k*nzlen     
         
         lx2=lx1 + nxlen-1
         ly2=ly1 + nylen-1
         lz2=lz1 + nzlen-1
         
         ! subcube version cannot run in parallel - it will abort
         call isoave1(Q,work1,work2,lx1,lx2,ly1,ly2,lz1,lz2)
      else
         if (nxdecomp*nydecomp*nzdecomp==1) then
            ! no subcubes:
            call isoavep(Q,q1,q2,q3)
         else
            range(1,1)=dble(i)/nxdecomp
            range(1,2)=dble(i+1)/nxdecomp
            range(2,1)=dble(j)/nxdecomp
            range(2,2)=dble(j+1)/nxdecomp
            range(3,1)=dble(k)/nxdecomp
            range(3,2)=dble(k+1)/nxdecomp
            call isoavep_subcube(Q,q1,q2,q3,range,work1,work2,work3,work4)
         endif
      endif
      


      if (my_pe==io_pe) then
         call copen(fname,"w",fid,ierr)
         if (ierr/=0) then
            write(message,'(a,i5)') "output_model(): Error opening .sf file errno=",ierr
            call abort(message)
         endif
         call writeisoave(fid,time)
         call cclose(fid,ierr)
      endif

   enddo
   enddo
   enddo


   time=time+tinc
   if (time > max(tstop,tstart)) exit
   if (time < min(tstop,tstart)) exit
enddo

call close_mpi
end program anal



subroutine dataio(time,Q,work1,work2,readflag)
use params
implicit none
real*8  :: Q(nx,ny,nz,3)
real*8  :: work1(nx,ny,nz)
real*8  :: work2(nx,ny,nz)
real*8  :: time
integer :: readflag   ! = 1 to read data, 0 to write data

real*8 time2
character(len=80) message,sdata
character(len=280) basename,fname

   time2=time
   write(sdata,'(f10.4)') 10000.0000 + time
   fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) // ".u"
   call print_message(fname(1:len_trim(fname)))
   call singlefile_io(time2,Q(1,1,1,1),fname,work1,work2,readflag,io_pe)

   fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) // ".v"
   call print_message(fname(1:len_trim(fname)))
   call singlefile_io(time2,Q(1,1,1,2),fname,work1,work2,readflag,io_pe)

   fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) // ".w"
   call print_message(fname(1:len_trim(fname)))
   call singlefile_io(time2,Q(1,1,1,3),fname,work1,work2,readflag,io_pe)

end subroutine





subroutine convert_sk(Q,work1,work2)
use params
implicit none
real*8  :: Q(nx,ny,nz,3)
real*8  :: work1(nx,ny,nz)
real*8  :: work2(nx,ny,nz)
real*8  :: time

character(len=80) message,sdata
character(len=280) basename,fname
integer :: N,ix,iy,iz

if (ncpu_x*ncpu_y*ncpu_z > 1) call abort("convert_sk must be run serial")
! read in data from alien file format, store in Q
open(unit = 10, form = 'unformatted', status = 'old', &
     file = '/home/scratch/taylorm/check256_hapiq_t0.8_velfield.out')
N=256
Q=0
time=0

print *,'reading in SK data'
read(10)(((Q(nx1+ix, ny1+iy, nz1+iz,1), &
           Q(nx1+ix, ny1+iy, nz1+iz,2), &
           Q(nx1+ix, ny1+iy, nz1+iz,3), &
           ix = 0, N-1), iy=0,N-1),iz=0,N-1)

Q=Q/(2*pi)
!
! and be sure to scale viscosity by 1/(2pi)**2
! Takashi's data: mu=.006 which scales to .00015198 
!
print *,'writing out DNS format data'
call dataio(time,Q,work1,work2,0)

end subroutine
