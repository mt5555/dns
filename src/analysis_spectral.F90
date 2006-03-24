#include "macros.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! compute structure functions for many different directions
! in the periodic cube.
!
! also computes helical modes, distribution of helicity angle,  
! and spectra without helical modes
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
!  to compile and run:   make analysis_sk ; analysis_sk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program anal
use params
use mpi
use isoave
use pdf
use spectrum
use sforcing
implicit none

! set this to 0 to use 1 cpu version.  
! 1 cpu version uses significantly less memory then parallel version
! but doesn't compute some scalars like helicity dissipation rate.
integer :: always_use_parallel_code = 1;
integer :: use_serial

real*8,allocatable  :: Q(:,:,:,:)
real*8,allocatable  :: q1(:,:,:,:)
real*8,allocatable  :: q2(:,:,:,:)
real*8,allocatable  :: q3(:,:,:,:)
real*8,allocatable   :: work1(:,:,:)
real*8,allocatable   :: work2(:,:,:)
real*8,allocatable  :: work3(:,:,:)
real*8,allocatable  :: work4(:,:,:)

character(len=80) message,sdata,idata
character(len=280) basename,fname,tname
integer ierr,i,j,k,n,km,im,jm,icount
real*8 :: tstart,tstop,tinc,time,time2
real*8 :: u,v,w,x,y
real*8 :: kr,ke,ck,xfac,range(3,2),dummy,scale
integer :: lx1,lx2,ly1,ly2,lz1,lz2,nxlen,nylen,nzlen
integer :: nxdecomp,nydecomp,nzdecomp,csig,header_type
logical :: compute_cj,compute_scalar, compute_uvw,compute_pdfs
logical :: compute_hspec, compute_hfree
logical :: read_uvw
CPOINTER :: fid,fid1,fid2,fidcore


call init_mpi
call init_mpi_comm3d()
call init_model


! header_type of input data:
!    1              DNS code standard format
!    2              no headers             
!    3              Ensight headers
!    4              4 byte (fortran) header 
!

!cd
header_type=1; scale=1;
!header_type=4; scale=1/(2*pi)    ! for Takeshi's data
compute_pdfs=.false.
compute_cj=.false.
compute_scalar=.false.
compute_uvw=.false.
str_type=0
compute_hspec=.false.
read_uvw=.false.
compute_hfree=.true.		!extracting helicity-free modes

tstart=7.0
tstop=7.0
tinc=0.2
icount=0

nxdecomp=1
nydecomp=1
nzdecomp=1



! to read times from  file times.dat:
! tstart=-1; tinc=0; tname="times.dat"


! these lines are modifed by some sed scripts for automatic running
! of this code by putting in new values of tstart, tstop, tinc,
! nxdecomp,nydecomp,nzdecom, etc.
!SEDtstart
!SEDdecomp
!SEDcompcj
!SEDcompscalar
!SEDheadertype


if (scale/=1) then
   print *,'NOTE: scaling data by: ',scale
endif


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


if (use_serial==1) then
   allocate(Q(nx,ny,nz,ndim))
   allocate(work1(nx,ny,nz))
   allocate(work2(nx,ny,nz))
else
   ! parallel version requires extra data:
   allocate(Q(nx,ny,nz,ndim))
   allocate(q1(nx,ny,nz,ndim))
   allocate(q2(nx,ny,nz,ndim))
   allocate(work1(nx,ny,nz))
   allocate(work2(nx,ny,nz))
   if (nxdecomp*nydecomp*nzdecomp>1) then
      allocate(q3(nx,ny,nz,ndim))
      allocate(work1(nx,ny,nz))
      allocate(work2(nx,ny,nz))
      allocate(work3(nx,ny,nz))
      allocate(work4(nx,ny,nz))
   endif
endif

if (compute_cj .or. compute_hfree) then
   if (.not. allocated(q3))  allocate(q3(nx,ny,nz,ndim))
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if needed, initialize some constants.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Q=0

if (use_serial==1) then
   !  read in SK's data, and compute helical structure function
   !call convert_sk(Q,work1,work2);  stop;
endif

! run hfree test:  
! after code is debugged, remove this line and undefine TESTEXP in 
! spectrum.F90
!call compute_hfree_spec(Q,q1,q2,q3)
!stop



time=tstart
do
   icount=icount+1
   if (tstart<0) then
      ! read times from unit 83
      fname= rundir(1:len_trim(rundir)) // tname(1:len_trim(tname))
      if (icount==1)  open(83,file=fname)
      read(83,*,err=100,end=100) time
   endif
   

   if (compute_pdfs) then
      if (my_pe==io_pe) then
         write(sdata,'(f10.4)') 10000.0000 + time
         fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) // ".new.sf"
         print *,fname
         call copen(fname,"w",fid,ierr)
         if (ierr/=0) then
            write(message,'(a,i5)')"output_model(): Error opening .sf file errno=",ierr
            call abort(message)
         endif
         
         fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) // ".new.spdf"
         print *,fname
         call copen(fname,"w",fid2,ierr)
         if (ierr/=0) then
            write(message,'(a,i5)')"output_model(): Error opening .spdf file errno=",ierr
            call abort(message)
         endif

         fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) //".new.cores"
         print *,fname
         call copen(fname,"w",fidcore,ierr)
         if (ierr/=0) then
            write(message,'(a,i5)')&
            "output_model(): Error opening .cores file errno=",ierr
            call abort(message)
         endif
      endif
      
      if (.not. read_uvw) then
         call input_uvw(time,Q,q1,q2(1,1,1,1),q2(1,1,1,2),header_type)	
         Q=Q*scale;
	 read_uvw=.true.
      endif
      
!      uscale=.01 / pi2               ! (del_u)              units:  m/s
!      epsscale=.01 / pi2**(2./3.)    ! (mu*gradu**2)**1/3   units: 
      if (my_pe==io_pe) then
         print *,'PDF BINSIZE  (U,EPS):  ',uscale,epsscale
      endif
      call compute_all_pdfs(Q,q1)

      call output_pdf(time,fid,fid1,fid2,fidcore)
      if (my_pe==io_pe) call cclose(fid,ierr)
      if (my_pe==io_pe) call cclose(fid2,ierr)
      if (my_pe==io_pe) call cclose(fidcore,ierr)
   endif
   
   
   
   if (compute_hspec) then
      if (.not. read_uvw) then	
         if (use_serial==1) then
            stop 'compute_hspec needs q1,q2 allocated'
            call input_uvw(time,Q,dummy,work1,work2,header_type)
            Q=Q*scale;
         else
            call input_uvw(time,Q,q1,q2(1,1,1,1),q2(1,1,1,2),header_type)	
            Q=Q*scale;
         endif
         read_uvw=.true.	
      endif
      call compute_helicity_spectrum(Q,q2,q1,0)
      call output_helicity_spec(time,time)
   endif
   
   if (compute_hfree) then
      if (.not. read_uvw) then	
         if (use_serial==1) then
            stop 'compute_hfree needs q1,q2 allocated'
         else
            call input_uvw(time,Q,q1,q2(1,1,1,1),q2(1,1,1,2),header_type)
            Q=Q*scale;
         endif
         read_uvw=.true.	
      endif

#if 0
      q1=0
      ! compute u_y - v_x
      call der(Q(1,1,1,2),work1,dummy,work2,DX_ONLY,1)
      q1(:,:,:,3) = work1
      call der(Q(1,1,1,1),work1,dummy,work2,DX_ONLY,2)
      q1(:,:,:,3) = q1(:,:,:,3)-work1

      ! compute w_x - u_z
      call der(Q(1,1,1,1),work1,dummy,work2,DX_ONLY,3)
      q1(:,:,:,2) = work1
      call der(Q(1,1,1,3),work1,dummy,work2,DX_ONLY,1)
      q1(:,:,:,2) = q1(:,:,:,2)-work1

      ! compute v_z - w_y
      call der(Q(1,1,1,3),work1,dummy,work2,DX_ONLY,2)
      q1(:,:,:,1) = work1
      call der(Q(1,1,1,2),work1,dummy,work2,DX_ONLY,3)
      q1(:,:,:,1) = q1(:,:,:,1)-work1

      ! the vorcity is stored in q1, velocity in Q
      ! compute work1 = enstrophy
      ! compute work2 = helicity
      work1=0
      work2=0
      do k=nz1,nz2
      do j=ny1,ny2
      do i=nx1,nx2
         work1(i,j,k)=work1(i,j,k) + q1(i,j,k,1)**2 +  &
                                     q1(i,j,k,2)**2 +  &
                                     q1(i,j,k,3)**2  
         work2(i,j,k)=work2(i,j,k) + q1(i,j,k,1)*Q(i,j,k,1) + &
                                     q1(i,j,k,2)*Q(i,j,k,2) + &
                                     q1(i,j,k,3)*Q(i,j,k,3) 
      enddo
      enddo
      enddo

      ! now compute u - helicity * vorticity / enstrophy
      do k=nz1,nz2
      do j=ny1,ny2
      do i=nx1,nx2
         q2(i,j,k,1) = Q(i,j,k,1) - work2(i,j,k)*q1(i,j,k,1)/work1(i,j,k)
         q2(i,j,k,2) = Q(i,j,k,2) - work2(i,j,k)*q1(i,j,k,2)/work1(i,j,k)
         q2(i,j,k,3) = Q(i,j,k,3) - work2(i,j,k)*q1(i,j,k,3)/work1(i,j,k)
      enddo
      enddo
      enddo
      
      ! the remaining piece
      
         q1 = Q - q2
      
         
#endif

      call compute_hfree_spec(Q,q1,q2,q3)
      call output_hfree_spec(time,time)
   endif
   
   
   if (compute_uvw) then
      if (.not. read_uvw) then	
         if (use_serial==1) then
            call input_uvw(time,Q,dummy,work1,work2,header_type)
            Q=Q*scale;
         else
            call input_uvw(time,Q,q1,q2(1,1,1,1),q2(1,1,1,2),header_type)
            Q=Q*scale;
         endif
         read_uvw=.true.	
      endif
      
      do i=0,nxdecomp-1
         do j=0,nydecomp-1
            do k=0,nzdecomp-1  
               
               if (my_pe==io_pe) then
                  write(sdata,'(f10.4)') 10000.0000 + time
                  write(idata,'(i1)') str_type
                  if (str_type==0) then
                     fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) // ".new.isostr"
                  else
                     fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) // ".isostr" // idata(1:1)
                  endif
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
                     call isoavep(Q,q1,q1,q2,3,csig)
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
                     write(message,'(a,i5)') "output_model(): Error opening .isostr file errno=",ierr
                     call abort(message)
                  endif
                  call writeisoave(fid,time)
                  call cclose(fid,ierr)
               endif
               
            enddo
         enddo
      enddo
      
   endif
   


   
   if (compute_cj) then
      if (my_pe==io_pe) then
         write(sdata,'(f10.4)') 10000.0000 + time
         fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) // ".new.isow2s2"
         print *,fname
      endif
      
      if (.not. read_uvw) then
         call print_message("calling input_uvw")
         call input_uvw(time,Q,q1,q2(1,1,1,1),q2(1,1,1,2),header_type)	
         Q=Q*scale;
	 read_uvw=.true. ! dont set to .true.: we trash Q below:
      endif
      call print_message("calling compute_w2s2")
      call compute_w2s2(Q,q1,q2,q3)
      read_uvw=.false.  ! call to compute_w2s2 has trashed data in Q 
      call print_message("calling isoavep")
      call isoavep(Q,q1,q1,q2,2,csig)
      
      if (my_pe==io_pe) then
         call copen(fname,"w",fid,ierr)
         if (ierr/=0) then
            write(message,'(a,i5)') "output_model(): Error opening .isow2s2 file errno=",ierr
            call abort(message)
         endif
         call writeisoave_w2s2(fid,time)
         call cclose(fid,ierr)
      endif
   endif
   
   
   if (compute_scalar) then
      if (my_pe==io_pe) then
         write(sdata,'(f10.4)') 10000.0000 + time
         fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // sdata(2:10) // ".new.iso1"
         print *,fname
      endif
      
      call singlefile_io(time,Q(1,1,1,1),fname,q1(1,1,1,1),q1(1,1,1,2),1,io_pe)
      
      call isoavep(Q,q1,q1,q2,1,csig)
      if (my_pe==io_pe) then
         call copen(fname,"w",fid,ierr)
         if (ierr/=0) then
            write(message,'(a,i5)') "output_model(): Error opening .iso1 file errno=",ierr
            call abort(message)
         endif
         call writeisoave_scalar(fid,time)
         call cclose(fid,ierr)
      endif
   endif

   ! reset our flag, so we will read in the nxt data set
   read_uvw=.false.   

   if (tstart>=0) then   
      time=time+tinc
      if (io_pe==my_pe) print *,'time, tstart, tstop: ',time,tstart,tstop 
      if (time > max(tstop,tstart)+.005) exit
      if (time < min(tstop,tstart)-.005) exit
   endif
enddo
100 continue
call close_mpi
end program anal





subroutine compute_w2s2(Q,gradu,gradv,gradw)
use params
implicit none
real*8  :: Q(nx,ny,nz,3)
real*8  :: gradu(nx,ny,nz,3)
real*8  :: gradv(nx,ny,nz,3)
real*8  :: gradw(nx,ny,nz,3)

!local
integer :: i,j,k,n,m1,m2
real*8 :: vor(3),uij,uji,dummy

! compute vorticity and strain:  q1=gradu, q1=grad
do n=1,3
   call der(Q(1,1,1,1),gradu(1,1,1,n),dummy,gradw,DX_ONLY,n)
   call der(Q(1,1,1,2),gradv(1,1,1,n),dummy,gradw,DX_ONLY,n)
enddo
do n=1,3
   call der(Q(1,1,1,3),gradw(1,1,1,n),dummy,Q(1,1,1,1),DX_ONLY,n)
enddo

! 
! Q(:,:,:,1) = vor**2
! Q(:,:,:,2) = S**2
!
Q=0
do k=nz1,nz2
   do j=ny1,ny2
      do i=nx1,nx2
         vor(1)=gradw(i,j,k,2)-gradv(i,j,k,3)
         vor(2)=gradu(i,j,k,3)-gradw(i,j,k,1)
         vor(3)=gradv(i,j,k,1)-gradu(i,j,k,2)
         Q(i,j,k,1)=vor(1)**2+vor(2)**2+vor(3)**2
         
         do m1=1,3
            do m2=1,3
               if (m1==1) uij=gradu(i,j,k,m2)
               if (m1==2) uij=gradv(i,j,k,m2)
               if (m1==3) uij=gradw(i,j,k,m2)
               if (m2==1) uji=gradu(i,j,k,m1)
               if (m2==2) uji=gradv(i,j,k,m1)
               if (m2==3) uji=gradw(i,j,k,m1)
               !S(m1,m2)= .5*(uij+uji)
               Q(i,j,k,2) = Q(i,j,k,2) + ( .5*(uij+uji) ) **2     
            enddo
         enddo
      enddo
   enddo
enddo
end subroutine compute_w2s2





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
