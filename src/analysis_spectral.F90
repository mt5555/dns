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
! compute structure functions for many different directions
! in the periodic cube.
!
! computes helical modes, distribution of helicity angle,  
! and spectra without helical modes 
!
! computes spectra of Craya-Herring projection modes (project_ch)
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
!  to compile and run:   make analysis_spectral ; analysis_spectral
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

real*8,allocatable  :: Q(:,:,:,:)
real*8,allocatable  :: Qhat(:,:,:,:)
real*8,allocatable  :: QI(:,:,:,:)
real*8,allocatable  :: QR(:,:,:,:)
real*8,allocatable  :: q1(:,:,:,:)
real*8,allocatable  :: q2(:,:,:,:)
real*8,allocatable  :: q3(:,:,:,:)
real*8,allocatable   :: work1(:,:,:)
real*8,allocatable   :: work2(:,:,:)

character(len=80) message,sdata,idata
character(len=280) basename,fname,tname
integer ierr,i,j,k,n,km,im,jm,icount
real*8 :: tstart,tstop,tinc,time,time2
real*8 :: u,v,w,x,y
real*8 :: kr,ke,ck,xfac,range(3,2),dummy,scale
integer :: lx1,lx2,ly1,ly2,lz1,lz2,nxlen,nylen,nzlen
integer :: nxdecomp,nydecomp,nzdecomp,csig,header_type
integer :: nints_e=16
real*8  :: ints_e(16)
logical :: compute_hspec, compute_hfree, compute_pv2spec, compute_pv2HA
logical :: compute_scalarsbous
logical :: read_uvw
logical :: project_ch
CPOINTER :: fid,fid1,fid2,fidcore,fid3


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
header_type=1; scale=1;           ! DNS standard data
!header_type=4; scale=1/(2*pi)    ! for Takeshi's data
compute_hspec=.false.
read_uvw=.false.
compute_hfree=.false.		!extracting helicity-free modes
project_ch=.true.         !Craya-Herring projection and spectra
compute_pv2spec = .false.  !potential enstrophy spectra .pv2spec,.normpvspec
compute_pv2HA = .false.    !compute Hussein Aluie's potential enstrophy spectra
compute_scalarsbous = .false. !compute .scalars-bous files

tstart=4.0
tstop=6.0
tinc=.1

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



!call writepoints(); stop


allocate(Q(nx,ny,nz,n_var))
allocate(q1(nx,ny,nz,n_var))
allocate(q2(nx,ny,nz,n_var))
allocate(work1(nx,ny,nz))
allocate(work2(nx,ny,nz))
if (nxdecomp*nydecomp*nzdecomp>1) then
   allocate(q3(nx,ny,nz,n_var))
endif
if (compute_hfree) then
   if (.not. allocated(q3))  allocate(q3(nx,ny,nz,ndim))
endif
if (compute_pv2spec) then
   if (.not. allocated(q3))  allocate(q3(nx,ny,nz,n_var))
endif
if (compute_scalarsbous) then
   if (.not. allocated(Qhat))  allocate(Qhat(g_nz2,nx_2dz,ny_2dz,n_var))
   if (.not. allocated(q3))  allocate(q3(nx,ny,nz,n_var))
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if needed, initialize some constants.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Q=0

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

   if (compute_hspec) then
      
      if (.not. read_uvw) then	
         call input_uvw(time,Q,q1,q2(1,1,1,1),q2(1,1,1,2),header_type)	
         Q=Q*scale;
         read_uvw=.true.	
      endif

      
      call compute_helicity_spectrum(Q,q2,q1,0)
      call output_helicity_spec(time,time)


      ! enable this block of code to recompute spectra.
      ! NOTE: this will output .spec, .pspec, .kspec, .hspec, .cospec, .spec2d
      ! and it will erase

      ! first read in the passive scalars too:
      call input_passive(runname,time,Q,work1,work2)

      call compute_spec(time,Q,q1,work1,work2)
      call output_spec(time,time)
      call compute_spec_2d(time,Q,q1,work1,work2)
      call output_2d_spec(time,time)


   endif
   
   if (compute_hfree) then
      if (.not. read_uvw) then	
         call input_uvw(time,Q,q1,q2(1,1,1,1),q2(1,1,1,2),header_type)
         Q=Q*scale;
         read_uvw=.true.	
      endif

      !this piece computes the non-helical part in physical space. Set #if 1 to run
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
   
   
   if (project_ch) then
      if (.not. read_uvw) then	
         call input_uvw(time,Q,q1,q2(1,1,1,1),q2(1,1,1,2),header_type)
         call input_passive(runname,time,Q,work1,work2)
         Q=Q*scale;
         read_uvw=.true.	
      endif
      if (.not. r_spec) then  ! r_spec reader will print stats, so we can skip this:
         call print_stats(Q,q1,work1,work2)
      endif

      call compute_project_ch(Q,q1,q2,work1,work2)
      call output_project_ch(time,time)
   endif

   
   if (compute_pv2spec) then
      if (.not. read_uvw) then	
         call input_uvw(time,Q,q1,q2(1,1,1,1),q2(1,1,1,2),header_type)
         call input_passive(runname,time,Q,work1,work2)
         Q=Q*scale;
         read_uvw=.true.	
      endif
      if (.not. r_spec) then  ! r_spec reader will print stats, so we can skip this:
         call print_stats(Q,q1,work1,work2)
      endif

      call compute_pv2_spec(time,Q,q1,q2,q3,work1,work2)
      call output_pv2_spec(time,time)
   endif

   
   if (compute_pv2HA) then
      if (.not. read_uvw) then	
         call input_uvw(time,Q,q1,q2(1,1,1,1),q2(1,1,1,2),header_type)
         call input_passive(runname,time,Q,work1,work2)
         Q=Q*scale;
         read_uvw=.true.	
      endif
      if (.not. r_spec) then  ! r_spec reader will print stats, so we can skip this:
         call print_stats(Q,q1,work1,work2)
      endif
      
      call compute_pv2_HA(Q,q1,work1,work2)
   endif
      
   if (compute_scalarsbous) then
      if (.not. read_uvw) then 
         call input_uvw(time,Q,q1,q2(1,1,1,1),q2(1,1,1,2),header_type)
         call input_passive(runname,time,Q,work1,work2)
         Q=Q*scale;
         read_uvw=.true.
      endif
      if (.not. r_spec) then  ! r_spec reader will print stats, so we can skip this:
         call print_stats(Q,q1,work1,work2)
      endif
      
      call compute_expensive_scalars_aspect(Q,Qhat,q1,q2,q3,work1,work2,nints_e,ints_e)


! output post-processing .scalars-bous file
      if (my_pe==io_pe) then
         write(message,'(f10.4)') 10000.0000 + time
         message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))//"-new" // message(2:10) // ".scalars-bous"
         call copen(message,"w",fid,ierr)
         write(6,*) "Opening scalars-bous file"
         if (ierr/=0) then
            write(message,'(a,i5)') "diag_output(): Error opening new .scalars-bous file errno=",ierr
            call abortdns(message)
         endif
         x=nints_e; call cwrite8(fid,x,1)
         call cwrite8(fid,time,1)
         call cwrite8(fid,ints_e,nints_e)
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

if (ncpu_x*ncpu_y*ncpu_z > 1) call abortdns("convert_sk must be run serial")
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


