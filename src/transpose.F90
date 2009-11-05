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

#undef A2AOVERLAP

module transpose
use params
use mpi
implicit none



integer :: io_nodes(0:ncpu_z),nio,inc,comm_io
real*8,private :: tmx1,tmx2
real*8,private,allocatable :: sendbuf(:),recbuf(:)

integer           :: mpi_maxio=-1
character(len=4)  :: mpi_stripe="64"
character(len=12) :: mpi_stride="8388608"
contains

subroutine transpose_init

#ifdef USE_MPI
integer max_size

! compute send/rec buffer size for x,y and z transpose
max_size=1
if (ncpu_z>1) then
   max_size=max(max_size,nx_2dz*nslabz*ny_2dz)
endif
if (ncpu_x>1) then
   max_size=max(max_size,nslabx*nslabz*ny_2dx)
endif
if (ncpu_y>1) then
   max_size=max(max_size,nslaby*nslabz*nx_2dy)
endif

allocate (sendbuf(max_size))
allocate (recbuf(max_size))



#endif
end subroutine




subroutine mpi_io_init(init)
!
! set up all the I/O processors
! and create communicator, comm_io
!
! init=1   first time, create comm_io
! init=0   not first time - free comm_io, create a new comm_io
!
integer init

integer i,dest_pe3(3),key,color,ierr

#ifdef USE_MPI_IO
if (init==0) then
   ! release the old communicator, construct a new one:
   call mpi_comm_free(comm_io,ierr)
endif

if (mpi_maxio<0) then
   if (g_nx <= 512) then
      mpi_maxio=4
      mpi_stripe="4"
   else if (g_nx <= 1024) then
      mpi_maxio=8
      mpi_stripe="8"
   else
      ! defaults:  
      mpi_maxio=8
      mpi_stripe="8"
#ifdef OSF1
      mpi_maxio=32
      mpi_stripe="32"
#endif
   endif
else
   ! calling program set mpi_maxio, so dont change any values
endif


nio=1
nio=min(mpi_maxio,ncpu_z)

! no more than 50% of the CPUS should do I/O
! on Q, this should be 25%, but sometimes we run 64 cpus on 32 nodes.
if (ncpu_x*ncpu_y*ncpu_z>4) nio=min(nio,(ncpu_x*ncpu_y*ncpu_z)/4)

inc=ncpu_z/nio
if (io_pe==my_pe) then
   print *,'number of mpi_io cpus: ',nio,do_mpi_io
   print *,'output_size (bytes): ',output_size
   print *,'input_size  (bytes): ',input_size
endif

if (nio>1) then
   do i=0,ncpu_z-1
      dest_pe3(1)=0
      dest_pe3(2)=0
      dest_pe3(3)=inc*(i/inc)
      call cart_rank(comm_3d,dest_pe3,io_nodes(i),ierr)
   enddo
   call mpi_barrier(comm_3d,ierr)
!print *,'my_pe=',my_pe,' my_z=',my_z, 'my io_pe: ',io_nodes(my_z)
   call mpi_barrier(comm_3d,ierr)
   if (io_nodes(my_z)<0) call abortdns("error in MPI-IO decomposition")
else
   io_nodes=0
endif


! am i an io pe?
key=0
color=0
if (io_nodes(my_z)==my_pe) color=1
key=0

! everyone with color=1 joins a new group, comm_sforcing
call mpi_comm_split(comm_3d,color,key,comm_io,ierr);
!if (color==0) call mpi_comm_free(comm_io,ierr)

#endif
end subroutine



#undef PBLOCK
#ifdef PBLOCK
  PBLOCK code deleted - see old versions in CVS if interested
#else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! transform data in a cartesian decomposition into a 2D decomposition
! with z as the leading index.
!
! input: p
! ouput: pt, and the dimensions of pt: n1,n1d,n2,n2d,n3,n3d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_to_z(p,pt,n1,n1d,n2,n2d,n3,n3d)
!use params
real*8 p(nx,ny,nz)
real*8 pt(g_nz2,nx_2dz,ny_2dz)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc,iproc2,splitx,iproc_x,iproc_y
integer i,j,k,jj,l,ii
#ifdef USE_MPI
integer ierr,dest_pe,request(2),statuses(MPI_STATUS_SIZE,2)
integer dest_pe3(3),tag
#endif

call wallclock(tmx1)

!
! tranpose from reference decompostion to z-pencil decomposition
!
! reference decompostion:  processors:   ncpu_x x ncpu_y x ncpu_z
! grid size (each process)               nslabx x nslaby x nslabz 
!
! In the z-pencil domain, the parallel decomposition is:
!    ncpu_x  x  (ncpu_y*ncpu_z)   x  1
!
!  the array is:  p(g_nz2,nx_2dz,ny_2dz)
!  where p(k,i,j)   i = x direction index
!                   j = y direction index
!                   k = z direction index
!  
! Consider the case with ncpu_z = N  (mesh size in Z-direction) 
! Then y direction ny_2dz=1 and we cant do Fourier computations!
!
! So we need to modify so that the decomposition is:
!     2*ncpu_x  x  (ncpu_y*ncpu_z/2)   x  1
!
! Original algorithm:
! each cube (which depending on ncpu_x,y and z, may actually be a pencil or 
! slab) is broken, along the y axis, into ncpu_z slabs of
! size ny_2dz = (ny2-ny1+1)/ncpu_z
!
! Now we want to modify this so that each cube is broken along the y axis
! into (ncpu_z/2) slabs, and along the x axis into 2.
! so it sends half as much x data and twice as much y data 
! 
!  
!  
!

n1=g_nz
n1d=g_nz2   	
n2=nx_2dz
n2d=nx_2dz
n3=ny_2dz
n3d=ny_2dz

splitx=nslabx/nx_2dz  

do iproc2=0,ncpu_z-1  ! loop over each slab
#ifdef A2AOVERLAP
   iproc=mod(ncpu_z+iproc2-my_z,ncpu_z)
#else
   iproc=iproc2
#endif
   ! each cube is broken into ncpu_z pieces, each of which 
   ! is sent to a different process
   iproc_y = iproc/splitx        ! split y direction into this many pieces
   iproc_x = mod(iproc,splitx)   ! split x direction into this many pieces
   

   if (iproc==my_z) then
      do j=1,ny_2dz  ! loop over points in a single slab
         jj=ny1 + iproc_y*ny_2dz +j -1
         !ASSERT("transpose_to_z jj failure 1",jj<=ny2)
         ASSERT("transpose_to_z jj failure 2",jj>=ny1)
         do k=nz1,nz2
         do i=1,nx_2dz
            ii=nx1 + iproc_x*nx_2dz +i -1
            if (jj>ny2) then
               pt(k+iproc*nslabz-nz1+1,i,j)=0
            else
               pt(k+iproc*nslabz-nz1+1,i,j)=p(ii,jj,k)
            endif
         enddo
#if 0
         do i=nx1,nx2
            if (jj>ny2) then
               pt(k+iproc*nslabz-nz1+1,i-nx1+1,j)=0
            else
               pt(k+iproc*nslabz-nz1+1,i-nx1+1,j)=p(i,jj,k)
            endif
         enddo
#endif
         enddo
      enddo
   else
#ifdef USE_MPI
      dest_pe3(1)=my_x
      dest_pe3(2)=my_y
      dest_pe3(3)=iproc
      call cart_rank(comm_3d,dest_pe3,dest_pe,ierr)

      l=ny_2dz*nslabz*nx_2dz
      tag=my_z
      call mpi_irecv(recbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(1),ierr)
      ASSERT("transpose_to_z: MPI_IRecv failure 1",ierr==0)


      l=0
      do j=1,ny_2dz  ! loop over points in a single slab
         jj=ny1 + iproc_y*ny_2dz +j -1
         !ASSERT("transpose_to_z jj failure 1",jj<=ny2)
         ASSERT("transpose_to_z jj failure 2",jj>=ny1)
         do k=nz1,nz2
         do i=1,nx_2dz
            ii=nx1 + iproc_x*nx_2dz +i -1
            l=l+1
            if (jj>ny2) then
               sendbuf(l)=0
            else
               sendbuf(l)=p(ii,jj,k)
            endif
         enddo
         enddo
      enddo

!     send/rec buffer to (my_x,my_y,iproc)
      tag=iproc
      call mpi_isend(sendbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(2),ierr)
      ASSERT("transpose_to_z: MPI_ISend failure 1",ierr==0)
      call mpi_waitall(2,request,statuses,ierr) 	
      ASSERT("transpose_to_z: MPI_waitalll failure 1",ierr==0)

      l=0
      do j=1,ny_2dz  ! loop over points in a single slab
         do k=nz1,nz2
         do i=1,nx_2dz
	    l=l+1
            pt(k+iproc*nslabz-nz1+1,i,j)=recbuf(l)
         enddo
         enddo
      enddo
#endif
   endif
enddo

call wallclock(tmx2) 
tims(6)=tims(6)+(tmx2-tmx1)          
end subroutine


#endif





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! transform data in a cartesian decomposition into a 2D decomposition
! with z as the leading index.
!
! input: p
! ouput: pt, and the dimensions of pt: n1,n1d,n2,n2d,n3,n3d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_to_z_from_refyxz(p,pt,n1,n1d,n2,n2d,n3,n3d)
!use params
real*8 p(ny,nx,nz)
real*8 pt(g_nz2,nx_2dz,ny_2dz)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc,iproc2,splitx,iproc_x,iproc_y
integer i,j,k,jj,l,ii
#ifdef USE_MPI
integer ierr,dest_pe,request(2),statuses(MPI_STATUS_SIZE,2)
integer dest_pe3(3),tag
#endif

call wallclock(tmx1)

!
! tranpose from reference decompostion to z-pencil decomposition
!
! reference decompostion:  processors:   ncpu_x x ncpu_y x ncpu_z
! grid size (each process)               nslabx x nslaby x nslabz 
!
! In the z-pencil domain, the parallel decomposition is:
!    ncpu_x  x  (ncpu_y*ncpu_z)   x  1
!
!  the array is:  p(g_nz2,nx_2dz,ny_2dz)
!  where p(k,i,j)   i = x direction index
!                   j = y direction index
!                   k = z direction index
!  
! Consider the case with ncpu_z = N  (mesh size in Z-direction) 
! Then y direction ny_2dz=1 and we cant do Fourier computations!
!
! So we need to modify so that the decomposition is:
!     2*ncpu_x  x  (ncpu_y*ncpu_z/2)   x  1
!
! Original algorithm:
! each cube (which depending on ncpu_x,y and z, may actually be a pencil or 
! slab) is broken, along the y axis, into ncpu_z slabs of
! size ny_2dz = (ny2-ny1+1)/ncpu_z
!
! Now we want to modify this so that each cube is broken along the y axis
! into (ncpu_z/2) slabs, and along the x axis into 2.
! so it sends half as much x data and twice as much y data 
! 
!  
!  
!

n1=g_nz
n1d=g_nz2   	
n2=nx_2dz
n2d=nx_2dz
n3=ny_2dz
n3d=ny_2dz

splitx=nslabx/nx_2dz  

do iproc2=0,ncpu_z-1  ! loop over each slab
#ifdef A2AOVERLAP
   iproc=mod(ncpu_z+iproc2-my_z,ncpu_z)
#else
   iproc=iproc2
#endif
   ! each cube is broken into ncpu_z pieces, each of which 
   ! is sent to a different process
   iproc_y = iproc/splitx        ! split y direction into this many pieces
   iproc_x = mod(iproc,splitx)   ! split x direction into this many pieces
   

   if (iproc==my_z) then
      do j=1,ny_2dz  ! loop over points in a single slab
         jj=ny1 + iproc_y*ny_2dz +j -1
         !ASSERT("transpose_to_z jj failure 1",jj<=ny2)
         ASSERT("transpose_to_z jj failure 2",jj>=ny1)
         do k=nz1,nz2
         do i=1,nx_2dz
            ii=nx1 + iproc_x*nx_2dz +i -1
            if (jj>ny2) then
               pt(k+iproc*nslabz-nz1+1,i,j)=0
            else
               pt(k+iproc*nslabz-nz1+1,i,j)=p(jj,ii,k)
            endif
         enddo
         enddo
      enddo
   else
#ifdef USE_MPI
      dest_pe3(1)=my_x
      dest_pe3(2)=my_y
      dest_pe3(3)=iproc
      call cart_rank(comm_3d,dest_pe3,dest_pe,ierr)

      l=ny_2dz*nslabz*nx_2dz
      tag=my_z
      call mpi_irecv(recbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(1),ierr)
      ASSERT("transpose_to_z: MPI_IRecv failure 1",ierr==0)


      l=0
      do j=1,ny_2dz  ! loop over points in a single slab
         jj=ny1 + iproc_y*ny_2dz +j -1
         !ASSERT("transpose_to_z jj failure 1",jj<=ny2)
         ASSERT("transpose_to_z jj failure 2",jj>=ny1)
         do k=nz1,nz2
         do i=1,nx_2dz
            ii=nx1 + iproc_x*nx_2dz +i -1
            l=l+1
            if (jj>ny2) then
               sendbuf(l)=0
            else
               sendbuf(l)=p(jj,ii,k)
            endif
         enddo
         enddo
      enddo

!     send/rec buffer to (my_x,my_y,iproc)
      tag=iproc
      call mpi_isend(sendbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(2),ierr)
      ASSERT("transpose_to_z: MPI_ISend failure 1",ierr==0)
      call mpi_waitall(2,request,statuses,ierr) 	
      ASSERT("transpose_to_z: MPI_waitalll failure 1",ierr==0)

      l=0
      do j=1,ny_2dz  ! loop over points in a single slab
         do k=nz1,nz2
         do i=1,nx_2dz
	    l=l+1
            pt(k+iproc*nslabz-nz1+1,i,j)=recbuf(l)
         enddo
         enddo
      enddo
#endif
   endif
enddo

call wallclock(tmx2) 
tims(20)=tims(20)+(tmx2-tmx1)          
end subroutine












!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! transform data to a cartesian decomposition from a 2D decomposition
! with z as the leading index.
!
! input: pt, and the dimensions of pt: n1,n1d,n2,n2d,n3,n3d
! ouput: p 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_from_z(pt,p,n1,n1d,n2,n2d,n3,n3d)
!use params
real*8 p(nx,ny,nz)
real*8 pt(g_nz2,nx_2dz,ny_2dz)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc,iproc2,splitx,iproc_x,iproc_y
integer i,j,k,jj,l,ii
#ifdef USE_MPI
integer ierr,dest_pe,request(2),statuses(MPI_STATUS_SIZE,2)
integer dest_pe3(3),tag
#endif

!
! each cube is broken, along the y axis, into ncpu_z slabs of
! size ny_2dz = (ny2-ny1+1)/ncpu_z
!
! in the z direction, the dimension is nslabz = (nz2-nz1+1)

call wallclock(tmx1)

! If any of these fail, then pt was probably not computed
! via a call to transpose_to_z().
ASSERT("transpose_from_z dimension failure 2",n1==g_nz)
ASSERT("transpose_from_z dimension failure 3",n1d==g_nz2)
ASSERT("transpose_from_z dimension failure 4",n2==nx_2dz)
ASSERT("transpose_from_z dimension failure 5",n2d==nx_2dz)
ASSERT("transpose_from_z dimension failure 6",n3==ny_2dz)
ASSERT("transpose_from_z dimension failure 7",n3d==ny_2dz)

splitx=nslabx/nx_2dz  
do iproc2=0,ncpu_z-1  ! loop over each slab
#ifdef A2AOVERLAP
   iproc=mod(ncpu_z+iproc2-my_z,ncpu_z)
#else
   iproc=iproc2
#endif
   iproc_y = iproc/splitx        ! split y direction into this many pieces
   iproc_x = mod(iproc,splitx)   ! split x direction into this many pieces


   if (iproc==my_z) then
      do j=1,ny_2dz  ! loop over points in a single slab
         jj=ny1 + iproc_y*ny_2dz +j -1
         !ASSERT("transpose_from_z jj failure 1",jj<=ny2)
         ASSERT("transpose_from_z jj failure 2",jj>=ny1)
         if (jj<=ny2) then
         do k=nz1,nz2
         do i=1,nx_2dz
            ii=nx1 + iproc_x*nx_2dz +i -1
            p(ii,jj,k)=pt(k+iproc*nslabz-nz1+1,i,j)
         enddo
         enddo
         endif
      enddo
   else

#ifdef USE_MPI
      dest_pe3(1)=my_x
      dest_pe3(2)=my_y
      dest_pe3(3)=iproc
      call cart_rank(comm_3d,dest_pe3,dest_pe,ierr)

      l=ny_2dz*nslabz*nx_2dz
      tag=my_z
      call mpi_irecv(recbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(1),ierr)
      ASSERT("transpose_from_z: MPI_IRecv failure 1",ierr==0)

      l=0
      do j=1,ny_2dz  ! loop over points in a single slab
         do k=nz1,nz2
         do i=1,nx_2dz
            l=l+1
            sendbuf(l)=pt(k+iproc*nslabz-nz1+1,i,j)
         enddo
         enddo
      enddo

!     send/rec
      tag=iproc
      call mpi_isend(sendbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(2),ierr)
      ASSERT("transpose_from_z: MPI_ISend failure 1",ierr==0)
      call mpi_waitall(2,request,statuses,ierr) 	
      ASSERT("transpose_from_z: MPI_waitalll failure 1",ierr==0)


      l=0
      do j=1,ny_2dz  ! loop over points in a single slab
         jj=ny1 + iproc_y*ny_2dz +j -1
         !ASSERT("transpose_from_z jj failure 1",jj<=ny2)
         ASSERT("transpose_from_z jj failure 2",jj>=ny1)
         do k=nz1,nz2
         do i=1,nx_2dz
            ii=nx1 + iproc_x*nx_2dz +i -1
	    l=l+1
            if (jj<=ny2) p(ii,jj,k)=recbuf(l)
         enddo
         enddo
      enddo
#endif
   endif
enddo
call wallclock(tmx2) 
tims(7)=tims(7)+(tmx2-tmx1)          

end subroutine








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! transform data to a cartesian decomposition from a 2D decomposition
! with z as the leading index.
!
! input: pt, and the dimensions of pt: n1,n1d,n2,n2d,n3,n3d
! ouput: p 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_from_z_to_refyxz(pt,p,n1,n1d,n2,n2d,n3,n3d)
!use params
real*8 p(ny,nx,nz)
real*8 pt(g_nz2,nx_2dz,ny_2dz)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc,iproc2,splitx,iproc_x,iproc_y
integer i,j,k,jj,l,ii
#ifdef USE_MPI
integer ierr,dest_pe,request(2),statuses(MPI_STATUS_SIZE,2)
integer dest_pe3(3),tag
#endif

!
! each cube is broken, along the y axis, into ncpu_z slabs of
! size ny_2dz = (ny2-ny1+1)/ncpu_z
!
! in the z direction, the dimension is nslabz = (nz2-nz1+1)

call wallclock(tmx1)

! If any of these fail, then pt was probably not computed
! via a call to transpose_to_z().
ASSERT("transpose_from_z dimension failure 2",n1==g_nz)
ASSERT("transpose_from_z dimension failure 3",n1d==g_nz2)
ASSERT("transpose_from_z dimension failure 4",n2==nx_2dz)
ASSERT("transpose_from_z dimension failure 5",n2d==nx_2dz)
ASSERT("transpose_from_z dimension failure 6",n3==ny_2dz)
ASSERT("transpose_from_z dimension failure 7",n3d==ny_2dz)

splitx=nslabx/nx_2dz  
do iproc2=0,ncpu_z-1  ! loop over each slab
#ifdef A2AOVERLAP
   iproc=mod(ncpu_z+iproc2-my_z,ncpu_z)
#else
   iproc=iproc2
#endif
   iproc_y = iproc/splitx        ! split y direction into this many pieces
   iproc_x = mod(iproc,splitx)   ! split x direction into this many pieces


   if (iproc==my_z) then
      do j=1,ny_2dz  ! loop over points in a single slab
         jj=ny1 + iproc_y*ny_2dz +j -1
         !ASSERT("transpose_from_z jj failure 1",jj<=ny2)
         ASSERT("transpose_from_z jj failure 2",jj>=ny1)
         if (jj<=ny2) then
         do k=nz1,nz2
         do i=1,nx_2dz
            ii=nx1 + iproc_x*nx_2dz +i -1
            p(jj,ii,k)=pt(k+iproc*nslabz-nz1+1,i,j)
         enddo
         enddo
         endif
      enddo
   else

#ifdef USE_MPI
      dest_pe3(1)=my_x
      dest_pe3(2)=my_y
      dest_pe3(3)=iproc
      call cart_rank(comm_3d,dest_pe3,dest_pe,ierr)

      l=ny_2dz*nslabz*nx_2dz
      tag=my_z
      call mpi_irecv(recbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(1),ierr)
      ASSERT("transpose_from_z: MPI_IRecv failure 1",ierr==0)

      l=0
      do j=1,ny_2dz  ! loop over points in a single slab
         do k=nz1,nz2
         do i=1,nx_2dz
            l=l+1
            sendbuf(l)=pt(k+iproc*nslabz-nz1+1,i,j)
         enddo
         enddo
      enddo

!     send/rec
      tag=iproc
      call mpi_isend(sendbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(2),ierr)
      ASSERT("transpose_from_z: MPI_ISend failure 1",ierr==0)
      call mpi_waitall(2,request,statuses,ierr) 	
      ASSERT("transpose_from_z: MPI_waitalll failure 1",ierr==0)


      l=0
      do j=1,ny_2dz  ! loop over points in a single slab
         jj=ny1 + iproc_y*ny_2dz +j -1
         !ASSERT("transpose_from_z jj failure 1",jj<=ny2)
         ASSERT("transpose_from_z jj failure 2",jj>=ny1)
         do k=nz1,nz2
         do i=1,nx_2dz
            ii=nx1 + iproc_x*nx_2dz +i -1
	    l=l+1
            if (jj<=ny2) p(jj,ii,k)=recbuf(l)
         enddo
         enddo
      enddo
#endif
   endif
enddo
call wallclock(tmx2) 
tims(21)=tims(21)+(tmx2-tmx1)          

end subroutine






























!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! transform data from a cartesian decomposition to a 2D decomposition
! with x as the leading index.  
! slabs are chopped up in the y direction.
!
! input: pt, and the dimensions of pt: n1,n1d,n2,n2d,n3,n3d
! ouput: p 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_to_x(p,pt,n1,n1d,n2,n2d,n3,n3d)
!use params
real*8 p(nx,ny,nz)
real*8 pt(g_nx2,nslabz,ny_2dx)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc,iproc2
integer i,j,k,jj,l
#ifdef USE_MPI
integer ierr,dest_pe,request(2),statuses(MPI_STATUS_SIZE,2)
integer dest_pe3(3),tag
#endif

call wallclock(tmx1)

!
! each cube is broken, along the z axis, into ncpu_x slabs of
! size ny_2dx
!
! in the x direction, the dimension is nslabx = (nx2-nx1+1)


! If any of these fail, then pt was probably not computed
! via a call to transpose_to_*().
n1=g_nx
n1d=g_nx2
n2=nslabz
n2d=nslabz
n3=ny_2dx
n3d=ny_2dx


do iproc2=0,ncpu_x-1  ! loop over each slab
#ifdef A2AOVERLAP
   iproc=mod(ncpu_x+iproc2-my_x,ncpu_x)
#else
   iproc=iproc2
#endif


   if (iproc==my_x) then
      do j=1,ny_2dx  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2dx +j -1
         !ASSERT("transpose_to_x2 jj failure 1",jj<=ny2)
         ASSERT("transpose_to_x2 jj failure 2",jj>=ny1)
         do i=nx1,nx2
         do k=nz1,nz2
            if (jj<=ny2) then
               pt(i+iproc*nslabx-nx1+1,k-nz1+1,j)=p(i,jj,k)
            else
               pt(i+iproc*nslabx-nx1+1,k-nz1+1,j)=0
            endif
         enddo
         enddo
      enddo
   else
#ifdef USE_MPI
      dest_pe3(1)=iproc
      dest_pe3(2)=my_y
      dest_pe3(3)=my_z
      call cart_rank(comm_3d,dest_pe3,dest_pe,ierr)

      l=ny_2dx*nslabx*nslabz      
      tag=my_x
      call mpi_irecv(recbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(1),ierr)
      ASSERT("transpose_to_x2: MPI_IRecv failure 1",ierr==0)

      l=0
      do j=1,ny_2dx  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2dx +j -1
         !ASSERT("transpose_to_x2 jj failure 1",jj<=ny2)
         ASSERT("transpose_to_x2 jj failure 2",jj>=ny1)
         do i=nx1,nx2
         do k=nz1,nz2
            l=l+1
            if (jj<=ny2) then
               sendbuf(l)=p(i,jj,k)
            else
               sendbuf(l)=0
            endif
         enddo
         enddo
      enddo

!     send/rec buffer 
      tag=iproc
      call mpi_isend(sendbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(2),ierr)
      ASSERT("transpose_to_x2: MPI_ISend failure 1",ierr==0)
      call mpi_waitall(2,request,statuses,ierr) 	
      ASSERT("transpose_to_x2: MPI_waitalll failure 1",ierr==0)


      l=0
      do j=1,ny_2dx  ! loop over points in a single slab
         do i=nx1,nx2
         do k=nz1,nz2
            l=l+1
            pt(i+iproc*nslabx-nx1+1,k-nz1+1,j)=recbuf(l)
         enddo
         enddo
      enddo
#endif
   endif
enddo

call wallclock(tmx2) 
tims(8)=tims(8)+(tmx2-tmx1)          
end subroutine






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! transform data from a cartesian decomposition to a 2D decomposition
! with x as the leading index.  
! slabs are chopped up in the y direction.
!
! input: pt, and the dimensions of pt: n1,n1d,n2,n2d,n3,n3d
! ouput: p 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_to_x_from_refyxz(p,pt,n1,n1d,n2,n2d,n3,n3d)
!use params
real*8 p(ny,nx,nz)
real*8 pt(g_nx2,nslabz,ny_2dx)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc,iproc2
integer i,j,k,jj,l
#ifdef USE_MPI
integer ierr,dest_pe,request(2),statuses(MPI_STATUS_SIZE,2)
integer dest_pe3(3),tag
#endif

call wallclock(tmx1)

!
! each cube is broken, along the z axis, into ncpu_x slabs of
! size ny_2dx
!
! in the x direction, the dimension is nslabx = (nx2-nx1+1)


! If any of these fail, then pt was probably not computed
! via a call to transpose_to_*().
n1=g_nx
n1d=g_nx2
n2=nslabz
n2d=nslabz
n3=ny_2dx
n3d=ny_2dx


do iproc2=0,ncpu_x-1  ! loop over each slab
#ifdef A2AOVERLAP
   iproc=mod(ncpu_x+iproc2-my_x,ncpu_x)
#else
   iproc=iproc2
#endif


   if (iproc==my_x) then
      do j=1,ny_2dx  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2dx +j -1
         !ASSERT("transpose_to_x2 jj failure 1",jj<=ny2)
         ASSERT("transpose_to_x2 jj failure 2",jj>=ny1)
         do i=nx1,nx2
         do k=nz1,nz2
            if (jj<=ny2) then
               pt(i+iproc*nslabx-nx1+1,k-nz1+1,j)=p(jj,i,k)
            else
               pt(i+iproc*nslabx-nx1+1,k-nz1+1,j)=0
            endif
         enddo
         enddo
      enddo
   else
#ifdef USE_MPI
      dest_pe3(1)=iproc
      dest_pe3(2)=my_y
      dest_pe3(3)=my_z
      call cart_rank(comm_3d,dest_pe3,dest_pe,ierr)

      l=ny_2dx*nslabx*nslabz      
      tag=my_x
      call mpi_irecv(recbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(1),ierr)
      ASSERT("transpose_to_x2: MPI_IRecv failure 1",ierr==0)

      l=0
      do j=1,ny_2dx  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2dx +j -1
         !ASSERT("transpose_to_x2 jj failure 1",jj<=ny2)
         ASSERT("transpose_to_x2 jj failure 2",jj>=ny1)
         do i=nx1,nx2
         do k=nz1,nz2
            l=l+1
            if (jj<=ny2) then
               sendbuf(l)=p(jj,i,k)
            else
               sendbuf(l)=0
            endif
         enddo
         enddo
      enddo

!     send/rec buffer 
      tag=iproc
      call mpi_isend(sendbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(2),ierr)
      ASSERT("transpose_to_x2: MPI_ISend failure 1",ierr==0)
      call mpi_waitall(2,request,statuses,ierr) 	
      ASSERT("transpose_to_x2: MPI_waitalll failure 1",ierr==0)


      l=0
      do j=1,ny_2dx  ! loop over points in a single slab
         do i=nx1,nx2
         do k=nz1,nz2
            l=l+1
            pt(i+iproc*nslabx-nx1+1,k-nz1+1,j)=recbuf(l)
         enddo
         enddo
      enddo
#endif
   endif
enddo

call wallclock(tmx2) 
tims(22)=tims(22)+(tmx2-tmx1)          
end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! transform data to a cartesian decomposition from a 2D decomposition
! with x as the leading index.  same as transpose_from_x, but
! slabs are chopped up in the y direction istead of z direction.
!
! input: pt, and the dimensions of pt: n1,n1d,n2,n2d,n3,n3d
! ouput: p 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_from_x(pt,p,n1,n1d,n2,n2d,n3,n3d)
!use params

real*8 p(nx,ny,nz)
real*8 pt(g_nx2,nslabz,ny_2dx)
integer n1,n1d,n2,n2d,n3,n3d


!local variables
integer iproc,iproc2
integer i,j,k,jj,l
#ifdef USE_MPI
integer ierr,dest_pe,request(2),statuses(MPI_STATUS_SIZE,2)
integer dest_pe3(3),tag
#endif

call wallclock(tmx1)

!
! each cube is broken, along the z axis, into ncpu_x slabs of
! size ny_2dx
!
! in the x direction, the dimension is nslabx = (nx2-nx1+1)


! If any of these fail, then pt was probably not computed
! via a call to transpose_to_*().
ASSERT("transpose_from_x2 dimension failure 2",n1==g_nx)
ASSERT("transpose_from_x2 dimension failure 3",n1d==g_nx2)
ASSERT("transpose_from_x2 dimension failure 4",n2==nslabz)
ASSERT("transpose_from_x2 dimension failure 5",n2d==nslabz)
ASSERT("transpose_from_x2 dimension failure 6",n3==ny_2dx)
ASSERT("transpose_from_x2 dimension failure 7",n3d==ny_2dx)


do iproc2=0,ncpu_x-1  ! loop over each slab
#ifdef A2AOVERLAP
   iproc=mod(ncpu_x+iproc2-my_x,ncpu_x)
#else
   iproc=iproc2
#endif

   if (iproc==my_x) then
      do j=1,ny_2dx  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2dx +j -1
         !ASSERT("transpose_from_x2 jj failure 1",jj<=ny2)
         ASSERT("transpose_from_x2 jj failure 2",jj>=ny1)
         if (jj<=ny2) then
         do i=nx1,nx2
         do k=nz1,nz2
            p(i,jj,k)=pt(i+iproc*nslabx-nx1+1,k-nz1+1,j)
         enddo
         enddo
         endif
      enddo
   else
#ifdef USE_MPI
      dest_pe3(1)=iproc
      dest_pe3(2)=my_y
      dest_pe3(3)=my_z
      call cart_rank(comm_3d,dest_pe3,dest_pe,ierr)

      l=ny_2dx*nslabx*nslabz 
      tag=my_x
      call mpi_irecv(recbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(1),ierr)
      ASSERT("transpose_from_x2: MPI_IRecv failure 1",ierr==0)

      l=0
      do j=1,ny_2dx  ! loop over points in a single slab
         do i=nx1,nx2
         do k=nz1,nz2
            l=l+1
            sendbuf(l)=pt(i+iproc*nslabx-nx1+1,k-nz1+1,j)
         enddo
         enddo
      enddo

!     send/rec buffer 
      tag=iproc
      call mpi_isend(sendbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(2),ierr)
      ASSERT("transpose_from_x2: MPI_ISend failure 1",ierr==0)
      call mpi_waitall(2,request,statuses,ierr) 	
      ASSERT("transpose_from_x2: MPI_waitalll failure 1",ierr==0)


      l=0
      do j=1,ny_2dx  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2dx +j -1
         !ASSERT("transpose_from_x2 jj failure 1",jj<=ny2)
         ASSERT("transpose_from_x2 jj failure 2",jj>=ny1)
         do i=nx1,nx2
         do k=nz1,nz2
            l=l+1
            if (jj<=ny2) p(i,jj,k)=recbuf(l)
         enddo
         enddo
      enddo
#endif
   endif
enddo

call wallclock(tmx2) 
tims(9)=tims(9)+(tmx2-tmx1)          
end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! transform data to a cartesian decomposition from a 2D decomposition
! with x as the leading index.  same as transpose_from_x, but
! slabs are chopped up in the y direction istead of z direction.
!
! input: pt, and the dimensions of pt: n1,n1d,n2,n2d,n3,n3d
! ouput: p 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_from_x_to_refyxz(pt,p,n1,n1d,n2,n2d,n3,n3d)
!use params

real*8 p(ny,nx,nz)
real*8 pt(g_nx2,nslabz,ny_2dx)
integer n1,n1d,n2,n2d,n3,n3d


!local variables
integer iproc,iproc2
integer i,j,k,jj,l
#ifdef USE_MPI
integer ierr,dest_pe,request(2),statuses(MPI_STATUS_SIZE,2)
integer dest_pe3(3),tag
#endif

call wallclock(tmx1)

!
! each cube is broken, along the z axis, into ncpu_x slabs of
! size ny_2dx
!
! in the x direction, the dimension is nslabx = (nx2-nx1+1)


! If any of these fail, then pt was probably not computed
! via a call to transpose_to_*().
ASSERT("transpose_from_x2 dimension failure 2",n1==g_nx)
ASSERT("transpose_from_x2 dimension failure 3",n1d==g_nx2)
ASSERT("transpose_from_x2 dimension failure 4",n2==nslabz)
ASSERT("transpose_from_x2 dimension failure 5",n2d==nslabz)
ASSERT("transpose_from_x2 dimension failure 6",n3==ny_2dx)
ASSERT("transpose_from_x2 dimension failure 7",n3d==ny_2dx)


do iproc2=0,ncpu_x-1  ! loop over each slab
#ifdef A2AOVERLAP
   iproc=mod(ncpu_x+iproc2-my_x,ncpu_x)
#else
   iproc=iproc2
#endif

   if (iproc==my_x) then
      do j=1,ny_2dx  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2dx +j -1
         !ASSERT("transpose_from_x2 jj failure 1",jj<=ny2)
         ASSERT("transpose_from_x2 jj failure 2",jj>=ny1)
         if (jj<=ny2) then
         do i=nx1,nx2
         do k=nz1,nz2
            p(jj,i,k)=pt(i+iproc*nslabx-nx1+1,k-nz1+1,j)
         enddo
         enddo
         endif
      enddo
   else
#ifdef USE_MPI
      dest_pe3(1)=iproc
      dest_pe3(2)=my_y
      dest_pe3(3)=my_z
      call cart_rank(comm_3d,dest_pe3,dest_pe,ierr)

      l=ny_2dx*nslabx*nslabz 
      tag=my_x
      call mpi_irecv(recbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(1),ierr)
      ASSERT("transpose_from_x2: MPI_IRecv failure 1",ierr==0)

      l=0
      do j=1,ny_2dx  ! loop over points in a single slab
         do i=nx1,nx2
         do k=nz1,nz2
            l=l+1
            sendbuf(l)=pt(i+iproc*nslabx-nx1+1,k-nz1+1,j)
         enddo
         enddo
      enddo

!     send/rec buffer 
      tag=iproc
      call mpi_isend(sendbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(2),ierr)
      ASSERT("transpose_from_x2: MPI_ISend failure 1",ierr==0)
      call mpi_waitall(2,request,statuses,ierr) 	
      ASSERT("transpose_from_x2: MPI_waitalll failure 1",ierr==0)


      l=0
      do j=1,ny_2dx  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2dx +j -1
         !ASSERT("transpose_from_x2 jj failure 1",jj<=ny2)
         ASSERT("transpose_from_x2 jj failure 2",jj>=ny1)
         do i=nx1,nx2
         do k=nz1,nz2
            l=l+1
            if (jj<=ny2) p(jj,i,k)=recbuf(l)
         enddo
         enddo
      enddo
#endif
   endif
enddo

call wallclock(tmx2) 
tims(23)=tims(23)+(tmx2-tmx1)          
end subroutine


















!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! transform data in a cartesian decomposition into a 2D decomposition
! with y as the leading index.
!
! input: p
! ouput: pt, and the dimensions of pt: n1,n1d,n2,n2d,n3,n3d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_to_y(p,pt,n1,n1d,n2,n2d,n3,n3d)
!use params

real*8 p(nx,ny,nz)
real*8 pt(g_ny2,nslabz,nx_2dy)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc,iproc2
integer i,j,k,ii,l
#ifdef USE_MPI
!real*8 sendbuf(nslaby*nslabz*nx_2dy)
!real*8 recbuf(nslaby*nslabz*nx_2dy)
integer ierr,dest_pe,request(2),statuses(MPI_STATUS_SIZE,2)
integer dest_pe3(3),tag
#endif

call wallclock(tmx1)


!
! each cube is broken, along the x axis, into ncpu_y slabs of
! size nx_2dy
!
! in the y direction, the dimension is nslaby = (ny2-ny1+1)


n1=g_ny
n1d=g_ny2
n2=nslabz
n2d=nslabz
n3=nx_2dy
n3d=nx_2dy



do iproc2=0,ncpu_y-1  ! loop over each slab
#ifdef A2AOVERLAP
   iproc=mod(ncpu_y+iproc2-my_y,ncpu_y)
#else
   iproc=iproc2
#endif


   if (iproc==my_y) then
      do i=1,nx_2dy  ! loop over points in a single slab
         ii=nx1 + iproc*nx_2dy +i -1
         !ASSERT("transpose_to_y ii failure 1",ii<=nx2)
         ASSERT("transpose_to_y ii failure 2",ii>=nx1)
         do j=ny1,ny2
         do k=nz1,nz2
            if (ii>nx2) then
               pt(j+iproc*nslaby-ny1+1,k-nz1+1,i)=0
            else
               pt(j+iproc*nslaby-ny1+1,k-nz1+1,i)=p(ii,j,k)
            endif
         enddo
         enddo
      enddo
   else
#ifdef USE_MPI
      dest_pe3(1)=my_x
      dest_pe3(2)=iproc
      dest_pe3(3)=my_z
      call cart_rank(comm_3d,dest_pe3,dest_pe,ierr)

      l=nx_2dy*nslaby*nslabz
      tag=my_y
      call mpi_irecv(recbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(1),ierr)
      ASSERT("transpose_to_y: MPI_IRecv failure 1",ierr==0)


      l=0
      do i=1,nx_2dy  ! loop over points in a single slab
         ii=nx1 + iproc*nx_2dy +i -1
         !ASSERT("transpose_to_y ii failure 1",ii<=nx2)
         ASSERT("transpose_to_y ii failure 2",ii>=nx1)
         do j=ny1,ny2
         do k=nz1,nz2
            l=l+1
            if (ii>nx2) then
               sendbuf(l)=0
            else
               sendbuf(l)=p(ii,j,k)
            endif
         enddo
         enddo
      enddo

!     send/rec buffer 
      tag=iproc
      call mpi_isend(sendbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(2),ierr)
      ASSERT("transpose_to_y: MPI_ISend failure 1",ierr==0)
      call mpi_waitall(2,request,statuses,ierr) 	
      ASSERT("transpose_to_y: MPI_waitalll failure 1",ierr==0)

      l=0
      do i=1,nx_2dy  ! loop over points in a single slab
         do j=ny1,ny2
         do k=nz1,nz2
            l=l+1
            pt(j+iproc*nslaby-ny1+1,k-nz1+1,i)=recbuf(l)
         enddo
         enddo
      enddo
#endif
   endif
enddo

call wallclock(tmx2) 
tims(10)=tims(10)+(tmx2-tmx1)          
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! transform data in to a cartesian decomposition from a 2D decomposition
! with y as the leading index.
!
! input: pt and its dimensions: n1,n1d,n2,n2d,n3,n3d
! ouput: p
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_from_y(pt,p,n1,n1d,n2,n2d,n3,n3d)
!use params

real*8 p(nx,ny,nz)
real*8 pt(g_ny2,nslabz,nx_2dy)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc,iproc2
integer i,j,k,ii,l
#ifdef USE_MPI
!real*8 sendbuf(nslaby*nslabz*nx_2dy)
!real*8 recbuf(nslaby*nslabz*nx_2dy)
integer ierr,dest_pe,request(2),statuses(MPI_STATUS_SIZE,2)
integer dest_pe3(3),tag
#endif

call wallclock(tmx1)

!
! each cube is broken, along the x axis, into ncpu_y slabs of
! size nx_2dy
!
! in the y direction, the dimension is nslaby = (ny2-ny1+1)




! If any of these fail, then pt was probably not computed
! via a call to transpose_to_y().
ASSERT("transpose_from_y dimension failure 2",n1==g_ny)
ASSERT("transpose_from_y dimension failure 3",n1d==g_ny2)
ASSERT("transpose_from_y dimension failure 4",n2==nslabz)
ASSERT("transpose_from_y dimension failure 5",n2d==nslabz)
ASSERT("transpose_from_y dimension failure 6",n3==nx_2dy)
ASSERT("transpose_from_y dimension failure 7",n3d==nx_2dy)



do iproc2=0,ncpu_y-1  ! loop over each slab
#ifdef A2AOVERLAP
   iproc=mod(ncpu_y+iproc2-my_y,ncpu_y)
#else
   iproc=iproc2
#endif
   if (iproc==my_y) then
      do i=1,nx_2dy  ! loop over points in a single slab
         ii=nx1 + iproc*nx_2dy +i-1
         !ASSERT("transpose_from_y ii failure 1",ii<=nx2)
         ASSERT("transpose_from_y ii failure 2",ii>=nx1)
         if (ii<=nx2) then
         do j=ny1,ny2
         do k=nz1,nz2
            p(ii,j,k)=pt(j+iproc*nslaby-ny1+1,k-nz1+1,i)
         enddo
         enddo
         endif
      enddo
   else
#ifdef USE_MPI
      dest_pe3(1)=my_x
      dest_pe3(2)=iproc
      dest_pe3(3)=my_z
      call cart_rank(comm_3d,dest_pe3,dest_pe,ierr)

      l=nx_2dy*nslaby*nslabz
      tag=my_y
      call mpi_irecv(recbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(1),ierr)
      ASSERT("transpose_from_y: MPI_IRecv failure 1",ierr==0)

      l=0
      do i=1,nx_2dy  ! loop over points in a single slab
         do j=ny1,ny2
         do k=nz1,nz2
            l=l+1
            sendbuf(l)=pt(j+iproc*nslaby-ny1+1,k-nz1+1,i)
         enddo
         enddo
      enddo

!     send/rec buffer 
      tag=iproc
      call mpi_isend(sendbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(2),ierr)
      ASSERT("transpose_from_y: MPI_ISend failure 1",ierr==0)
      call mpi_waitall(2,request,statuses,ierr) 	
      ASSERT("transpose_from_y: MPI_waitalll failure 1",ierr==0)


      l=0
      do i=1,nx_2dy  ! loop over points in a single slab
         ii=nx1 + iproc*nx_2dy + i -1
         !ASSERT("transpose_from_y ii failure 1",ii<=nx2)
         ASSERT("transpose_from_y ii failure 2",ii>=nx1)
         do j=ny1,ny2
         do k=nz1,nz2
            l=l+1
            if (ii<=nx2) p(ii,j,k)=recbuf(l)
         enddo
         enddo
      enddo
#endif
   endif
enddo

call wallclock(tmx2) 
tims(11)=tims(11)+(tmx2-tmx1)          
end subroutine









subroutine output1(p,pt,buf,fid,fpe_main,offset)
!use params
!use mpi
CPOINTER fid
integer :: fpe_main             ! cpu to do the I/O
real*8 :: p(nx,ny,nz)
real*8 :: pt(g_nx2,nslabz,ny_2dx)
real*8 :: buf(o_nx,ny_2dx)

! local vars
real*8 saved_edge(o_nx)
integer :: sending_pe,ierr,tag,z_pe,y_pe,x_pe
integer i,j,k,l,extra_k,kuse,dest_pe3(3),fpe
integer n1,n1d,n2,n2d,n3,n3d
integer :: ny_2dx_actual ,offset
logical :: first_seek
#ifdef USE_MPI
integer request,statuses(MPI_STATUS_SIZE)
#endif
#ifdef USE_MPI_IO
integer (KIND=MPI_OFFSET_KIND) ::  zpos,ypos
#endif

ny_2dx_actual = ny_2dx

first_seek=.true.


call transpose_to_x(p,pt,n1,n1d,n2,n2d,n3,n3d)
ASSERT("output1 dimension failure 2",n1==g_nx)
ASSERT("output1 dimension failure 3",n1d==g_nx2)
ASSERT("output1 dimension failure 4",n2==nslabz)
ASSERT("output1 dimension failure 5",n2d==nslabz)
ASSERT("output1 dimension failure 6",n3==ny_2dx)
ASSERT("output1 dimension failure 7",n3d==ny_2dx)

if (o_nz>g_nz .and. g_bdy_z1/=PERIODIC) then
   call abortdns("output1: cannot handle o_nz=g_nz+1 with non-periodic b.c.")
endif

do z_pe=0,ncpu_z-1
extra_k=0
fpe=fpe_main
if (do_mpi_io) fpe=io_nodes(z_pe)
if (z_pe==ncpu_z-1 .and. o_nz>g_nz) extra_k=1
do k=1,nslabz+extra_k
do y_pe=0,ncpu_y-1
do x_pe=0,ncpu_x-1


   ! for non-perfect load balanced cases, the last few columns
   ! of data will not be real data.
   ! Find the largest value jx<=ny_2dx such that the corresponding
   ! value of j in the 3D decompositoin is i<=ny2
   ! set ny_2dx_actual = jx. 
   ! formula:  
   ! 
   ! ny2 == j == -1 + ny1 + jx + my_x*ny_2dx    
   ! nslaby = jx + my_x*ny_2dx
   ! jx = nslaby - my_x*ny_2dx



   ! output pt(1:g_nx,k,1:ny_2dx) from cpus: x_pe,y_pe,z_pe
   l=o_nx*ny_2dx_actual

   if (k>nslabz) then 
      ! this is the periodic z-direction case. refetch data at z=0
      kuse=1
      dest_pe3(1)=x_pe
      dest_pe3(2)=y_pe
      dest_pe3(3)=0
   else
      kuse=k
      dest_pe3(1)=x_pe
      dest_pe3(2)=y_pe
      dest_pe3(3)=z_pe
   endif

#ifdef USE_MPI
   tag=1
   call cart_rank(comm_3d,dest_pe3,sending_pe,ierr)
#else
   sending_pe=my_pe
#endif

   if (sending_pe==my_pe) then

      buf(1:g_nx,:)=pt(1:g_nx,kuse,:)
      if (my_pe == fpe) then
         ! dont send message to self
      else
#ifdef USE_MPI
         tag=1
         call mpi_isend(buf,l,MPI_REAL8,fpe,tag,comm_3d,request,ierr)
         ASSERT("output1: MPI_ISend failure",ierr==0)
         call mpi_waitall(1,request,statuses,ierr) 	
         ASSERT("output1: MPI_waitalll failure",ierr==0)
#endif
      endif
   endif

   if (my_pe==fpe) then
      if (sending_pe==my_pe) then
         ! dont recieve message from self
      else
#ifdef USE_MPI
         call mpi_irecv(buf,l,MPI_REAL8,sending_pe,tag,comm_3d,request,ierr)
         ASSERT("output1: MPI_IRecv failure",ierr==0)
         call mpi_waitall(1,request,statuses,ierr) 	
         ASSERT("output1: MPI_waitalll failure",ierr==0)
#endif
      endif

      !
      !  we now have the slab of data at position:
      !  zpos = z_pe*nslabz + k 
      !  ypos = y_pe*nslaby + x_pe*ny_2dx_actual
      !
      !  so the offset is:  zpos*(o_nx*o_ny) + ypos*o_nx
      !
      if (o_nx>g_nx) then
         if (g_bdy_x1==PERIODIC) then
            buf(o_nx,:)=buf(1,:)  ! append to the end, x-direction
         else
            buf(o_nx,:)=0
         endif
      endif
#ifdef USE_MPI_IO
      if (do_mpi_io) then
         zpos = z_pe*nslabz + k-1 
         ypos = y_pe*nslaby + x_pe*ny_2dx_actual
         zpos = output_size*(offset + zpos*o_nx*o_ny+ypos*o_nx)
         if (first_seek) then
            call mpi_file_seek(fid,zpos,MPI_SEEK_SET,ierr)
            first_seek=.false.
         endif
      endif
#endif
      call mwrite8(fid,buf,o_nx*ny_2dx_actual)


      if (o_ny>g_ny) then
         if (y_pe==0 .and. x_pe==0) then
            if (g_bdy_y1==PERIODIC) then
               ! save this edge to append to the (y_pe=ncpu_y,x_pe-ncpu_x) edge
               saved_edge=buf(:,1)
            else
               saved_edge=0
            endif
         endif
         if (y_pe==ncpu_y-1 .and. x_pe==ncpu_x-1) then     ! append to the end, y-direction
            if (do_mpi_io) then
#ifdef USE_MPI_IO
!              no longer needed - we are automatically at the right location
!              zpos = z_pe*nslabz + k -1
!              ypos = y_pe*nslaby + (1+x_pe)*ny_2dx_actual
!              zpos = output_size*(offset + zpos*o_nx*o_ny+ypos*o_nx)
!              call mpi_file_seek(fid,zpos,MPI_SEEK_SET,ierr)
!              call mpi_file_write(fid,saved_edge,o_nx,MPI_REAL8,statuses,ierr)
#endif
            endif
            call mwrite8(fid,saved_edge,o_nx)
         endif
      endif

   endif
enddo
enddo
enddo
enddo

end subroutine



!
!
! dealias_nx = number of coefficients = 2+2*im_max
! if im_max = g_nx/2-1               full spectrum
! for im_max < g_nx/2:               all wave numbers up to
!                                    and including im_max
!
!  for no dealiasing:     im_max = g_nx/2-1
!  for 2/3 dealias rule,  im_max = g_nx/3
!
!integer,parameter :: dealias_nx = 2+2*(g_nx/3) = 2 + 2*im_max
!integer,parameter :: dealias_ny = 2+2*(g_ny/3) = 2 + 2*jm_max
!integer,parameter :: dealias_nz = 2+2*(g_nz/3) = 2 + 2*km_max
!
!
!  coefficients are stored:  0 6 1 1 2 2 3 3 4 4 5 5   for g_nx=12
!  if the spectrum is being truncated, then we need to 
!  reshuffle the coefficients to look like:  (for g_nx_truncated=8)
!           0  4 1 1 2 2 3 3     where 4 is the cosine mode.
!
!  this is tricky, so for now just set it to zero.  If dealiasing
!  is being used, it will already be set to 0.  If dealiasing is
!  disabled (dealias==0) and im_max < g_nx/2, go in and set to 
!  zero by hand before outputting.  (not yet implemented)
!
!  

subroutine output1_spec(p,pt,buf,fid,fpe,im_max,jm_max,km_max)
!use params
!use mpi

CPOINTER fid
integer :: fpe,im_max,jm_max,km_max       ! cpu to do the I/O
real*8 :: p(nx,ny,nz)
real*8 :: pt(g_nx2,nslabz,ny_2dx)
real*8 :: buf(2+2*im_max,ny_2dx)

! local vars
integer sending_pe,ierr,tag,z_pe,y_pe,x_pe
#ifdef USE_MPI
integer request(2),statuses(MPI_STATUS_SIZE,2)
#endif
integer i,j,k,l,extra_k,dest_pe3(3)
integer n1,n1d,n2,n2d,n3,n3d,jj
integer dealias_nx,dealias_ny,dealias_nz
logical truncation

dealias_nx=min(2+2*im_max,g_nx)
dealias_ny=min(2+2*jm_max,g_ny)
dealias_nz=min(2+2*km_max,g_nz)

truncation=(dealias_nx < g_nx .or. dealias_ny < g_ny .or. dealias_nz < g_nz)
if (dealias==0 .and. truncation) then
   call abortdns("output1_spec: spectral truncation not accurate")
endif

#if 0
do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
!   if (imcord(i)==g_nx/2) print *,imcord(i),jmcord(j),kmcord(k),p(i,j,k)
!   if (jmcord(j)==g_ny/2) print *,imcord(i),jmcord(j),kmcord(k),p(i,j,k)
!   if (kmcord(k)==g_nz/2) print *,imcord(i),jmcord(j),kmcord(k),p(i,j,k)
   if (imcord(i)==g_nx/2) p(i,j,k)=0
   if (jmcord(j)==g_ny/2) p(i,j,k)=0
   if (kmcord(k)==g_nz/2) p(i,j,k)=0
enddo
enddo
enddo
#endif



call transpose_to_x(p,pt,n1,n1d,n2,n2d,n3,n3d)
ASSERT("output1 dimension failure 2",n1==g_nx)
ASSERT("output1 dimension failure 3",n1d==g_nx2)
ASSERT("output1 dimension failure 4",n2==nslabz)
ASSERT("output1 dimension failure 5",n2d==nslabz)
ASSERT("output1 dimension failure 6",n3==ny_2dx)
ASSERT("output1 dimension failure 7",n3d==ny_2dx)


do z_pe=0,ncpu_z-1
do k=1,nslabz
if (z_pe*nslabz+k <= dealias_nz)  then
do y_pe=0,ncpu_y-1
do x_pe=0,ncpu_x-1


   dest_pe3(1)=x_pe
   dest_pe3(2)=y_pe
   dest_pe3(3)=z_pe

#ifdef USE_MPI
   tag=1
   call cart_rank(comm_3d,dest_pe3,sending_pe,ierr)
#else
   sending_pe=my_pe
#endif

   if (sending_pe==my_pe) then

      ! output pt(1:g_nx,k,1:ny_2dx) from cpus: x_pe,y_pe,z_pe
      jj=0
      do j=1,ny_2dx
         if (x_jmcord(j)>jm_max .and. x_jmcord(j)/=(g_ny/2)) exit
         jj=j
      enddo
      l=dealias_nx*jj
      if (l>0) buf(1:dealias_nx,1:jj)=pt(1:dealias_nx,k,1:jj)

      if (my_pe == fpe) then
         ! dont send message to self
      else
#ifdef USE_MPI
         tag=1
         call mpi_isend(l,1,MPI_INTEGER,fpe,tag,comm_3d,request(2),ierr)
         ASSERT("output1: MPI_ISend failure",ierr==0)

         if (l==0) then
            call mpi_waitall(1,request,statuses,ierr) 	
            ASSERT("output1: MPI_waitalll failure",ierr==0)
         else
            call mpi_isend(buf,l,MPI_REAL8,fpe,tag,comm_3d,request(1),ierr)
            ASSERT("output1: MPI_ISend failure",ierr==0)
            call mpi_waitall(2,request,statuses,ierr) 	
            ASSERT("output1: MPI_waitalll failure",ierr==0)
         endif
#endif
      endif
   endif

   if (my_pe==fpe) then
      if (sending_pe==my_pe) then
         ! dont recieve message from self
      else
#ifdef USE_MPI
         call mpi_irecv(l,1,MPI_INTEGER,sending_pe,tag,comm_3d,request,ierr)
         ASSERT("output1: MPI_IRecv failure",ierr==0)
         call mpi_waitall(1,request,statuses,ierr) 	
         ASSERT("output1: MPI_waitalll failure",ierr==0)
         if (l>0) then
            call mpi_irecv(buf,l,MPI_REAL8,sending_pe,tag,comm_3d,request,ierr)
            ASSERT("output1: MPI_IRecv failure",ierr==0)
            call mpi_waitall(1,request,statuses,ierr) 	
            ASSERT("output1: MPI_waitalll failure",ierr==0)
         endif
#endif
      endif
      
      if (l>0) call mwrite8(fid,buf,l)
   endif
enddo
enddo
endif
enddo
enddo

end subroutine






subroutine input1_spec(p,pt,buf,fid,fpe,im_max,jm_max,km_max)
!use params
!use mpi

CPOINTER fid
integer :: fpe,im_max,jm_max,km_max           ! cpu to do the I/O
real*8 :: p(nx,ny,nz)
real*8 :: pt(g_nx2,nslabz,ny_2dx)
real*8 :: buf(2+2*im_max,ny_2dx)

! local vars
integer sending_pe,ierr,tag,z_pe,y_pe,x_pe
#ifdef USE_MPI
integer request(2),statuses(MPI_STATUS_SIZE,2)
#endif
integer i,j,k,l,extra_k,dest_pe3(3)
integer n1,n1d,n2,n2d,n3,n3d,jj
integer dealias_nx,dealias_ny,dealias_nz
logical truncation

dealias_nx=min(2+2*im_max,g_nx)
dealias_ny=min(2+2*jm_max,g_ny)
dealias_nz=min(2+2*km_max,g_nz)

!truncation=(dealias_nx < g_nx .or. dealias_ny < g_ny .or. dealias_nz < g_nz)


pt=0
do z_pe=0,ncpu_z-1
do k=1,nslabz
if (z_pe*nslabz+k <= dealias_nz)  then
do y_pe=0,ncpu_y-1
do x_pe=0,ncpu_x-1


   dest_pe3(1)=x_pe
   dest_pe3(2)=y_pe
   dest_pe3(3)=z_pe

#ifdef USE_MPI
   tag=1
   call cart_rank(comm_3d,dest_pe3,sending_pe,ierr)
#else
   sending_pe=my_pe
#endif

   if (sending_pe==my_pe) then

      ! output pt(1:g_nx,k,1:ny_2dx) from cpus: x_pe,y_pe,z_pe
      jj=0
      do j=1,ny_2dx
         if (x_jmcord(j)>jm_max  .and. x_jmcord(j)/=(g_ny/2)) exit
         jj=j
      enddo
      l=dealias_nx*jj


      if (my_pe == fpe) then
         ! dont send message to self
      else
#ifdef USE_MPI
         tag=1
         call mpi_isend(l,1,MPI_INTEGER,fpe,tag,comm_3d,request,ierr)
         ASSERT("output1: MPI_ISend failure",ierr==0)
         call mpi_waitall(1,request,statuses,ierr) 	
         ASSERT("output1: MPI_waitalll failure",ierr==0)

         if (l>0) then
            call mpi_irecv(buf,l,MPI_REAL8,fpe,tag,comm_3d,request,ierr)
            ASSERT("output1: MPI_ISend failure",ierr==0)
            call mpi_waitall(1,request,statuses,ierr) 	
            ASSERT("output1: MPI_waitalll failure",ierr==0)
            pt(1:dealias_nx,k,1:jj)=buf(1:dealias_nx,1:jj)
         endif
#endif
      endif
   endif

   if (my_pe==fpe) then
      if (sending_pe==my_pe) then
         ! dont recieve message from self
         if (l>0) then
            call mread8(fid,buf,l)
            pt(1:dealias_nx,k,1:jj)=buf(1:dealias_nx,1:jj)
         endif
      else
#ifdef USE_MPI
         call mpi_irecv(l,1,MPI_INTEGER,sending_pe,tag,comm_3d,request,ierr)
         ASSERT("output1: MPI_IRecv failure",ierr==0)
         call mpi_waitall(1,request,statuses,ierr) 	
         ASSERT("output1: MPI_waitalll failure",ierr==0)

         if (l>0) then
            call mread8(fid,buf,l)
            call mpi_isend(buf,l,MPI_REAL8,sending_pe,tag,comm_3d,request,ierr)
            ASSERT("output1: MPI_IRecv failure",ierr==0)
            call mpi_waitall(1,request,statuses,ierr) 	
            ASSERT("output1: MPI_waitalll failure",ierr==0)
         endif
#endif
      endif
   endif


enddo
enddo
endif
enddo
enddo


n1=g_nx
n1d=g_nx2
n2=nslabz
n2d=nslabz
n3=ny_2dx
n3d=ny_2dx
call transpose_from_x(pt,p,n1,n1d,n2,n2d,n3,n3d)


end subroutine







!
!  parallel input from a single file
!  if random==.true., then fill p with random numbers between -1.. 1
!  This is used to get random I.C. that are independent of the
!  parallel decomposition
!
!
subroutine input1(p,pt,buf,fid,fpe_main,random,offset)
!use params
!use mpi
integer :: fpe_main         ! cpu to do all the file I/O

CPOINTER fid
real*8 :: p(nx,ny,nz)
real*8 :: pt(g_nx2,nslabz,ny_2dx)
real*8 :: buf(o_nx,ny_2dx)
integer :: offset,fpe
logical :: first_seek

! local vars
real*8 saved_edge(o_nx)
integer destination_pe,ierr,tag,z_pe,y_pe,x_pe
logical :: random


#ifdef USE_MPI
integer request,statuses(MPI_STATUS_SIZE)
#endif
#ifdef USE_MPI_IO
integer (KIND=MPI_OFFSET_KIND) ::  zpos,ypos
#endif
integer i,j,k,l,extra_k,kuse,dest_pe3(3)
integer n1,n1d,n2,n2d,n3,n3d
integer :: ny_2dx_actual 
ny_2dx_actual = ny_2dx
first_seek=.true.

if (o_nz>g_nz .and. g_bdy_z1/=PERIODIC) then
   call abortdns("output1: cannot handle o_nz=g_nz+1 with non-periodic b.c.")
endif


do z_pe=0,ncpu_z-1
extra_k=0
fpe=fpe_main
if (do_mpi_io) fpe=io_nodes(z_pe)
if (z_pe==ncpu_z-1 .and. o_nz>g_nz) extra_k=1
do k=1,nslabz+extra_k
do y_pe=0,ncpu_y-1
do x_pe=0,ncpu_x-1

   ! for non-perfect load balanced cases, the last few columns
   ! of data will not be real data.
   ! Find the largest value ix<=ny_2dx such that the corresponding
   ! value of i in the 3D decompositoin is i<=nx2
   ! set ny_2dx_actual = ix. 

   ! output pt(1:g_nx,k,1:ny_2dx) from cpus: x_pe,y_pe,z_pe
   l=o_nx*ny_2dx_actual

   if (k>nslabz) then 
      ! this is the periodic z-direction case. refetch data at z=0
      kuse=1
      dest_pe3(1)=x_pe
      dest_pe3(2)=y_pe
      dest_pe3(3)=0
   else
      kuse=k
      dest_pe3(1)=x_pe
      dest_pe3(2)=y_pe
      dest_pe3(3)=z_pe
   endif

#ifdef USE_MPI
   tag=1
   call cart_rank(comm_3d,dest_pe3,destination_pe,ierr)
#else
   destination_pe=my_pe
#endif


   if (my_pe==fpe) then
      if (random) then
         call random_data(buf,o_nx*ny_2dx_actual)
      else
#ifdef USE_MPI_IO
         if (do_mpi_io) then
            zpos = z_pe*nslabz + k-1 
            ypos = y_pe*nslaby + x_pe*ny_2dx_actual
            zpos = input_size*(offset + zpos*o_nx*o_ny+ypos*o_nx)
            if (first_seek) then
               call mpi_file_seek(fid,zpos,MPI_SEEK_SET,ierr)
               first_seek=.false.
            endif
!            call mpi_file_read(fid,buf,o_nx*ny_2dx_actual,MPI_REAL8,statuses,ierr)
         endif
#endif
         call mread8(fid,buf,o_nx*ny_2dx_actual)
      endif

      if (o_ny>g_ny) then
      if (y_pe==ncpu_y-1 .and. x_pe==ncpu_x-1) then    
         ! read and discard periodic duplicate points
         if (random) then
            call random_data(saved_edge,o_nx)
         else
#if 0
            if (do_mpi_io) then
               zpos = z_pe*nslabz + k -1
               ypos = y_pe*nslaby + (1+x_pe)*ny_2dx_actual
               zpos = input_size*(offset + zpos*o_nx*o_ny+ypos*o_nx)
               ! call mpi_file_seek(fid,zpos,MPI_SEEK_SET,ierr)
               ! call mpi_file_read(fid,saved_edge,o_nx,MPI_REAL8,statuses,ierr)
            endif
#endif
            call mread8(fid,saved_edge,o_nx)      
         endif
      endif
      endif

   endif

   if (my_pe==fpe) then
      if (destination_pe==my_pe) then
         ! dont send message to self
      else
#ifdef USE_MPI
         call mpi_isend(buf,l,MPI_REAL8,destination_pe,tag,comm_3d,request,ierr)
         ASSERT("output1: MPI_IRecv failure",ierr==0)
         call mpi_waitall(1,request,statuses,ierr) 	
         ASSERT("output1: MPI_waitalll failure",ierr==0)
#endif
      endif
   endif


   if (destination_pe==my_pe) then
      if (my_pe == fpe) then
         ! dont recieve message from self
      else
#ifdef USE_MPI
         tag=1
         call mpi_irecv(buf,l,MPI_REAL8,fpe,tag,comm_3d,request,ierr)
         ASSERT("output1: MPI_ISend failure",ierr==0)
         call mpi_waitall(1,request,statuses,ierr) 	
         ASSERT("output1: MPI_waitalll failure",ierr==0)
#endif
      endif
      pt(1:g_nx,kuse,:)=buf(1:g_nx,:)
   endif

      
enddo
enddo
enddo
enddo


n1=g_nx
n1d=g_nx2
n2=nslabz
n2d=nslabz
n3=ny_2dx
n3d=ny_2dx
call transpose_from_x(pt,p,n1,n1d,n2,n2d,n3,n3d)


end subroutine








subroutine transpose_from_z_3d(Qhat,q1)
!use params
implicit none
real*8 :: Qhat(g_nz2,nx_2dz,ny_2dz,n_var)
real*8 :: q1(nx,ny,nz,n_var)

! local
integer :: n1,n1d,n2,n2d,n3,n3d,n

n1=g_nz
n1d=g_nz2
n2=nx_2dz
n2d=nx_2dz
n3=ny_2dz
n3d=ny_2dz
do n=1,ndim
   call transpose_from_z(Qhat(1,1,1,n),q1(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo
end subroutine transpose_from_z_3d



subroutine transpose_from_x_3d(Qin,Qout)
!use params
implicit none
real*8 :: Qin(g_nx2,nslabz,ny_2dx,n_var)
real*8 :: Qout(nx,ny,nz,n_var)
! local
integer :: n1,n1d,n2,n2d,n3,n3d,n
n1=g_nx
n1d=g_nx2
n2=nslabz
n2d=nslabz
n3=ny_2dx
n3d=ny_2dx
do n=1,n_var
   call transpose_from_x(Qin(1,1,1,n),Qout(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo
end subroutine transpose_from_x_3d



subroutine transpose_to_x_3d(Qin,Qout)
!use params
implicit none
real*8 :: Qin(nx,ny,nz,n_var)
real*8 :: Qout(g_nx2,nslabz,ny_2dx,n_var)
! local
integer :: n1,n1d,n2,n2d,n3,n3d,n
do n=1,n_var
   call transpose_to_x(Qin(1,1,1,n),Qout(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo
end subroutine transpose_to_x_3d











end module


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  use these routines for all read/writes that may be using MPI_IO,
!  if MPI_IO support was compiled in
!
!  if do_mpi_io==.true:    use MPI calls  (open file with MPI-IO)
!                .false.   use cread/cwrite calls  (open file with copen)
!
!  these routines will also automaticaly convert to single precision if
!  output_size=4
!  if input_size=4:  not yet coded.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mread8(fid,buf,len)
use params
use mpi
implicit none
CPOINTER fid
integer :: len,ierr
real*8 :: buf(len)
character(len=80) :: message
call mread8e(fid,buf,len,ierr)
if (ierr /= len) then
   write(message,*) "ERROR: mread8() elements requested: ",len," elements read: ",ierr
   call abortdns(message)
endif
end subroutine



subroutine mread8e(fid,buf,len,ierr)
use params
use mpi
implicit none
CPOINTER fid
integer :: len,ierr
real*8 :: buf(len)
real*4, allocatable ::  buf4(:)
#ifdef USE_MPI_IO
integer statuses(MPI_STATUS_SIZE)
#endif

if (input_size==4) then
   allocate(buf4(len))
endif

if (do_mpi_io) then 
#ifdef USE_MPI_IO
   if (input_size==4) then	
      call mpi_file_read(fid,buf4,len,MPI_REAL4,statuses,ierr)
   else
      call mpi_file_read(fid,buf,len,MPI_REAL8,statuses,ierr)
   endif
   if (ierr==0) ierr=len  ! return length read, if OK
#else
   call abortdns("MPI_IO support not compiled in")	
#endif
else
   if (input_size==4) then
      call cread4e(fid,buf4,len,ierr)
   else
      call cread8e(fid,buf,len,ierr)
   endif
endif

if (input_size==4) then
   buf(1:len)=buf4(1:len)
   deallocate(buf4)
endif


end subroutine




subroutine mread1e(fid,buf,len,ierr)
use params
use mpi
implicit none
CPOINTER fid
integer :: len,ierr
real*8 :: buf(len)
#ifdef USE_MPI_IO
integer statuses(MPI_STATUS_SIZE)
#endif

if (do_mpi_io) then 
#ifdef USE_MPI_IO
   call mpi_file_read(fid,buf,len,MPI_BYTE,statuses,ierr)
   if (ierr==0) ierr=len  ! return length read, if OK
#else
   call abortdns("MPI_IO support not compiled in")	
#endif
else
   call cread1e(fid,buf,len,ierr)
endif

end subroutine



subroutine mwrite8(fid,buf,len)
use params
use mpi
implicit none
CPOINTER fid
integer :: len
real*8 :: buf(len)
real*4, allocatable ::  buf4(:)
integer :: ierr
#ifdef USE_MPI_IO
integer statuses(MPI_STATUS_SIZE)
#endif

if (output_size==4) then
   allocate(buf4(len))
   buf4(1:len)=buf(1:len)
endif

if (do_mpi_io) then 
#ifdef USE_MPI_IO
   if (output_size==4) then
      call mpi_file_write(fid,buf4,len,MPI_REAL4,statuses,ierr)
   else
      call mpi_file_write(fid,buf,len,MPI_REAL8,statuses,ierr)
   endif
#else
   call abortdns("MPI_IO support not compiled in")	
#endif
else
   if (output_size==4) then
      call cwrite4(fid,buf4,len)
   else
      call cwrite8(fid,buf,len)
   endif
endif

if (output_size==4) deallocate(buf4)
end subroutine




subroutine mwrite1(fid,buf,len)
use params
use mpi
implicit none
CPOINTER fid
integer :: len
real*8 :: buf(len)
integer :: ierr
#ifdef USE_MPI_IO
integer statuses(MPI_STATUS_SIZE)
#endif

if (do_mpi_io) then 
#ifdef USE_MPI_IO
   call mpi_file_write(fid,buf,len,MPI_BYTE,statuses,ierr)
#else
   call abortdns("MPI_IO support not compiled in")	
#endif
else
   call cwrite1(fid,buf,len)
endif

end subroutine
