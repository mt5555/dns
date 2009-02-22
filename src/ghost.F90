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

module ghost
implicit none


real*8,private :: tmx1,tmx2
integer,private :: nghost=2

real*8,private,allocatable :: recbuf0(:)
real*8,private,allocatable :: recbuf1(:)
#ifdef USE_MPI
real*8,private,allocatable :: sendbuf0(:)
real*8,private,allocatable :: sendbuf1(:)
#endif



logical,private,save :: firstcall=.true.

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! check dimensions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ghost_init()
use params
implicit none

integer :: n
!
! dont use bx1,bx2 in these checks, because they will fail
! at a real boundary if offset_bdy is true.  But this is not
! a problem because ghost update does not update ghost cells at
! real boundaries
!

if ( (nx1-1)<nghost ) then 
    call abortdns("nx1 too small for number of ghost cells requested")
endif
if ( (nx-nx2)<nghost) then 
    call abortdns("nx too small for number of ghost cells requested")
endif
if ( (ny1-1)<nghost) then 
    call abortdns("ny1 too small for number of ghost cells requested")
endif
if ( (ny-ny2)<nghost) then 
    call abortdns("ny too small for number of ghost cells requested")
endif
if ( (nz1-1)<nghost .and. ndim==3 ) then 
    call abortdns("nz1 too small for number of ghost cells requested")
endif
if ( (nz-nz2)<nghost .and.  ndim==3 ) then 
    call abortdns("nz too small for number of ghost cells requested")
endif

firstcall=.false.


n = max(ny*nz,nx*ny,nx*nz)*ndim*nghost
allocate (recbuf0(n))
allocate (recbuf1(n))
#ifdef USE_MPI
n=1
if (ncpu_z>1) then
   n=max(n,nx*ny*nghost*ndim)
endif
if (ncpu_x>1) then
   n=max(n,ny*nz*nghost*ndim)
endif
if (ncpu_y>1) then
   n=max(n,nz*nx*nghost*ndim)
endif

allocate (sendbuf0(n))
allocate (sendbuf1(n))
#endif


end subroutine






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! update ghostcells of variable p
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ghost_update_x(p,nvar)
use params
use mpi
implicit none
integer :: nvar
real*8 :: p(nx,ny,nz,nvar)

!
! nx1,nx2 :-> 1,nx
! nslabx -> nx
!

!local variables
integer :: i,j,k,n,l,nmesg,x0,x1,y0,y1,z0,z1
logical :: ghost0=.true.,ghost1=.true.


#ifdef USE_MPI
integer :: ierr,dest_pe0,dest_pe1,request(12),statuses(MPI_STATUS_SIZE,12)
integer :: dest_pe3(3),tag
#endif

if (firstcall) call ghost_init
call wallclock(tmx1)
nmesg=0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! compute X direction ghost cells
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x0=my_x-1
if (x0<0) then
   if (bdy_x1==PERIODIC) then
      x0=ncpu_x-1
   else
      x0=my_x
   endif
endif

x1=my_x+1
if (x1>=ncpu_x) then
   if (bdy_x2==PERIODIC) then
      x1=0
   else
      x1=my_x
   endif
endif


#ifdef USE_MPI
dest_pe3(1)=x0
dest_pe3(2)=my_y
dest_pe3(3)=my_z
call cart_rank(comm_3d,dest_pe3,dest_pe0,ierr)
dest_pe3(1)=x1
call cart_rank(comm_3d,dest_pe3,dest_pe1,ierr)
#endif




! get information of left edge ghost cells of X direction:
if (x0==my_x) then
   l=0
   do n=1,nvar
   do k=bz1,bz2
   do j=by1,by2
   do i=(nghost-1),0,-1
      l=l+1
      if (bdy_x1==PERIODIC) then
         !periodic:  nx2-1,nx2-0  -->   nx1-2,nx1-1
         recbuf0(l)=p(nx2-i,j,k,n)
      else if (bdy_x1==REFLECT) then
         !reflection:  nx1+2,nx1+1  -->   nx1-2,nx1-1
         recbuf0(l)=p(nx1+i+1,j,k,n)
      else if (bdy_x1==REFLECT_ODD) then
         !reflection:  nx1+2,nx1+1  -->   nx1-2,nx1-1
         recbuf0(l)=-p(nx1+i+1,j,k,n)
      else
         ! b.c., so dont touch ghost cells
         ghost0=.false.
         goto 90
      endif
   enddo
   enddo
   enddo
   enddo
else
#ifdef USE_MPI
   l=0
   do n=1,nvar
   do k=bz1,bz2
   do j=by1,by2
   do i=0,nghost-1
      ! regular ghost cell update
      ! nx1,nx1+1 on this processor goes to nx2+1,nx2+2 on dest_pe0
      l=l+1
      sendbuf0(l)=p(nx1+i,j,k,n)
   enddo
   enddo
   enddo
   enddo
   tag=2
   nmesg=nmesg+1
   call mpi_irecv(recbuf0,l,MPI_REAL8,dest_pe0,tag,comm_3d,request(nmesg),ierr)

   tag=1
   nmesg=nmesg+1
   call mpi_isend(sendbuf0,l,MPI_REAL8,dest_pe0,tag,comm_3d,request(nmesg),ierr)
#endif
endif
90 continue



! get information for right edge ghost cells of X direction:
if (x1==my_x) then
   l=0
   do n=1,nvar
   do k=bz1,bz2
   do j=by1,by2
   do i=1,nghost
      l=l+1
      if (bdy_x2==PERIODIC) then
         !periodic:  nx1,nx1+1  -->   nx2+1,nx2+2
         recbuf1(l)=p(nx1+i-1,j,k,n)
      else if (bdy_x2==REFLECT) then
         !reflection:  nx2-1,nx2-2  -->   nx2+1,nx2+2
         recbuf1(l)=p(nx2-i,j,k,n)
      else if (bdy_x2==REFLECT_ODD) then
         !reflection:  nx2-1,nx2-2  -->   nx2+1,nx2+2
         recbuf1(l)=-p(nx2-i,j,k,n)
      else
         ! b.c., so dont touch ghost cells
         ghost1=.false.
         goto 100
      endif
   enddo
   enddo
   enddo
   enddo
else
#ifdef USE_MPI
   l=0
   do n=1,nvar
   do k=bz1,bz2
   do j=by1,by2
   do i=nghost-1,0,-1
      ! regular ghost cell update
      ! nx2-1,nx2 on this processor goes to nx1-2,nx1-1 on dest_pe1
      l=l+1
      sendbuf1(l)=p(nx2-i,j,k,n)
   enddo
   enddo
   enddo
   enddo

   tag=1
   nmesg=nmesg+1
   call mpi_irecv(recbuf1,l,MPI_REAL8,dest_pe1,tag,comm_3d,request(nmesg),ierr)

   tag=2
   nmesg=nmesg+1
   call mpi_isend(sendbuf1,l,MPI_REAL8,dest_pe1,tag,comm_3d,request(nmesg),ierr)
#endif
endif
100 continue




#ifdef USE_MPI
if (nmesg>0) then
   call mpi_waitall(nmesg,request,statuses,ierr) 	
   ASSERT("ghost cell update:  MPI_waitalll failure 1",ierr==0)
endif
#endif


l=0
do n=1,nvar
do k=bz1,bz2
do j=by1,by2
do i=1,nghost
   l=l+1
   ! left edge:  p(nx1-2,,) then p(nx1-1,,)
   if (ghost0) p(nx1-(nghost-i+1),j,k,n)=recbuf0(l)
   ! right edge:  p(nx2+1,,) then p(nx2+2,,)
   if (ghost1) p(nx2+i,j,k,n)=recbuf1(l)
enddo
enddo
enddo
enddo




call wallclock(tmx2)
tims(13)=tims(13)+(tmx2-tmx1)          

end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! update ghostcells of variable p
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ghost_update_y(p,nvar)
use params
use mpi
implicit none
integer :: nvar
real*8 :: p(nx,ny,nz,nvar)


!local variables
integer :: i,j,k,n,l,nmesg,x0,x1,y0,y1,z0,z1
logical :: ghost0=.true.,ghost1=.true.

#ifdef USE_MPI
integer :: ierr,dest_pe0,dest_pe1,request(12),statuses(MPI_STATUS_SIZE,12)
integer :: dest_pe3(3),tag
#endif

if (firstcall) call ghost_init
call wallclock(tmx1)
nmesg=0




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! compute Y direction ghost cells
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
y0=my_y-1
if (y0<0) then
   if (bdy_y1==PERIODIC) then
      y0=ncpu_y-1
   else
      y0=my_y
   endif
endif

y1=my_y+1
if (y1>=ncpu_y) then
   if (bdy_y2==PERIODIC) then
      y1=0
   else
      y1=my_y
   endif
endif


#ifdef USE_MPI
dest_pe3(1)=my_x
dest_pe3(2)=y0
dest_pe3(3)=my_z
call cart_rank(comm_3d,dest_pe3,dest_pe0,ierr)
dest_pe3(2)=y1
call cart_rank(comm_3d,dest_pe3,dest_pe1,ierr)
#endif




! get information of left edge ghost cells of Y direction:
if (y0==my_y) then
   l=0
   do n=1,nvar
   do k=bz1,bz2
   do j=(nghost-1),0,-1
   do i=bx1,bx2
      l=l+1
      if (bdy_y1==PERIODIC) then
         recbuf0(l)=p(i,ny2-j,k,n)
      else if (bdy_y1==REFLECT) then
         recbuf0(l)=p(i,ny1+j+1,k,n)
      else if (bdy_y1==REFLECT_ODD) then
         recbuf0(l)=-p(i,ny1+j+1,k,n)
      else
         ghost0=.false.
         goto 90
      endif
   enddo
   enddo
   enddo
   enddo
else
#ifdef USE_MPI
   l=0
   do n=1,nvar
   do k=bz1,bz2
   do j=0,nghost-1
   do i=bx1,bx2
      ! regular ghost cell update
      l=l+1
      sendbuf0(l)=p(i,ny1+j,k,n)
   enddo
   enddo
   enddo
   enddo
   tag=20
   nmesg=nmesg+1
   call mpi_irecv(recbuf0,l,MPI_REAL8,dest_pe0,tag,comm_3d,request(nmesg),ierr)

   tag=10
   nmesg=nmesg+1
   call mpi_isend(sendbuf0,l,MPI_REAL8,dest_pe0,tag,comm_3d,request(nmesg),ierr)
#endif
endif
90 continue

! get information for right edge ghost cells of X direction:
if (y1==my_y) then
   l=0
   do n=1,nvar
   do k=bz1,bz2
   do j=1,nghost
   do i=bx1,bx2
      l=l+1
      if (bdy_y2==PERIODIC) then
         recbuf1(l)=p(i,ny1+j-1,k,n)
      else if (bdy_y2==REFLECT) then
         recbuf1(l)=p(i,ny2-j,k,n)
      else if (bdy_y2==REFLECT_ODD) then
         recbuf1(l)=-p(i,ny2-j,k,n)
      else
         ghost1=.false.
         goto 100
      endif
   enddo
   enddo
   enddo
   enddo
else
#ifdef USE_MPI
   l=0
   do n=1,nvar
   do k=bz1,bz2
   do j=nghost-1,0,-1
   do i=bx1,bx2
      ! regular ghost cell update
      l=l+1
      sendbuf1(l)=p(i,ny2-j,k,n)
   enddo
   enddo
   enddo
   enddo

   tag=10
   nmesg=nmesg+1
   call mpi_irecv(recbuf1,l,MPI_REAL8,dest_pe1,tag,comm_3d,request(nmesg),ierr)

   tag=20
   nmesg=nmesg+1
   call mpi_isend(sendbuf1,l,MPI_REAL8,dest_pe1,tag,comm_3d,request(nmesg),ierr)
#endif
endif
100 continue





#ifdef USE_MPI
if (nmesg>0) then
   call mpi_waitall(nmesg,request,statuses,ierr) 	
   ASSERT("ghost cell update:  MPI_waitalll failure 1",ierr==0)
endif
#endif


l=0
do n=1,nvar  
do k=bz1,bz2
do j=1,nghost
do i=bx1,bx2
   l=l+1
   ! left edge: 
   if (ghost0) p(i,ny1-(nghost-j+1),k,n)=recbuf0(l)
   ! right edge: 
   if (ghost1) p(i,ny2+j,k,n)=recbuf1(l)
enddo
enddo
enddo
enddo



call wallclock(tmx2)
tims(13)=tims(13)+(tmx2-tmx1)          

end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! update ghostcells of variable p
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ghost_update_z(p,nvar)
use params
use mpi
implicit none
integer :: nvar
real*8 :: p(nx,ny,nz,nvar)


!local variables
integer :: i,j,k,n,l,nmesg,x0,x1,y0,y1,z0,z1
logical :: ghost0=.true.,ghost1=.true.

#ifdef USE_MPI
integer :: ierr,dest_pe0,dest_pe1,request(12),statuses(MPI_STATUS_SIZE,12)
integer :: dest_pe3(3),tag
#endif

if (firstcall) call ghost_init
call wallclock(tmx1)
nmesg=0





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! compute Z direction ghost cells
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (ndim==3) then
z0=my_z-1
if (z0<0) then
   if (bdy_z1==PERIODIC) then
      z0=ncpu_z-1
   else
      z0=my_z
   endif
endif

z1=my_z+1
if (z1>=ncpu_z) then
   if (bdy_z2==PERIODIC) then
      z1=0
   else
      z1=my_z
   endif
endif


#ifdef USE_MPI
dest_pe3(1)=my_x
dest_pe3(2)=my_y
dest_pe3(3)=z0
call cart_rank(comm_3d,dest_pe3,dest_pe0,ierr)
dest_pe3(3)=z1
call cart_rank(comm_3d,dest_pe3,dest_pe1,ierr)
#endif




! get information of left edge ghost cells of Z direction:
if (z0==my_z) then
   l=0
   do n=1,nvar
   do k=(nghost-1),0,-1
   do j=by1,by2
   do i=bx1,bx2
      l=l+1
      if (bdy_z1==PERIODIC) then
         recbuf0(l)=p(i,j,nz2-k,n)
      else if (bdy_z1==REFLECT) then
         recbuf0(l)=p(i,j,nz1+k+1,n)
      else if (bdy_z1==REFLECT_ODD) then
         recbuf0(l)=-p(i,j,nz1+k+1,n)
      else
         ghost0=.false.
         goto 90
      endif
   enddo
   enddo
   enddo
   enddo
else
#ifdef USE_MPI
   l=0
   do n=1,nvar
   do k=0,nghost-1
   do j=by1,by2
   do i=bx1,bx2
      ! regular ghost cell update
      l=l+1
      sendbuf0(l)=p(i,j,nz1+k,n)
   enddo
   enddo
   enddo
   enddo
   tag=200
   nmesg=nmesg+1
   call mpi_irecv(recbuf0,l,MPI_REAL8,dest_pe0,tag,comm_3d,request(nmesg),ierr)

   tag=100
   nmesg=nmesg+1
   call mpi_isend(sendbuf0,l,MPI_REAL8,dest_pe0,tag,comm_3d,request(nmesg),ierr)
#endif
endif
90 continue

! get information for right edge ghost cells of X direction:
if (z1==my_z) then
   l=0
   do n=1,nvar
   do k=1,nghost
   do j=by1,by2
   do i=bx1,bx2
      l=l+1
      if (bdy_z2==PERIODIC) then
         recbuf1(l)=p(i,j,nz1+k-1,n)
      else if (bdy_z2==REFLECT) then
         recbuf1(l)=p(i,j,nz2-k,n)
      else if (bdy_z2==REFLECT_ODD) then
         recbuf1(l)=-p(i,j,nz2-k,n)
      else
         ghost1=.false.
         goto 100
      endif
   enddo
   enddo
   enddo
   enddo
else
#ifdef USE_MPI
   l=0
   do n=1,nvar
   do k=nghost-1,0,-1
   do j=by1,by2
   do i=bx1,bx2
      ! regular ghost cell update
      l=l+1
      sendbuf1(l)=p(i,j,nz2-k,n)
   enddo
   enddo
   enddo
   enddo

   tag=100
   nmesg=nmesg+1
   call mpi_irecv(recbuf1,l,MPI_REAL8,dest_pe1,tag,comm_3d,request(nmesg),ierr)

   tag=200
   nmesg=nmesg+1
   call mpi_isend(sendbuf1,l,MPI_REAL8,dest_pe1,tag,comm_3d,request(nmesg),ierr)
#endif
endif
endif
100 continue




#ifdef USE_MPI
if (nmesg>0) then
   call mpi_waitall(nmesg,request,statuses,ierr) 	
   ASSERT("ghost cell update:  MPI_waitalll failure 1",ierr==0)
endif
#endif

if (ndim==3) then
l=0
do n=1,nvar
do k=1,nghost
do j=by1,by2
do i=bx1,bx2
   l=l+1
   ! left edge: 
   if (ghost0) p(i,j,nz1-(nghost-k+1),n)=recbuf0(l)
   ! right edge: 
   if (ghost1) p(i,j,nz2+k,n)=recbuf1(l)
enddo
enddo
enddo
enddo
endif



call wallclock(tmx2)
tims(13)=tims(13)+(tmx2-tmx1)          

end subroutine

end module








