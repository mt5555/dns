#include "macros.h"

module ghost
implicit none


real*8,private :: tmx1,tmx2
integer,private :: nghost=2

logical,save :: firstcall=.true.

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! check dimensions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ghost_init()
use params
implicit none

!
! dont use bx1,bx2 in these checks, because they will fail
! at a real boundary if offset_bdy is true.  But this is not
! a problem because ghost update does not update ghost cells at
! real boundaries
!

if ( (nx1-1)<nghost ) then 
    call abort("nx1 too small for number of ghost cells requested")
endif
if ( (nx-nx2)<nghost) then 
    call abort("nx too small for number of ghost cells requested")
endif
if ( (ny1-1)<nghost) then 
    call abort("ny1 too small for number of ghost cells requested")
endif
if ( (ny-ny2)<nghost) then 
    call abort("ny too small for number of ghost cells requested")
endif
if ( (nz1-1)<nghost .and. ndim==3 ) then 
    call abort("nz1 too small for number of ghost cells requested")
endif
if ( (nz-nz2)<nghost .and.  ndim==3 ) then 
    call abort("nz too small for number of ghost cells requested")
endif

firstcall=.false.

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
real*8 :: recbufx0(ny*nz*nvar*nghost)
real*8 :: recbufx1(ny*nz*nvar*nghost)
logical :: ghost0=.true.,ghost1=.true.

#ifdef USE_MPI
real*8 :: sendbufx0(ny*nz*nvar*nghost)
real*8 :: sendbufx1(ny*nz*nvar*nghost)
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
call mpi_cart_rank(comm_3d,dest_pe3,dest_pe0,ierr)
dest_pe3(1)=x1
call mpi_cart_rank(comm_3d,dest_pe3,dest_pe1,ierr)
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
         recbufx0(l)=p(nx2-i,j,k,n)
      else if (bdy_x1==REFLECT) then
         !reflection:  nx1+2,nx1+1  -->   nx1-2,nx1-1
         recbufx0(l)=p(nx1+i+1,j,k,n)
      else if (bdy_x1==REFLECT_ODD) then
         !reflection:  nx1+2,nx1+1  -->   nx1-2,nx1-1
         recbufx0(l)=-p(nx1+i+1,j,k,n)
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
      sendbufx0(l)=p(nx1+i,j,k,n)
   enddo
   enddo
   enddo
   enddo
   tag=2
   nmesg=nmesg+1
   call MPI_IRecv(recbufx0,l,MPI_REAL8,dest_pe0,tag,comm_3d,request(nmesg),ierr)

   tag=1
   nmesg=nmesg+1
   call MPI_ISend(sendbufx0,l,MPI_REAL8,dest_pe0,tag,comm_3d,request(nmesg),ierr)
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
         recbufx1(l)=p(nx1+i-1,j,k,n)
      else if (bdy_x2==REFLECT) then
         !reflection:  nx2-1,nx2-2  -->   nx2+1,nx2+2
         recbufx1(l)=p(nx2-i,j,k,n)
      else if (bdy_x2==REFLECT_ODD) then
         !reflection:  nx2-1,nx2-2  -->   nx2+1,nx2+2
         recbufx1(l)=-p(nx2-i,j,k,n)
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
      sendbufx1(l)=p(nx2-i,j,k,n)
   enddo
   enddo
   enddo
   enddo

   tag=1
   nmesg=nmesg+1
   call MPI_IRecv(recbufx1,l,MPI_REAL8,dest_pe1,tag,comm_3d,request(nmesg),ierr)

   tag=2
   nmesg=nmesg+1
   call MPI_ISend(sendbufx1,l,MPI_REAL8,dest_pe1,tag,comm_3d,request(nmesg),ierr)
#endif
endif
100 continue




#ifdef USE_MPI
call MPI_waitall(nmesg,request,statuses,ierr) 	
ASSERT("ghost cell update:  MPI_waitalll failure 1",ierr==0)
#endif


l=0
do n=1,nvar
do k=bz1,bz2
do j=by1,by2
do i=1,nghost
   l=l+1
   ! left edge:  p(nx1-2,,) then p(nx1-1,,)
   if (ghost0) p(nx1-(nghost-i+1),j,k,n)=recbufx0(l)
   ! right edge:  p(nx2+1,,) then p(nx2+2,,)
   if (ghost1) p(nx2+i,j,k,n)=recbufx1(l)
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
real*8 :: recbufy0(nx*nz*nvar*nghost)
real*8 :: recbufy1(nx*nz*nvar*nghost)
logical :: ghost0=.true.,ghost1=.true.

#ifdef USE_MPI
real*8 :: sendbufy0(nx*nz*nvar*nghost)
real*8 :: sendbufy1(nx*nz*nvar*nghost)
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
call mpi_cart_rank(comm_3d,dest_pe3,dest_pe0,ierr)
dest_pe3(2)=y1
call mpi_cart_rank(comm_3d,dest_pe3,dest_pe1,ierr)
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
         recbufy0(l)=p(i,ny2-j,k,n)
      else if (bdy_y1==REFLECT) then
         recbufy0(l)=p(i,ny1+j+1,k,n)
      else if (bdy_y1==REFLECT_ODD) then
         recbufy0(l)=-p(i,ny1+j+1,k,n)
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
      sendbufy0(l)=p(i,ny1+j,k,n)
   enddo
   enddo
   enddo
   enddo
   tag=20
   nmesg=nmesg+1
   call MPI_IRecv(recbufy0,l,MPI_REAL8,dest_pe0,tag,comm_3d,request(nmesg),ierr)

   tag=10
   nmesg=nmesg+1
   call MPI_ISend(sendbufy0,l,MPI_REAL8,dest_pe0,tag,comm_3d,request(nmesg),ierr)
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
         recbufy1(l)=p(i,ny1+j-1,k,n)
      else if (bdy_y2==REFLECT) then
         recbufy1(l)=p(i,ny2-j,k,n)
      else if (bdy_y2==REFLECT_ODD) then
         recbufy1(l)=-p(i,ny2-j,k,n)
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
      sendbufy1(l)=p(i,ny2-j,k,n)
   enddo
   enddo
   enddo
   enddo

   tag=10
   nmesg=nmesg+1
   call MPI_IRecv(recbufy1,l,MPI_REAL8,dest_pe1,tag,comm_3d,request(nmesg),ierr)

   tag=20
   nmesg=nmesg+1
   call MPI_ISend(sendbufy1,l,MPI_REAL8,dest_pe1,tag,comm_3d,request(nmesg),ierr)
#endif
endif
100 continue





#ifdef USE_MPI
call MPI_waitall(nmesg,request,statuses,ierr) 	
ASSERT("ghost cell update:  MPI_waitalll failure 1",ierr==0)
#endif


l=0
do n=1,nvar  
do k=bz1,bz2
do j=1,nghost
do i=bx1,bx2
   l=l+1
   ! left edge: 
   if (ghost0) p(i,ny1-(nghost-j+1),k,n)=recbufy0(l)
   ! right edge: 
   if (ghost1) p(i,ny2+j,k,n)=recbufy1(l)
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
real*8 :: recbufz0(nx*ny*nvar*nghost)
real*8 :: recbufz1(nx*ny*nvar*nghost)
logical :: ghost0=.true.,ghost1=.true.

#ifdef USE_MPI
real*8 :: sendbufz0(nx*ny*nvar*nghost)
real*8 :: sendbufz1(nx*ny*nvar*nghost)
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
call mpi_cart_rank(comm_3d,dest_pe3,dest_pe0,ierr)
dest_pe3(3)=z1
call mpi_cart_rank(comm_3d,dest_pe3,dest_pe1,ierr)
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
         recbufz0(l)=p(i,j,nz2-k,n)
      else if (bdy_z1==REFLECT) then
         recbufz0(l)=p(i,j,nz1+k+1,n)
      else if (bdy_z1==REFLECT_ODD) then
         recbufz0(l)=-p(i,j,nz1+k+1,n)
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
      sendbufz0(l)=p(i,j,nz1+k,n)
   enddo
   enddo
   enddo
   enddo
   tag=200
   nmesg=nmesg+1
   call MPI_IRecv(recbufz0,l,MPI_REAL8,dest_pe0,tag,comm_3d,request(nmesg),ierr)

   tag=100
   nmesg=nmesg+1
   call MPI_ISend(sendbufz0,l,MPI_REAL8,dest_pe0,tag,comm_3d,request(nmesg),ierr)
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
         recbufz1(l)=p(i,j,nz1+k-1,n)
      else if (bdy_z2==REFLECT) then
         recbufz1(l)=p(i,j,nz2-k,n)
      else if (bdy_z2==REFLECT_ODD) then
         recbufz1(l)=-p(i,j,nz2-k,n)
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
      sendbufz1(l)=p(i,j,nz2-k,n)
   enddo
   enddo
   enddo
   enddo

   tag=100
   nmesg=nmesg+1
   call MPI_IRecv(recbufz1,l,MPI_REAL8,dest_pe1,tag,comm_3d,request(nmesg),ierr)

   tag=200
   nmesg=nmesg+1
   call MPI_ISend(sendbufz1,l,MPI_REAL8,dest_pe1,tag,comm_3d,request(nmesg),ierr)
#endif
endif
endif
100 continue




#ifdef USE_MPI
call MPI_waitall(nmesg,request,statuses,ierr) 	
ASSERT("ghost cell update:  MPI_waitalll failure 1",ierr==0)
#endif

if (ndim==3) then
l=0
do n=1,nvar
do k=1,nghost
do j=by1,by2
do i=bx1,bx2
   l=l+1
   ! left edge: 
   if (ghost0) p(i,j,nz1-(nghost-k+1),n)=recbufz0(l)
   ! right edge: 
   if (ghost1) p(i,j,nz2+k,n)=recbufz1(l)
enddo
enddo
enddo
enddo
endif



call wallclock(tmx2)
tims(13)=tims(13)+(tmx2-tmx1)          

end subroutine

end module








