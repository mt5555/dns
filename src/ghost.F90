#include "macros.h"

module ghost
use params
use mpi
implicit none


real*8,private :: tmx1,tmx2

integer,private :: nghost=2


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! update ghostcells of variable p
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ghost_update(p,nvar)
use params
implicit none
integer :: nvar
real*8 :: p(nx,ny,nz,nvar)


!local variables
integer :: i,j,k,n,l,nmesg,x1,x2,y1,y2,z1,z2
real*8 :: recbufx1(nslaby*nslabz*nvar)
real*8 :: recbufx2(nslaby*nslabz*nvar)

#ifdef USE_MPI
real*8 :: sendbufx1(nslaby*nslabz*nvar)
real*8 :: sendbufx2(nslaby*nslabz*nvar)
integer :: ndata
integer :: ierr,dest_pe1,dest_pe2,request(12),statuses(MPI_STATUS_SIZE,12)
integer :: dest_pe3(3),tag
#endif

call wallclock(tmx1)
nmesg=0

!
! compute X direction ghost cells
!
x1=my_x-1
if (x1<0) then
   ! periodic
   x1=ncpu_x
   ! reflection
   ! x1=my_x
endif

x2=my_x+1
if (x2>0) then
   ! periodic
   x2=0
   ! reflection
   ! x2=my_x
endif


#ifdef USE_MPI
ndata=nvar*nslaby*nslabz*nghost
tag=1
dest_pe3(1)=x1
dest_pe3(2)=my_y
dest_pe3(3)=my_z
call mpi_cart_rank(comm_3d,dest_pe3,dest_pe1,ierr)
dest_pe3(1)=x2
call mpi_cart_rank(comm_3d,dest_pe3,dest_pe2,ierr)
#endif



! get information of left edge ghost cells of X direction:
if (x1==my_x) then
   l=0
   do n=1,nvar
   do k=nz1,nz2
   do j=ny1,ny2
   do i=(nghost-1),0,-1
      l=l+1
      !periodic:  nx2-1,nx2-0  -->   nx1-2,nx1-1
      recbufx1(l)=p(nx2-i,j,k,n)
      !reflection:  nx1+2,nx1+1  -->   nx1-2,nx1-1
      !recbufx1(l)=p(nx1+i+1,j,k,n)
   enddo
   enddo
   enddo
   enddo
else
#ifdef USE_MPI
   l=0
   do n=1,nvar
   do k=nz1,nz2
   do j=ny1,ny2
   do i=0,nghost-1
      ! regular ghost cell update
      ! nx1,nx1+1 on this processor goes to nx2+1,nx2+2 on dest_pe1
      l=l+1
      sendbufx1(l)=p(nx1+i,j,k,n)
   enddo
   enddo
   enddo
   enddo
   nmesg=nmesg+1
   call MPI_IRecv(recbufx1,ndata,MPI_REAL8,dest_pe1,tag,comm_3d,request(nmesg),ierr)
   nmesg=nmesg+1
   call MPI_ISend(sendbufx1,ndata,MPI_REAL8,dest_pe1,tag,comm_3d,request(nmesg),ierr)
#endif
endif


! get information for right edge ghost cells of X direction:
if (x2==my_x) then
   l=0
   do n=1,nvar
   do k=nz1,nz2
   do j=ny1,ny2
   do i=1,nghost
      l=l+1
      !periodic:  nx1,nx1+1  -->   nx2+1,nx2+2
      recbufx2(l)=p(nx1+i-1,j,k,n)
      !reflection:  nx2-1,nx2-2  -->   nx2+1,nx2+2
      !recbufx2(l)=p(nx2-i,j,k,n)
   enddo
   enddo
   enddo
   enddo
else
#ifdef USE_MPI
   l=0
   do n=1,nvar
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nghost-1,0,-1
      ! regular ghost cell update
      ! nx2-1,nx2 on this processor goes to nx1-2,nx1-1 on dest_pe2
      l=l+1
      sendbufx2(l)=p(nx2-i,j,k,n)
   enddo
   enddo
   enddo
   enddo
   nmesg=nmesg+1
   call MPI_IRecv(recbufx2,ndata,MPI_REAL8,dest_pe2,tag,comm_3d,request(nmesg),ierr)
   nmesg=nmesg+1
   call MPI_ISend(sendbufx2,ndata,MPI_REAL8,dest_pe2,tag,comm_3d,request(nmesg),ierr)
#endif
endif











#ifdef USE_MPI
call MPI_waitall(nmesg,request,statuses,ierr) 	
ASSERT("ghost cell update:  MPI_waitalll failure 1",ierr==0)
#endif

l=0
do n=1,nvar
do k=nz1,nz2
do j=ny1,ny2
do i=1,nghost
   l=l+1
   ! left edge:  p(nx1-2,,) then p(nx1-1,,)
   p(nx1-(nghost-i+1),j,k,n)=recbufx1(l)
   ! right edge:  p(nx2+1,,) then p(nx2+2,,)
   p(nx2+i,j,k,n)=recbufx2(l)
enddo
enddo
enddo
enddo



end subroutine

end module
