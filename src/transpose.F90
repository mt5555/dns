#include "macros.h"

module transpose
use params
use mpi
implicit none


real*8,private :: tmx1,tmx2
real*8,private,allocatable :: sendbuf(:),recbuf(:)

contains

subroutine transpose_init
#ifdef USE_MPI
integer max_size
max_size=1
if (ncpu_z>1) then
   max_size=max(max_size,nslabx*nslabz*ny_2dz)
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



#undef PBLOCK
#ifdef PBLOCK

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
use params
real*8 p(nx,ny,nz)
real*8 pt(g_nz2,nslabx,ny_2dz)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc,iproc2
integer i,j,k,jj,l
integer,parameter :: maxbuf=min(PBLOCK,ncpu_z)   ! maximum number of MPI buffers
integer :: n          ! loop over message blocks
integer :: ncount     ! number of message blocks
integer :: nb         ! loop over MPI buffers
integer :: nbcount    ! actual number of MPI buffers

#ifdef USE_MPI
real*8 :: sendbuf(nslabx*nslabz*ny_2dz,maxbuf)
real*8 :: recbuf(nslabx*nslabz*ny_2dz,maxbuf)
integer :: ierr,dest_pe,request(maxbuf,2),statuses(MPI_STATUS_SIZE,maxbuf,2)
integer :: dest_pe3(3),tag
#endif

call wallclock(tmx1)

!
! each cube is broken, along the y axis, into mpidims(3) slabs of
! size ny_2dz = (ny2-ny1+1)/mpidims(3)
!
! in the z direction, the dimension is nslabz = (nz2-nz1+1)
!

n1=g_nz
n1d=g_nz2   	
n2=nslabx
n2d=nslabx
n3=ny_2dz
n3d=ny_2dz


! do the message to one's self
iproc=my_z
do j=1,ny_2dz  ! loop over points in a single slab
   jj=ny1 + iproc*ny_2dz +j -1
   !ASSERT("transpose_to_z jj failure 1s",jj<=ny2)
   ASSERT("transpose_to_z jj failure 2s",jj>=ny1)
   do k=nz1,nz2
   do i=nx1,nx2
      if (jj>ny2) then
         pt(k+iproc*nslabz-nz1+1,i-nx1+1,j)=0
      else
         pt(k+iproc*nslabz-nz1+1,i-nx1+1,j)=p(i,jj,k)
      endif
   enddo
   enddo
enddo


#ifdef USE_MPI
iproc2=0
do  ! loop over blocks of messages of size nbcount

   ncount=min(maxbuf, mpidims(3)-iproc2)
   nb=0
   do n=1,ncount 
      iproc=iproc2+n-1
      if (iproc/=my_z) then
         nb=nb+1
         dest_pe3(1)=my_x
         dest_pe3(2)=my_y
         dest_pe3(3)=iproc
         call mpi_cart_rank(comm_3d,dest_pe3,dest_pe,ierr)
         ASSERT("transpose_to_z: MPI_cart_rank failure",ierr==0)

         tag=my_z
         L=nslabz*nslabx*ny_2dz
         call MPI_IRecv(recbuf(1,nb),L,MPI_REAL8,dest_pe,tag,comm_3d,request(nb,1),ierr)
         ASSERT("transpose_to_z: MPI_IRecv failure 1",ierr==0)

         L=0
         do j=1,ny_2dz  ! loop over points in a single slab
            jj=ny1 + iproc*ny_2dz +j -1
            !ASSERT("transpose_to_z jj failure 1",jj<=ny2)
            ASSERT("transpose_to_z jj failure 2",jj>=ny1)
            do k=nz1,nz2
               do i=nx1,nx2
                  L=L+1
                  if (jj>ny2) then
                     sendbuf(L,nb)=0
                  else
                     sendbuf(L,nb)=p(i,jj,k)
                  endif
               enddo
            enddo
         enddo
         
         tag=iproc
         call MPI_ISend(sendbuf(1,nb),L,MPI_REAL8,dest_pe,tag,comm_3d,request(nb,2),ierr)
         ASSERT("transpose_to_z: MPI_ISend failure 1",ierr==0)
      endif
   enddo


   ! wait on all the receives
   nbcount=nb
   if (nbcount>0) then
      ! wait for all receives to finish
      call MPI_waitall(nbcount,request(1,1),statuses,ierr) 	
      ASSERT("transpose_to_z: MPI_waitalll failure 1",ierr==0)

      ! copy data out of receive buffers   
      nb=0
      do n=1,ncount
         iproc=iproc2+n-1
         if (iproc/=my_z) then
            nb=nb+1
            L=0
            do j=1,ny_2dz  ! loop over points in a single slab
               do k=nz1,nz2
                  do i=nx1,nx2
                     L=L+1
                     pt(k+iproc*nslabz-nz1+1,i-nx1+1,j)=recbuf(L,nb)
                  enddo
               enddo
            enddo
         endif
      enddo

      ! wait for all sends to finish
      call MPI_waitall(nbcount,request(1,2),statuses,ierr) 	
      ASSERT("transpose_to_z: MPI_waitalll failure 1",ierr==0)
   endif


   iproc2=iproc2+ncount
   ASSERT("iproc2 error 1 ",iproc2<=mpidims(3))
   if (iproc2==mpidims(3)) exit
enddo

#endif


call wallclock(tmx2) 
tims(6)=tims(6)+(tmx2-tmx1)          
end subroutine

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
use params
real*8 p(nx,ny,nz)
real*8 pt(g_nz2,nslabx,ny_2dz)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc
integer i,j,k,jj,l
#ifdef USE_MPI
!real*8 sendbuf(nslabx*nslabz*ny_2dz)
!real*8 recbuf(nslabx*nslabz*ny_2dz)
integer ierr,dest_pe,request(2),statuses(MPI_STATUS_SIZE,2)
integer dest_pe3(3),tag
#endif

call wallclock(tmx1)

!
! each cube is broken, along the y axis, into mpidims(3) slabs of
! size ny_2dz = (ny2-ny1+1)/mpidims(3)
!
! in the z direction, the dimension is nslabz = (nz2-nz1+1)
!

n1=g_nz
n1d=g_nz2   	
n2=nslabx
n2d=nslabx
n3=ny_2dz
n3d=ny_2dz



do iproc=0,mpidims(3)-1  ! loop over each slab

   if (iproc==my_z) then
      do j=1,ny_2dz  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2dz +j -1
         !ASSERT("transpose_to_z jj failure 1",jj<=ny2)
         ASSERT("transpose_to_z jj failure 2",jj>=ny1)
         do k=nz1,nz2
         do i=nx1,nx2
            if (jj>ny2) then
               pt(k+iproc*nslabz-nz1+1,i-nx1+1,j)=0
            else
               pt(k+iproc*nslabz-nz1+1,i-nx1+1,j)=p(i,jj,k)
            endif
         enddo
         enddo
      enddo
   else
#ifdef USE_MPI
      dest_pe3(1)=my_x
      dest_pe3(2)=my_y
      dest_pe3(3)=iproc
      call mpi_cart_rank(comm_3d,dest_pe3,dest_pe,ierr)
      ASSERT("transpose_to_z: MPI_cart_rank failure",ierr==0)

      l=0
      do j=1,ny_2dz  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2dz +j -1
         !ASSERT("transpose_to_z jj failure 1",jj<=ny2)
         ASSERT("transpose_to_z jj failure 2",jj>=ny1)
         do k=nz1,nz2
         do i=nx1,nx2
            l=l+1
            if (jj>ny2) then
               sendbuf(l)=0
            else
               sendbuf(l)=p(i,jj,k)
            endif
         enddo
         enddo
      enddo

!     send/rec buffer to (my_x,my_y,iproc)
      tag=my_z
      call MPI_IRecv(recbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(1),ierr)
      ASSERT("transpose_to_z: MPI_IRecv failure 1",ierr==0)
      tag=iproc
      call MPI_ISend(sendbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(2),ierr)
      ASSERT("transpose_to_z: MPI_ISend failure 1",ierr==0)
      call MPI_waitall(2,request,statuses,ierr) 	
      ASSERT("transpose_to_z: MPI_waitalll failure 1",ierr==0)

      l=0
      do j=1,ny_2dz  ! loop over points in a single slab
         do k=nz1,nz2
         do i=nx1,nx2
	    l=l+1
            pt(k+iproc*nslabz-nz1+1,i-nx1+1,j)=recbuf(l)
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
! transform data to a cartesian decomposition from a 2D decomposition
! with z as the leading index.
!
! input: pt, and the dimensions of pt: n1,n1d,n2,n2d,n3,n3d
! ouput: p 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_from_z(pt,p,n1,n1d,n2,n2d,n3,n3d)
use params
real*8 p(nx,ny,nz)
real*8 pt(g_nz2,nslabx,ny_2dz)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc
integer i,j,k,jj,l
#ifdef USE_MPI
!real*8 sendbuf(nslabx*nslabz*ny_2dz)
!real*8 recbuf(nslabx*nslabz*ny_2dz)
integer ierr,dest_pe,request(2),statuses(MPI_STATUS_SIZE,2)
integer dest_pe3(3),tag
#endif

!
! each cube is broken, along the y axis, into mpidims(3) slabs of
! size ny_2dz = (ny2-ny1+1)/mpidims(3)
!
! in the z direction, the dimension is nslabz = (nz2-nz1+1)

call wallclock(tmx1)

! If any of these fail, then pt was probably not computed
! via a call to transpose_to_z().
ASSERT("transpose_from_z dimension failure 2",n1==g_nz)
ASSERT("transpose_from_z dimension failure 3",n1d==g_nz2)
ASSERT("transpose_from_z dimension failure 4",n2==nslabx)
ASSERT("transpose_from_z dimension failure 5",n2d==nslabx)
ASSERT("transpose_from_z dimension failure 6",n3==ny_2dz)
ASSERT("transpose_from_z dimension failure 7",n3d==ny_2dz)

do iproc=0,mpidims(3)-1  ! loop over each slab


   if (iproc==my_z) then
      do j=1,ny_2dz  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2dz +j -1
         !ASSERT("transpose_from_z jj failure 1",jj<=ny2)
         ASSERT("transpose_from_z jj failure 2",jj>=ny1)
         if (jj<=ny2) then
         do k=nz1,nz2
         do i=nx1,nx2
            p(i,jj,k)=pt(k+iproc*nslabz-nz1+1,i-nx1+1,j)
         enddo
         enddo
         endif
      enddo
   else

#ifdef USE_MPI
      dest_pe3(1)=my_x
      dest_pe3(2)=my_y
      dest_pe3(3)=iproc
      call mpi_cart_rank(comm_3d,dest_pe3,dest_pe,ierr)
      ASSERT("transpose_from_z: MPI_cart_rank failure 1",ierr==0)

      l=0
      do j=1,ny_2dz  ! loop over points in a single slab
         do k=nz1,nz2
         do i=nx1,nx2
            l=l+1
            sendbuf(l)=pt(k+iproc*nslabz-nz1+1,i-nx1+1,j)
         enddo
         enddo
      enddo

!     send/rec


      tag=my_z
      call MPI_IRecv(recbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(1),ierr)
      ASSERT("transpose_from_z: MPI_IRecv failure 1",ierr==0)
      tag=iproc
      call MPI_ISend(sendbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(2),ierr)
      ASSERT("transpose_from_z: MPI_ISend failure 1",ierr==0)
      call MPI_waitall(2,request,statuses,ierr) 	
      ASSERT("transpose_from_z: MPI_waitalll failure 1",ierr==0)


      l=0
      do j=1,ny_2dz  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2dz +j -1
         !ASSERT("transpose_from_z jj failure 1",jj<=ny2)
         ASSERT("transpose_from_z jj failure 2",jj>=ny1)
         do k=nz1,nz2
         do i=nx1,nx2
	    l=l+1
            if (jj<=ny2) p(i,jj,k)=recbuf(l)
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
! transform data from a cartesian decomposition to a 2D decomposition
! with x as the leading index.  
! slabs are chopped up in the y direction.
!
! input: pt, and the dimensions of pt: n1,n1d,n2,n2d,n3,n3d
! ouput: p 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_to_x(p,pt,n1,n1d,n2,n2d,n3,n3d)
use params
real*8 p(nx,ny,nz)
real*8 pt(g_nx2,nslabz,ny_2dx)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc
integer i,j,k,jj,l
#ifdef USE_MPI
!real*8 sendbuf(nslabx*nslabz*ny_2dx)
!real*8 recbuf(nslabx*nslabz*ny_2dx)
integer ierr,dest_pe,request(2),statuses(MPI_STATUS_SIZE,2)
integer dest_pe3(3),tag
#endif

call wallclock(tmx1)

!
! each cube is broken, along the z axis, into mpidims(1) slabs of
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


do iproc=0,mpidims(1)-1  ! loop over each slab
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
      call mpi_cart_rank(comm_3d,dest_pe3,dest_pe,ierr)
      ASSERT("transpose_to_x2: MPI_cart_rank failure 1",ierr==0)


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
      tag=my_x
      call MPI_IRecv(recbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(1),ierr)
      ASSERT("transpose_to_x2: MPI_IRecv failure 1",ierr==0)
      tag=iproc
      call MPI_ISend(sendbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(2),ierr)
      ASSERT("transpose_to_x2: MPI_ISend failure 1",ierr==0)
      call MPI_waitall(2,request,statuses,ierr) 	
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
! transform data to a cartesian decomposition from a 2D decomposition
! with x as the leading index.  same as transpose_from_x, but
! slabs are chopped up in the y direction istead of z direction.
!
! input: pt, and the dimensions of pt: n1,n1d,n2,n2d,n3,n3d
! ouput: p 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_from_x(pt,p,n1,n1d,n2,n2d,n3,n3d)
use params

real*8 p(nx,ny,nz)
real*8 pt(g_nx2,nslabz,ny_2dx)
integer n1,n1d,n2,n2d,n3,n3d


!local variables
integer iproc
integer i,j,k,jj,l
#ifdef USE_MPI
!real*8 sendbuf(nslabx*nslabz*ny_2dx)
!real*8 recbuf(nslabx*nslabz*ny_2dx)
integer ierr,dest_pe,request(2),statuses(MPI_STATUS_SIZE,2)
integer dest_pe3(3),tag
#endif

call wallclock(tmx1)

!
! each cube is broken, along the z axis, into mpidims(1) slabs of
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


do iproc=0,mpidims(1)-1  ! loop over each slab
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
      call mpi_cart_rank(comm_3d,dest_pe3,dest_pe,ierr)
      ASSERT("transpose_from_x: MPI_cart_rank failure 1",ierr==0)


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
      tag=my_x
      call MPI_IRecv(recbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(1),ierr)
      ASSERT("transpose_from_x2: MPI_IRecv failure 1",ierr==0)
      tag=iproc
      call MPI_ISend(sendbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(2),ierr)
      ASSERT("transpose_from_x2: MPI_ISend failure 1",ierr==0)
      call MPI_waitall(2,request,statuses,ierr) 	
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
! transform data in a cartesian decomposition into a 2D decomposition
! with y as the leading index.
!
! input: p
! ouput: pt, and the dimensions of pt: n1,n1d,n2,n2d,n3,n3d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_to_y(p,pt,n1,n1d,n2,n2d,n3,n3d)
use params

real*8 p(nx,ny,nz)
real*8 pt(g_ny2,nslabz,nx_2dy)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc
integer i,j,k,ii,l
#ifdef USE_MPI
!real*8 sendbuf(nslaby*nslabz*nx_2dy)
!real*8 recbuf(nslaby*nslabz*nx_2dy)
integer ierr,dest_pe,request(2),statuses(MPI_STATUS_SIZE,2)
integer dest_pe3(3),tag
#endif

call wallclock(tmx1)


!
! each cube is broken, along the x axis, into mpidims(2) slabs of
! size nx_2dy
!
! in the y direction, the dimension is nslaby = (ny2-ny1+1)


n1=g_ny
n1d=g_ny2
n2=nslabz
n2d=nslabz
n3=nx_2dy
n3d=nx_2dy



do iproc=0,mpidims(2)-1  ! loop over each slab
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
      call mpi_cart_rank(comm_3d,dest_pe3,dest_pe,ierr)
      ASSERT("transpose_to_y: MPI_cart_rank failure 1",ierr==0)



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
      tag=my_y
      call MPI_IRecv(recbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(1),ierr)
      ASSERT("transpose_to_y: MPI_IRecv failure 1",ierr==0)
      tag=iproc
      call MPI_ISend(sendbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(2),ierr)
      ASSERT("transpose_to_y: MPI_ISend failure 1",ierr==0)
      call MPI_waitall(2,request,statuses,ierr) 	
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
use params

real*8 p(nx,ny,nz)
real*8 pt(g_ny2,nslabz,nx_2dy)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc
integer i,j,k,ii,l
#ifdef USE_MPI
!real*8 sendbuf(nslaby*nslabz*nx_2dy)
!real*8 recbuf(nslaby*nslabz*nx_2dy)
integer ierr,dest_pe,request(2),statuses(MPI_STATUS_SIZE,2)
integer dest_pe3(3),tag
#endif

call wallclock(tmx1)

!
! each cube is broken, along the x axis, into mpidims(2) slabs of
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



do iproc=0,mpidims(2)-1  ! loop over each slab
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
      call mpi_cart_rank(comm_3d,dest_pe3,dest_pe,ierr)
      ASSERT("transpose_from_y: MPI_cart_rank failure 1",ierr==0)


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
      tag=my_y
      call MPI_IRecv(recbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(1),ierr)
      ASSERT("transpose_from_y: MPI_IRecv failure 1",ierr==0)
      tag=iproc
      call MPI_ISend(sendbuf,l,MPI_REAL8,dest_pe,tag,comm_3d,request(2),ierr)
      ASSERT("transpose_from_y: MPI_ISend failure 1",ierr==0)
      call MPI_waitall(2,request,statuses,ierr) 	
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









subroutine output1(p,pt,buf,fid,fpe)
use params
use mpi
CPOINTER fid
integer :: fpe           ! cpu to do the I/O
real*8 :: p(nx,ny,nz)
real*8 :: pt(g_nx2,nslabz,ny_2dx)
real*8 :: buf(o_nx,ny_2dx)

! local vars
real*8 saved_edge(o_nx)
integer sending_pe,ierr,tag,z_pe,y_pe,x_pe
#ifdef USE_MPI
integer request,statuses(MPI_STATUS_SIZE)
#endif
integer i,j,k,l,extra_k,kuse,dest_pe3(3)
integer n1,n1d,n2,n2d,n3,n3d
integer :: ny_2dx_actual 
ny_2dx_actual = ny_2dx

call transpose_to_x(p,pt,n1,n1d,n2,n2d,n3,n3d)
ASSERT("output1 dimension failure 2",n1==g_nx)
ASSERT("output1 dimension failure 3",n1d==g_nx2)
ASSERT("output1 dimension failure 4",n2==nslabz)
ASSERT("output1 dimension failure 5",n2d==nslabz)
ASSERT("output1 dimension failure 6",n3==ny_2dx)
ASSERT("output1 dimension failure 7",n3d==ny_2dx)


if (o_nz>g_nz .and. g_bdy_z1/=PERIODIC) then
   call abort("output1: cannot handle o_nz=g_nz+1 with non-periodic b.c.")
endif

do z_pe=0,ncpu_z-1
extra_k=0
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
   call mpi_cart_rank(comm_3d,dest_pe3,sending_pe,ierr)
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
         call MPI_ISend(buf,l,MPI_REAL8,fpe,tag,comm_3d,request,ierr)
         ASSERT("output1: MPI_ISend failure",ierr==0)
         call MPI_waitall(1,request,statuses,ierr) 	
         ASSERT("output1: MPI_waitalll failure",ierr==0)
#endif
      endif
   endif

   if (my_pe==fpe) then
      if (sending_pe==my_pe) then
         ! dont recieve message from self
      else
#ifdef USE_MPI
         call MPI_IRecv(buf,l,MPI_REAL8,sending_pe,tag,comm_3d,request,ierr)
         ASSERT("output1: MPI_IRecv failure",ierr==0)
         call MPI_waitall(1,request,statuses,ierr) 	
         ASSERT("output1: MPI_waitalll failure",ierr==0)
#endif
      endif
      
      if (o_nx>g_nx) then
         if (g_bdy_x1==PERIODIC) then
            buf(o_nx,:)=buf(1,:)  ! append to the end, x-direction
         else
            buf(o_nx,:)=0
         endif
      endif
      call cwrite8(fid,buf,o_nx*ny_2dx_actual)

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
            call cwrite8(fid,saved_edge,o_nx)
         endif
      endif

   endif
enddo
enddo
enddo
enddo

end subroutine








subroutine output1_spec(p,pt,buf,fid,fpe)
use params
use mpi

CPOINTER fid
integer :: fpe           ! cpu to do the I/O
real*8 :: p(nx,ny,nz)
real*8 :: pt(g_nx2,nslabz,ny_2dx)
real*8 :: buf(dealias_nx,ny_2dx)

! local vars
integer sending_pe,ierr,tag,z_pe,y_pe,x_pe
#ifdef USE_MPI
integer request(2),statuses(MPI_STATUS_SIZE,2)
#endif
integer i,j,k,l,extra_k,dest_pe3(3)
integer n1,n1d,n2,n2d,n3,n3d,jj



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
   call mpi_cart_rank(comm_3d,dest_pe3,sending_pe,ierr)
#else
   sending_pe=my_pe
#endif

   if (sending_pe==my_pe) then

      ! output pt(1:g_nx,k,1:ny_2dx) from cpus: x_pe,y_pe,z_pe
      jj=0
      do j=1,ny_2dx
         jj=j
         if (x_jmcord(j)>(g_ny/3)  .and. x_jmcord(j)/=(g_ny/2)) exit
      enddo
      l=dealias_nx*jj
      if (l>0) buf(1:dealias_nx,1:jj)=pt(1:dealias_nx,k,1:jj)

      if (my_pe == fpe) then
         ! dont send message to self
      else
#ifdef USE_MPI
         tag=1
         call MPI_ISend(l,1,MPI_INTEGER,fpe,tag,comm_3d,request(2),ierr)
         ASSERT("output1: MPI_ISend failure",ierr==0)

         if (l==0) then
            call MPI_waitall(1,request,statuses,ierr) 	
            ASSERT("output1: MPI_waitalll failure",ierr==0)
         else
            call MPI_ISend(buf,l,MPI_REAL8,fpe,tag,comm_3d,request(1),ierr)
            ASSERT("output1: MPI_ISend failure",ierr==0)
            call MPI_waitall(2,request,statuses,ierr) 	
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
         call MPI_IRecv(l,1,MPI_INTEGER,sending_pe,tag,comm_3d,request,ierr)
         ASSERT("output1: MPI_IRecv failure",ierr==0)
         call MPI_waitall(1,request,statuses,ierr) 	
         ASSERT("output1: MPI_waitalll failure",ierr==0)
         if (l>0) then
            call MPI_IRecv(buf,l,MPI_REAL8,sending_pe,tag,comm_3d,request,ierr)
            ASSERT("output1: MPI_IRecv failure",ierr==0)
            call MPI_waitall(1,request,statuses,ierr) 	
            ASSERT("output1: MPI_waitalll failure",ierr==0)
         endif
#endif
      endif
      
      if (l>0) call cwrite8(fid,buf,l)
   endif
enddo
enddo
endif
enddo
enddo

end subroutine






subroutine input1_spec(p,pt,buf,fid,fpe)
use params
use mpi

CPOINTER fid
integer :: fpe           ! cpu to do the I/O
real*8 :: p(nx,ny,nz)
real*8 :: pt(g_nx2,nslabz,ny_2dx)
real*8 :: buf(dealias_nx,ny_2dx)

! local vars
integer sending_pe,ierr,tag,z_pe,y_pe,x_pe
#ifdef USE_MPI
integer request(2),statuses(MPI_STATUS_SIZE,2)
#endif
integer i,j,k,l,extra_k,dest_pe3(3)
integer n1,n1d,n2,n2d,n3,n3d,jj


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
   call mpi_cart_rank(comm_3d,dest_pe3,sending_pe,ierr)
#else
   sending_pe=my_pe
#endif

   if (sending_pe==my_pe) then

      ! output pt(1:g_nx,k,1:ny_2dx) from cpus: x_pe,y_pe,z_pe
      jj=0
      do j=1,ny_2dx
         jj=j
         if (x_jmcord(j)>(g_ny/3)  .and. x_jmcord(j)/=(g_ny/2)) exit
      enddo
      l=dealias_nx*jj

      if (my_pe == fpe) then
         ! dont send message to self
      else
#ifdef USE_MPI
         tag=1
         call MPI_ISend(l,1,MPI_INTEGER,fpe,tag,comm_3d,request,ierr)
         ASSERT("output1: MPI_ISend failure",ierr==0)
         call MPI_waitall(1,request,statuses,ierr) 	
         ASSERT("output1: MPI_waitalll failure",ierr==0)

         if (l>0) then
            call MPI_IRecv(buf,l,MPI_REAL8,fpe,tag,comm_3d,request,ierr)
            ASSERT("output1: MPI_ISend failure",ierr==0)
            call MPI_waitall(1,request,statuses,ierr) 	
            ASSERT("output1: MPI_waitalll failure",ierr==0)
         endif
#endif
      endif
   endif

   if (my_pe==fpe) then
      if (sending_pe==my_pe) then
         ! dont recieve message from self
         if (l>0) call cread8(fid,buf,l)
      else
#ifdef USE_MPI
         call MPI_IRecv(l,1,MPI_INTEGER,sending_pe,tag,comm_3d,request,ierr)
         ASSERT("output1: MPI_IRecv failure",ierr==0)
         call MPI_waitall(1,request,statuses,ierr) 	
         ASSERT("output1: MPI_waitalll failure",ierr==0)

         if (l>0) then
            call cread8(fid,buf,l)
            call MPI_ISend(buf,l,MPI_REAL8,sending_pe,tag,comm_3d,request,ierr)
            ASSERT("output1: MPI_IRecv failure",ierr==0)
            call MPI_waitall(1,request,statuses,ierr) 	
            ASSERT("output1: MPI_waitalll failure",ierr==0)
         endif
#endif
      endif
   endif

   if (l>0) then 
      pt(1:dealias_nx,k,1:jj)=buf(1:dealias_nx,1:jj)
   endif

enddo
enddo
endif
enddo
enddo


call transpose_from_x(pt,p,n1,n1d,n2,n2d,n3,n3d)


end subroutine







!
!  parallel input from a single file
!  if random==.true., then fill p with random numbers between -1.. 1
!  This is used to get random I.C. that are independent of the
!  parallel decomposition
!
!
subroutine input1(p,pt,buf,fid,fpe,random)
use params
use mpi
integer :: fpe         ! cpu to do all the file I/O

CPOINTER fid
real*8 :: p(nx,ny,nz)
real*8 :: pt(g_nx2,nslabz,ny_2dx)
real*8 :: buf(o_nx,ny_2dx)

! local vars
real*8 saved_edge(o_nx)
integer destination_pe,ierr,tag,z_pe,y_pe,x_pe
logical :: random


#ifdef USE_MPI
integer request,statuses(MPI_STATUS_SIZE)
#endif
integer i,j,k,l,extra_k,kuse,dest_pe3(3)
integer n1,n1d,n2,n2d,n3,n3d
integer :: ny_2dx_actual 
ny_2dx_actual = ny_2dx

if (o_nz>g_nz .and. g_bdy_z1/=PERIODIC) then
   call abort("output1: cannot handle o_nz=g_nz+1 with non-periodic b.c.")
endif


do z_pe=0,ncpu_z-1
extra_k=0
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
   call mpi_cart_rank(comm_3d,dest_pe3,destination_pe,ierr)
#else
   destination_pe=my_pe
#endif


   if (my_pe==fpe) then
      if (random) then
         call random_data(buf,o_nx*ny_2dx_actual)
      else
         call cread8(fid,buf,o_nx*ny_2dx_actual)
      endif

      if (o_ny>g_ny) then
      if (y_pe==ncpu_y-1 .and. x_pe==ncpu_x-1) then    
         ! read and discard periodic duplicate points
         if (random) then
            call random_data(saved_edge,o_nx)
         else
            call cread8(fid,saved_edge,o_nx)      
         endif
      endif
      endif

   endif

   if (my_pe==fpe) then
      if (destination_pe==my_pe) then
         ! dont send message to self
      else
#ifdef USE_MPI
         call MPI_ISend(buf,l,MPI_REAL8,destination_pe,tag,comm_3d,request,ierr)
         ASSERT("output1: MPI_IRecv failure",ierr==0)
         call MPI_waitall(1,request,statuses,ierr) 	
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
         call MPI_IRecv(buf,l,MPI_REAL8,fpe,tag,comm_3d,request,ierr)
         ASSERT("output1: MPI_ISend failure",ierr==0)
         call MPI_waitall(1,request,statuses,ierr) 	
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











end module
