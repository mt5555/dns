#include "macros.h"

module transform
use params
implicit none


contains


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
real*8 pt(g_nz2,nslabx,ny_2d)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc
integer i,j,k,jj,l

!
! each cube is broken, along the y axis, into ncpu_z slabs of
! size ny_2d = (ny2-ny1+1)/ncpu_z
!
! in the z direction, the dimension is nslabz = (nz2-nz1+1)
!

n1=g_nz
n1d=g_nz2   	
n2=nslabx
n2d=nslabx
n3=ny_2d
n3d=ny_2d



do iproc=0,ncpu_z-1  ! loop over each slab
!   if (iproc==myproc_z) then
    if (.true.) then
      do j=1,ny_2d  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2d +j -1
         ASSERT("transpose_to_z jj failure 1",jj<=ny2)
         ASSERT("transpose_to_z jj failure 2",jj>=ny1)
         do k=nz1,nz2
         do i=nx1,nx2
            pt(k+iproc*nslabz-nz1+1,i-nx1+1,jj)=p(i,jj,k)
         enddo
         enddo
      enddo
   else
#ifdef MPI
      l=0
      do j=1,ny_2d  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2d +j -1
         ASSERT("transpose_to_z jj failure 1",jj<=ny2)
         ASSERT("transpose_to_z jj failure 2",jj>=ny1)
         do k=nz1,nz2
         do i=nx1,nx2
            l=l+1
            sendbuf(l)=p(i,jj,k)
         enddo
         enddo
      enddo

!     send buffer to (myproc_x,iproc,mproc_z)
!     rec  buffer from (myproc_x,iproc,mproc_z)

      l=0
      do j=1,ny_2d  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2d +j -1
         ASSERT("transpose_to_z jj failure 1",jj<=ny2)
         ASSERT("transpose_to_z jj failure 2",jj>=ny1)
         do k=nz1,nz2
         do i=nx1,nx2
	    l=l+1
            pt(k+iproc*nslabz-nz1+1,i-nx1+1,jj)=recbuf(l)
         enddo
         enddo
      enddo
#endif
   endif
enddo
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
use params
real*8 p(nx,ny,nz)
real*8 pt(g_nz2,nslabx,ny_2d)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc
integer i,j,k,jj,l

!
! each cube is broken, along the y axis, into ncpu_z slabs of
! size ny_2d = (ny2-ny1+1)/ncpu_z
!
! in the z direction, the dimension is nslabz = (nz2-nz1+1)



! If any of these fail, then pt was probably not computed
! via a call to transpose_to_z().
ASSERT("transpose_from_x dimension failure 2",n1==g_nz)
ASSERT("transpose_from_x dimension failure 3",n1d==g_nz2)
ASSERT("transpose_from_x dimension failure 4",n2==nslabx)
ASSERT("transpose_from_x dimension failure 5",n2d==nslabx)
ASSERT("transpose_from_x dimension failure 6",n3==ny_2d)
ASSERT("transpose_from_x dimension failure 7",n3d==ny_2d)



do iproc=0,ncpu_z-1  ! loop over each slab
!   if (iproc==myproc_z) then
    if (.true.) then
      do j=1,ny_2d  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2d +j -1
         ASSERT("transpose_from_z jj failure 1",jj<=ny2)
         ASSERT("transpose_from_z jj failure 2",jj>=ny1)
         do k=nz1,nz2
         do i=nx1,nx2
            p(i,jj,k)=pt(k+iproc*nslabz-nz1+1,i-nx1+1,jj)
         enddo
         enddo
      enddo
   else
#ifdef MPI
      l=0
      do j=1,ny_2d  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2d +j -1
         ASSERT("transpose_from_z jj failure 1",jj<=ny2)
         ASSERT("transpose_from_z jj failure 2",jj>=ny1)
         do k=nz1,nz2
         do i=nx1,nx2
            l=l+1
            sendbuf(l)=pt(k+iproc*nslabz-nz1+1,i-nx1+1,jj)
         enddo
         enddo
      enddo

!     send buffer to (myproc_x,iproc,mproc_z)
!     rec  buffer from (myproc_x,iproc,mproc_z)

      l=0
      do j=1,ny_2d  ! loop over points in a single slab
         jj=ny1 + iproc*ny_2d +j -1
         ASSERT("transpose_from_z jj failure 1",jj<=ny2)
         ASSERT("transpose_from_z jj failure 2",jj>=ny1)
         do k=nz1,nz2
         do i=nx1,nx2
	    l=l+1
            p(i,jj,k)=recbuf(l)
         enddo
         enddo
      enddo
#endif
   endif
enddo
end subroutine




























!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! transform data in a cartesian decomposition into a 2D decomposition
! with x as the leading index.
!
! input: p
! ouput: pt, and the dimensions of pt: n1,n1d,n2,n2d,n3,n3d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_to_x(p,pt,n1,n1d,n2,n2d,n3,n3d)
use params

real*8 p(nx,ny,nz)
real*8 pt(g_nx2,nslaby,nz_2d)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc
integer i,j,k,kk,l



!
! each cube is broken, along the z axis, into ncpu_x slabs of
! size nz_2d
!
! in the x direction, the dimension is nslabx = (nx2-nx1+1)


n1=g_nx
n1d=g_nx2
n2=nslaby
n2d=nslaby
n3=nz_2d
n3d=nz_2d



do iproc=0,ncpu_x-1  ! loop over each slab
!   if (iproc==myproc_x) then
    if (.true.) then
      do k=1,nz_2d  ! loop over points in a single slab
         kk=nz1 + iproc*nz_2d +k -1
         ASSERT("transpose_to_x kk failure 1",kk<=nz2)
         ASSERT("transpose_to_x kk failure 2",kk>=nz1)
         do i=nx1,nx2
         do j=ny1,ny2
            pt(i+iproc*nslabx-nx1+1,j-ny1+1,kk)=p(i,j,kk)
         enddo
         enddo
      enddo
   else
#ifdef MPI
      l=0
      do k=1,nz_2d  ! loop over points in a single slab
         kk=nz1 + iproc*nz_2d +k -1
         ASSERT("transpose_to_x kk failure 1",kk<=nz2)
         ASSERT("transpose_to_x kk failure 2",kk>=nz1)
         do i=nx1,nx2
         do j=ny1,ny2
            l=l+1
            sendbuf(l)=p(i,j,kk)
         enddo
         enddo
      enddo

!     send buffer to (iproc,myproc_y,mproc_z)
!     rec  buffer from (iproc,myproc_y,mproc_z)

      l=0
      do k=1,nz_2d  ! loop over points in a single slab
         kk=nz1 + iproc*nz_2d +k -1
         ASSERT("transpose_to_x kk failure 1",kk<=nz2)
         ASSERT("transpose_to_x kk failure 2",kk>=nz1)
         do i=nx1,nx2
         do j=ny1,ny2
            l=l+1
            pt(i+iproc*nslabx-nx1+1,j-ny1+1,kk)=recbuf(l)
         enddo
         enddo
      enddo
#endif
   endif
enddo

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! transform data to a cartesian decomposition from a 2D decomposition
! with x as the leading index.
!
! input: pt, and the dimensions of pt: n1,n1d,n2,n2d,n3,n3d
! ouput: p 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine transpose_from_x(pt,p,n1,n1d,n2,n2d,n3,n3d)
use params

real*8 p(nx,ny,nz)
real*8 pt(g_nx2,nslaby,nz_2d)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc
integer i,j,k,kk,l

!
! each cube is broken, along the z axis, into ncpu_x slabs of
! size nz_2d
!
! in the x direction, the dimension is nslabx = (nx2-nx1+1)


! If any of these fail, then pt was probably not computed
! via a call to transpose_to_y().
ASSERT("transpose_from_x dimension failure 2",n1==g_nx)
ASSERT("transpose_from_x dimension failure 3",n1d==g_nx2)
ASSERT("transpose_from_x dimension failure 4",n2==nslaby)
ASSERT("transpose_from_x dimension failure 5",n2d==nslaby)
ASSERT("transpose_from_x dimension failure 6",n3==nz_2d)
ASSERT("transpose_from_x dimension failure 7",n3d==nz_2d)


do iproc=0,ncpu_x-1  ! loop over each slab
!   if (iproc==myproc_x) then
    if (.true.) then
      do k=1,nz_2d  ! loop over points in a single slab
         kk=nz1 + iproc*nz_2d +k -1
         ASSERT("transpose_to_x kk failure 1",kk<=nz2)
         ASSERT("transpose_to_x kk failure 2",kk>=nz1)
         do i=nx1,nx2
         do j=ny1,ny2
            p(i,j,kk)=pt(i+iproc*nslabx-nx1+1,j-ny1+1,kk)
         enddo
         enddo
      enddo
   else
#ifdef MPI
      l=0
      do k=1,nz_2d  ! loop over points in a single slab
         kk=nz1 + iproc*nz_2d +k -1
         ASSERT("transpose_to_x kk failure 1",kk<=nz2)
         ASSERT("transpose_to_x kk failure 2",kk>=nz1)
         do i=nx1,nx2
         do j=ny1,ny2
            l=l+1
            sendbuf(l)=pt(i+iproc*nslabx-nx1+1,j-ny1+1,kk)
         enddo
         enddo
      enddo

!     send buffer to (iproc,myproc_y,mproc_z)
!     rec  buffer from (iproc,myproc_y,mproc_z)

      l=0
      do k=1,nz_2d  ! loop over points in a single slab
         kk=nz1 + iproc*nz_2d +k -1
         ASSERT("transpose_to_x kk failure 1",kk<=nz2)
         ASSERT("transpose_to_x kk failure 2",kk>=nz1)
         do i=nx1,nx2
         do j=ny1,ny2
            l=l+1
            p(i,j,kk)=recbuf(l)
         enddo
         enddo
      enddo
#endif
   endif
enddo

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
real*8 pt(g_ny2,nslabz,nx_2d)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc
integer i,j,k,ii,l

!
! each cube is broken, along the x axis, into ncpu_y slabs of
! size nx_2d
!
! in the y direction, the dimension is nslaby = (ny2-ny1+1)


n1=g_ny
n1d=g_ny2
n2=nslabz
n2d=nslabz
n3=nx_2d
n3d=nx_2d



do iproc=0,ncpu_y-1  ! loop over each slab
!   if (iproc==myproc_y) then
    if (.true.) then
      do i=1,nx_2d  ! loop over points in a single slab
         ii=nx1 + iproc*nx_2d +i -1
         ASSERT("transpose_to_y ii failure 1",ii<=nx2)
         ASSERT("transpose_to_y ii failure 2",ii>=nx1)
         do j=ny1,ny2
         do k=nz1,nz2
            pt(j+iproc*nslaby-ny1+1,k-nz1+1,ii)=p(ii,j,k)
         enddo
         enddo
      enddo
   else
#ifdef MPI
      l=0
      do i=1,nx_2d  ! loop over points in a single slab
         ii=nx1 + iproc*nx_2d +i -1
         ASSERT("transpose_to_y ii failure 1",ii<=nx2)
         ASSERT("transpose_to_y ii failure 2",ii>=nx1)
         do j=ny1,ny2
         do k=nz1,nz2
            l=l+1
            sendbuf(l)=p(ii,j,k)
         enddo
         enddo
      enddo

!     send buffer to (myproc_x,iproc,mproc_z)
!     rec  buffer from (myproc_x,iproc,mproc_z)

      l=0
      do i=1,nx_2d  ! loop over points in a single slab
         ii=nx1 + iproc*nx_2d +i -1
         ASSERT("transpose_to_y ii failure 1",ii<=nx2)
         ASSERT("transpose_to_y ii failure 2",ii>=nx1)
         do j=ny1,ny2
         do k=nz1,nz2
            l=l+1
            pt(j+iproc*nslaby-ny1+1,k-nz1+1,ii)=recbuf(l)
         enddo
         enddo
      enddo
#endif
   endif
enddo

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
real*8 pt(g_ny2,nslabz,nx_2d)
integer n1,n1d,n2,n2d,n3,n3d

!local variables
integer iproc
integer i,j,k,ii,l

!
! each cube is broken, along the x axis, into ncpu_y slabs of
! size nx_2d
!
! in the y direction, the dimension is nslaby = (ny2-ny1+1)




! If any of these fail, then pt was probably not computed
! via a call to transpose_to_y().
ASSERT("transpose_from_y dimension failure 2",n1==g_ny)
ASSERT("transpose_from_y dimension failure 3",n1d==g_ny2)
ASSERT("transpose_from_y dimension failure 4",n2==(nz2-nz1+1))
ASSERT("transpose_from_y dimension failure 5",n2d==nz)
ASSERT("transpose_from_y dimension failure 6",n3==nx_2d)
ASSERT("transpose_from_y dimension failure 7",n3d==nx_2d)



do iproc=0,ncpu_y-1  ! loop over each slab
!   if (iproc==myproc_y) then
    if (.true.) then
      do i=1,nx_2d  ! loop over points in a single slab
         ii=nx1 + iproc*nx_2d +i-1
         ASSERT("transpose_from_y ii failure 1",ii<=nx2)
         ASSERT("transpose_from_y ii failure 2",ii>=nx1)
         do j=ny1,ny2
         do k=nz1,nz2
            p(ii,j,k)=pt(j+iproc*nslaby-ny1+1,k-nz1+1,ii)
         enddo
         enddo
      enddo
   else
#ifdef MPI
      l=0
      do i=1,nx_2d  ! loop over points in a single slab
         ii=nx1 + iproc*nx_2d +i -1
         ASSERT("transpose_from_y ii failure 1",ii<=nx2)
         ASSERT("transpose_from_y ii failure 2",ii>=nx1)
         do j=ny1,ny2
         do k=nz1,nz2
            l=l+1
            sendbuf(l)=pt(j+iproc*nslaby-ny1+1,k-nz1+1,ii)
         enddo
         enddo
      enddo

!     send buffer to (myproc_x,iproc,mproc_z)
!     rec  buffer from (myproc_x,iproc,mproc_z)

      l=0
      do i=1,nx_2d  ! loop over points in a single slab
         ii=nx1 + iproc*nx_2d + i -1
         ASSERT("transpose_from_y ii failure 1",ii<=nx2)
         ASSERT("transpose_from_y ii failure 2",ii>=nx1)
         do j=ny1,ny2
         do k=nz1,nz2
            l=l+1
            p(ii,j,k)=recbuf(l)
         enddo
         enddo
      enddo
#endif
   endif
enddo

end subroutine

end module

