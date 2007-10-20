#include "macros.h"







subroutine compute_div(Q,p,work1,work2,divmx,divi)
use params
use mpi
implicit none
integer pe
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: divmx,divi

!local variablesn
real*8 :: p(nx,ny,nz)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: g_mx
integer i,j,k,ierr

divmx=0
divi=0
call divergence(p,Q,work1,work2)
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   divmx=max(divmx,abs(p(i,j,k)))
   divi=divi+p(i,j,k)
enddo
enddo
enddo
divi=divi/g_nx/g_ny/g_nz

#ifdef USE_MPI
   g_mx=divmx
   call mpi_allreduce(g_mx,divmx,1,MPI_REAL8,MPI_MAX,comm_3d,ierr)
   g_mx=divi
   call mpi_allreduce(g_mx,divi,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
end subroutine



subroutine zaverage(Q,Qz)
!
!  store the z-average (repeated on every z-slab) of Q in Qz
!
use params
use mpi
implicit none
integer pe
real*8 :: Q(nx,ny,nz)
real*8 :: Qz(nx,ny,nz)

! local
integer i,j,k,ierr
real*8 :: work(nx,ny)


if (ndim/=3) then
   call abortdns("zaverage()  requires ndim==3")
endif
if (ncpu_x*ncpu_y>1) then
   call abortdns("zaverage()  requies ncpu_x=ncpu_y=1")
endif
if (nx2-nx1+1 /= g_nx) then
   call abortdns("zaverage()  x dimension error")
endif
if (ny2-ny1+1 /= g_ny) then
   call abortdns("zaverage()  y dimension error")
endif


! their z-ave and then reduce to io_pe
work=0
do k=nz1,nz2
   do j=ny1,ny2
      do i=nx1,nx2
         work(i,j)=work(i,j)+Q(i,j,k)/g_nz
      enddo
   enddo
enddo
#ifdef USE_MPI
Qz(:,:,1)=work
call mpi_allreduce(Qz,work,nx*ny,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
do k=nz1,nz2
   do j=ny1,ny2
      do i=nx1,nx2
         Qz(i,j,k)=work(i,j)
      enddo
   enddo
enddo
end subroutine zaverage









