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
   divmx=max(divmx,p(i,j,k))
   divi=divi+p(i,j,k)
enddo
enddo
enddo
divi=divi/g_nx/g_ny/g_nz

#ifdef USE_MPI
   g_mx=divmx
   call MPI_allreduce(g_mx,divmx,1,MPI_REAL8,MPI_MAX,comm_3d,ierr)
   g_mx=divi
   call MPI_allreduce(g_mx,divi,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
end subroutine











