#include "macros.h"
subroutine bc_preloop
use params
implicit none



end subroutine




subroutine bc_onesided(w)
! on non-periodic or non-reflective boundarys:
!
! fiddle first ghost cell so that 4th order derivative at 1st interior point
! will be the same as a 2nd order centered scheme.  
! (assuming boundary values already set)

use params
implicit none
real*8 :: w(nx,ny)


!local
integer i,j
integer :: bx1,bx2,by1,by2
bx1=nx1
bx2=nx2
by1=ny1
by2=ny2

if (my_x==0 .and. bdy_x1==INFLOW0_ONESIDED) then
   !           nx1-1     nx1   nx1+1   nx1+2    nx1+3
   ! stencil ( 1/12     -2/3      0     2/3     -1/12)     /h
   !                    -1/2            1/2                /h
   ! multiply both sided by 12h:
   !           nx1-1     nx1   nx1+1   nx1+2    nx1+3
   ! stencil    1       -8      0       8       -1
   !                    -6              6         

   do j=ny1,ny2
      w(bx1-1,j)= 2*w(bx1,j)  - 2*w(bx1+2,j) +  w(bx1+3,j)
   enddo
endif
if (my_x==ncpu_x-1 .and. bdy_x2==INFLOW0_ONESIDED) then
   !           nx2-3   nx2-2   nx2-1   nx2     nx2+1
   ! stencil    1       -8      0       8       -1
   !                    -6              6         
   do j=ny1,ny2
      w(bx2+1,j)= 2*w(bx2,j)  -2*w(bx2-2,j)  +  w(bx2-3,j)
   enddo
endif


if (my_y==0 .and. bdy_y1==INFLOW0_ONESIDED) then
   do i=nx1,nx2
      w(i,by1-1)= 2*w(i,by1)  - 2*w(i,by1+2) +  w(i,by1+3)
   enddo
endif

if (my_y==ncpu_y-1 .and. bdy_y2==INFLOW0_ONESIDED) then
   do i=nx1,nx2
      w(i,by2+1)=  2*w(i,by2)  -2*w(i,by2-2)  +  w(i,by2-3)
   enddo
endif



end subroutine




subroutine bc_biotsavart(w,psi)
! on non-periodic or non-reflective boundarys:
! use biot-savar law to compute boundary data for PSI.
!
! except for y1 boundary, where we set PSI=0
!
!     Finds psi on boundary from w using Biot-Savart
!
!     psi = -(1/4pi)*sum w(i,j)*logterm(i,j, 0,k,w)*Delta A - ubar*y
!
!     is the streamfunction in a reference frame moving with the 
!     predetermined velocity ubar
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

use params
use mpi
implicit none
real*8 w(nx,ny)
real*8 psi(nx,ny)

! local
integer i,j,k,l,ierr
real*8 :: dela
real*8,external :: logterm
real*8 :: ubar=0
real*8 :: eps=1e-8
real*8 :: psi_b(max(g_ny,g_nx),2,2)   ! psi_b(k,1,1) = x1 boundary
                                      ! psi_b(k,1,2) = x2 boundary
                                      ! psi_b(k,2,1) = y1 boundary
                                      ! psi_b(k,2,2) = y2 boundary
real*8 :: temp(max(g_ny,g_nx),2,2)


! init to zero on boundary
psi_b=0

do j=ny1,ny2
do i=nx1,nx2
   if (abs(w(i,j)).ge.eps) then
      do k=1,g_nx  !,10
         psi_b(k,2,2) = psi_b(k,2,2) - w(i,j)*logterm(i,j,k,g_ny)
      enddo
      do k=2,g_ny  !,10
         psi_b(k,1,1) = psi_b(k,1,1) - w(i,j)*logterm(i,j,1,k)
         psi_b(k,1,2) = psi_b(k,1,2) - w(i,j)*logterm(i,j,g_nx,k)
      enddo
   endif
enddo
enddo


#if 0
C     INTERPOLATE TO INTERMEDIATE POINTS
      call intpsi(psi)
#endif


#ifdef USE_MPI
temp=psi_b
k=max(g_ny,g_nx)*4
call MPI_allreduce(temp,psi_b,k,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif


! add ubar correction
! only needs to be done on processer which owns boundary data
! but lets have everyone do it for now...
dela = delx*dely/(4*pi)
do k=1,g_nx
   ! y2 boundary
   psi_b(k,2,2) = psi_b(k,2,2)*dela - ubar*g_ycord(g_ny)
enddo
do k=2,g_ny
   ! x1,x2 boundary
   psi_b(k,1,1) = psi_b(k,1,1)*dela  - ubar*g_ycord(k)
   psi_b(k,1,2) = psi_b(k,1,2)*dela  - ubar*g_ycord(k)
enddo


! set PSI boundary data:

if (my_x==0 .and. bdy_x1==INFLOW0_ONESIDED) then
   do j=ny1,ny2
      l = j-ny1+1 + nslaby*my_y
      psi(nx1,j)= psi_b(L,1,1)
   enddo
endif
if (my_x==ncpu_x-1 .and. bdy_x2==INFLOW0_ONESIDED) then
   do j=ny1,ny2
      l = j-ny1+1 + nslaby*my_y
      psi(nx2,j)= psi_b(L,1,2)
   enddo
endif

if (my_y==0 .and. bdy_y1==INFLOW0_ONESIDED) then
   do i=nx1,nx2
      l = i-nx1+1 + nslabx*my_x
      psi(i,ny1)= psi_b(L,2,1)
   enddo
endif

if (my_y==ncpu_y-1   .and. bdy_y2==INFLOW0_ONESIDED) then
   do i=nx1,nx2
      l = i-nx1+1 + nslabx*my_x
      psi(i,ny2)= psi_b(L,2,2)
   enddo
endif






return
end subroutine



real*8 FUNCTION logterm(i,j,k,l)
!
!     Finds streamfunction at (k,l) induced by filament pair at (i,j)
!
use params
implicit none
integer i,j,k,l
real*8 difx,dify,sumy,denom1,denom2

difx = xcord(i) - g_xcord(k)
dify = ycord(j) - g_ycord(l)
sumy = ycord(j) + g_ycord(l)

denom1 = difx**2 + dify**2
denom2 = difx**2 + sumy**2

logterm = log(denom1/denom2)
end function 


