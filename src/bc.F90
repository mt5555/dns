#include "macros.h"


module bc
use params
implicit none


! boundary condition on PSI.  
! used by bc_biotsavart.  
! if runbs=1,   psi_b is computed using biot savart law
!    runbs=0,   boundary data psi_b copied into psi
!
! you can also use the routine set_biotsavart() to set psi_b
! based on the boundary values in an array PSI
!
integer,parameter,private ::  bsize=max(1+g_nx,1+g_ny)
real*8,private :: psi_b(bsize,2,2)   ! psi_b(k,1,1) = x1 boundary
                                              ! psi_b(k,1,2) = x2 boundary
                                              ! psi_b(k,2,1) = y1 boundary
                                              ! psi_b(k,2,2) = y2 boundary
                                              ! psi_b(k,2,2) = y2 boundary
real*8 :: psi_b_temp(bsize,2,2)
real*8 :: biotsavart_cutoff
real*8 :: biotsavart_ubar=0
integer :: biotsavart_apply=5   !apply biotsavart every 5th timestep

contains



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



if (REALBOUNDARY(bdy_x1)) then
   !           nx1-1     nx1   nx1+1   nx1+2    nx1+3
   ! stencil ( 1/12     -2/3      0     2/3     -1/12)     /h
   !                    -1/2            1/2                /h
   ! multiply both sided by 12h:
   !           nx1-1     nx1   nx1+1   nx1+2    nx1+3
   ! stencil    1       -8      0       8       -1
   !                    -6              6         

   do j=by1,by2
      w(bx1-1,j)= 2*w(bx1,j)  - 2*w(bx1+2,j) +  w(bx1+3,j)
   enddo
endif
if (REALBOUNDARY(bdy_x2)) then
   !           bx2-3   nx2-2   nx2-1   nx2     nx2+1
   ! stencil    1       -8      0       8       -1
   !                    -6              6         
   do j=by1,by2
      w(bx2+1,j)= 2*w(bx2,j)  -2*w(bx2-2,j)  +  w(bx2-3,j)
   enddo
endif


if (REALBOUNDARY(bdy_y1)) then
   do i=bx1,bx2
      w(i,by1-1)= 2*w(i,by1)  - 2*w(i,by1+2) +  w(i,by1+3)
   enddo
endif

if (REALBOUNDARY(bdy_y2)) then
   do i=bx1,bx2
      w(i,by2+1)=  2*w(i,by2)  -2*w(i,by2-2)  +  w(i,by2-3)
   enddo
endif



end subroutine








subroutine bcw_impose(w)
!
! on non-periodic or non-reflective boundarys:
!
! set w=0 on boundary for inflow
! interpolate for outflow.
!
! set boundary so that centered differences one away are the same as
! onesided differnces at that point.  
! we get the same anser if we do this for: 
!    2nd order d/dx centered = 2nd order d/dx onesided 
!    2nd order dd/dxx centered = 1st order dd/dxx onsided
!    2nd order extrapolation
!    
!  
use params
implicit none
real*8 w(nx,ny)
!real*8 psi(nx,ny)


!local
integer i,j
real*8 :: u,v

if REALBOUNDARY(bdy_x1) then
   !             nx1   nx1+1   nx1+2    nx1+3
   ! stencil (   -1             1              ) /2h
   !                    -3      4       -1     ) /2h
   !      
   !
   !         -(nx1) + (nx1+2)  = -3(nx1+1) + 4(nx1+2) + (nx1+3)
   !          nx1 = 3(nx1+1) - 3(nx1+2) - (nx1+3)
   do j=by1,by2
      
      if (bdy_x1==INFLOW0_ONESIDED) then
         w(bx1,j)= 3*w(bx1+1,j)-3*w(bx1+2,j)+w(bx1+3,j)
      else
         w(bx1,j)=0
      endif
   enddo
endif
if REALBOUNDARY(bdy_x2) then
   do j=by1,by2

      if (bdy_x2==INFLOW0_ONESIDED) then
         w(bx2,j)=0
      else
         w(bx2,j)=0
      endif
   enddo
endif


if REALBOUNDARY(bdy_y1) then
   do i=bx1,bx2

      if (bdy_y1==INFLOW0_ONESIDED) then
         stop 'should not happen'
      else
         w(i,by1)=0
      endif
   enddo
endif

if REALBOUNDARY(bdy_y2) then
   do i=bx1,bx2
      if (bdy_y2==INFLOW0_ONESIDED) then
         if (xcord(i)<=xscale/2) then
            w(i,by2)=0
         else
            w(i,by2)=  3*w(i,by2-1) - 3*w(i,by2-2) + w(i,by2-3)
         endif
      else
         w(i,by2)=0
      endif
   enddo
endif


call ghost_update_x_reshape(w,1)
call ghost_update_y_reshape(w,1)

end subroutine
















subroutine set_biotsavart(psi)
!
! on non-periodic or non-reflective boundarys:
! set the boundary data psi_b based on boundary data in psi.
!
use params
implicit none
real*8 psi(nx,ny)

integer :: i,j,k,l
! set PSI boundary data:

if (REALBOUNDARY(bdy_x1)) then
   do j=by1,by2
      l = j-by1+1 + nslaby*my_y
      psi_b(L,1,1) = psi(bx1,j) 
   enddo
endif


if (REALBOUNDARY(bdy_x2)) then
   do j=by1,by2
      l = j-by1+1 + nslaby*my_y
      psi_b(L,1,2) = psi(bx2,j) 
   enddo
endif

if (REALBOUNDARY(bdy_y1)) then
   do i=bx1,bx2
      l = i-bx1+1 + nslabx*my_x
      psi_b(L,2,1) = psi(i,by1) 
   enddo
endif


if (REALBOUNDARY(bdy_y2)) then
   do i=bx1,bx2
      l = i-bx1+1 + nslabx*my_x
      psi_b(L,2,2) = psi(i,by2) 
   enddo
endif
end subroutine      







subroutine bc_biotsavart(w,psi,runbs)
!
! on non-periodic or non-reflective boundarys:
! use biot-savar law to compute boundary data for PSI.
!
! (actually: hard coded to do x1,x2 and y2 boundaries only)
!
! except for y1 boundary, where we set PSI=0
!
!     Finds psi on boundary from w using Biot-Savart
!
!     psi = -(1/4pi)*sum w(i,j)*logterm(i,j, 0,k,w)*Delta A - ubar*y
!
!     is the streamfunction in a reference frame moving with the 
!     predetermined velocity ubar
!
use params
use mpi
implicit none
real*8 w(nx,ny)
real*8 psi(nx,ny)
integer :: runbs

! local
integer i,j,k,l,ierr,last,nskip
real*8 :: dela
real*8 tmx1,tmx2
logical interp

call wallclock(tmx1)

if (runbs==1) then
! init to zero on boundary
psi_b=0


! should we use interpolation to compute PSI?
interp=(biotsavart_cutoff>0)



if (interp) then
   nskip=10
   ! interpolate every nskip point
   do j=by1,by2
   if (mod(j,100)==0 .and. my_pe==io_pe) then
      print *,'interp: ',by2,j
   endif
   do i=bx1,bx2
   if (abs(w(i,j)).ge.biotsavart_cutoff) then
      do k=1,o_nx,nskip
         !psi_b(k,2,1) = psi_b(k,2,1) - w(i,j)*logterm(i,j,k,1)
         psi_b(k,2,2) = psi_b(k,2,2) - w(i,j)*logterm(i,j,k,o_ny)
      enddo
      last=k-nskip
      do k=last+1,o_nx
         !psi_b(k,2,1) = psi_b(k,2,1) - w(i,j)*logterm(i,j,k,1)
         psi_b(k,2,2) = psi_b(k,2,2) - w(i,j)*logterm(i,j,k,o_ny)
      enddo
      do k=11,o_ny,nskip
            psi_b(k,1,1) = psi_b(k,1,1) - w(i,j)*logterm(i,j,1,k)
            psi_b(k,1,2) = psi_b(k,1,2) - w(i,j)*logterm(i,j,o_nx,k)
      enddo
      last=k-nskip
      do k=last+1,o_ny
            psi_b(k,1,1) = psi_b(k,1,1) - w(i,j)*logterm(i,j,1,k)
            psi_b(k,1,2) = psi_b(k,1,2) - w(i,j)*logterm(i,j,o_nx,k)
      enddo
   endif
   enddo
   enddo

else
   do j=by1,by2
   do i=bx1,bx2
   if (abs(w(i,j)).ge.biotsavart_cutoff) then
      do k=1,o_nx
         !psi_b(k,2,1) = psi_b(k,2,1) - w(i,j)*logterm(i,j,k,1)
         psi_b(k,2,2) = psi_b(k,2,2) - w(i,j)*logterm(i,j,k,o_ny)
      enddo
      do k=2,o_ny 
         psi_b(k,1,1) = psi_b(k,1,1) - w(i,j)*logterm(i,j,1,k)
         psi_b(k,1,2) = psi_b(k,1,2) - w(i,j)*logterm(i,j,o_nx,k)
      enddo
!      print *,i,j,w(i,j)
   endif
   enddo
enddo
endif



#ifdef USE_MPI
psi_b_temp=psi_b
call MPI_allreduce(psi_b_temp,psi_b,4*bsize,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif

!print *,'psib',psi_b(o_ny,1,1)


if (interp) call intpsi(psi_b,nskip,bsize)


! add ubar correction
! only needs to be done on processer which owns boundary data
! but lets have everyone do it for now...
dela = delx*dely/(4*pi)
do k=1,o_nx
   ! y2 boundary
   !psi_b(k,2,1) = psi_b(k,2,1)*dela - biotsavart_ubar*g_ycord(1)
   psi_b(k,2,2) = psi_b(k,2,2)*dela - biotsavart_ubar*g_ycord(o_ny)
enddo
do k=2,o_ny
   ! x1,x2 boundary
   psi_b(k,1,1) = psi_b(k,1,1)*dela  - biotsavart_ubar*g_ycord(k)
   psi_b(k,1,2) = psi_b(k,1,2)*dela  - biotsavart_ubar*g_ycord(k)
enddo
endif




! set PSI boundary data:

if (REALBOUNDARY(bdy_x1)) then
   do j=by1,by2
      l = j-by1+1 + nslaby*my_y
      psi(bx1,j)= psi_b(L,1,1)
   enddo
endif
if (REALBOUNDARY(bdy_x2)) then
   do j=by1,by2
      l = j-by1+1 + nslaby*my_y
      psi(bx2,j)= psi_b(L,1,2)
   enddo
endif
if (REALBOUNDARY(bdy_y1)) then
   do i=bx1,bx2
      l = i-bx1+1 + nslabx*my_x
      psi(i,by1)= psi_b(L,2,1)
   enddo
endif
if (REALBOUNDARY(bdy_y2)) then
   do i=bx1,bx2
      l = i-bx1+1 + nslabx*my_x
      psi(i,by2)= psi_b(L,2,2)
   enddo
endif
      
call wallclock(tmx2)
tims(14)=tims(14)+(tmx2-tmx1)


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
if (denom1==0) then
   logterm=0
else
   logterm = log(denom1/denom2)
endif
end function 





subroutine intpsi(psi,nskip,n)
use params
implicit none
integer :: n,nskip
real*8 psi(n,2,2)

integer i,k

! x1 and x2 boundary
! psi(i,1,*) was set for i=1,o_nx,nskip

do k=1,2
do i=1,o_ny,nskip
   if (i+3*nskip>o_ny) exit
   if (i==1) then  
      ! do the [y0,y1] piece
      call interpn(psi(i,1,k),psi(i+nskip,1,k),psi(i+2*nskip,1,k),psi(i+3*nskip,1,k),psi(i+1,1,k),nskip,0)
   endif
   if (i+4*nskip>o_ny) then  
      ! do the [y2,y3] piece, this is the final iteration
      call interpn(psi(i,1,k),psi(i+nskip,1,k),psi(i+2*nskip,1,k),psi(i+3*nskip,1,k),psi(i+1+2*nskip,1,k),nskip,2)
   endif
   ! using [y0,y1,y2,y3], interpolate the [y1,y2] piece
   call interpn(psi(i,1,k),psi(i+nskip,1,k),psi(i+2*nskip,1,k),psi(i+3*nskip,1,k),psi(i+1+nskip,1,k),nskip,1)
enddo
enddo


! y2 boundary
! psi(i,2,2) was set for i=1,o_ny,10
k=2
do i=1,o_nx,nskip
   if (i+3*nskip>o_nx) exit
   if (i==1) then  
      ! do the [y0,y1] piece
      call interpn(psi(i,2,k),psi(i+nskip,2,k),psi(i+2*nskip,2,k),psi(i+3*nskip,2,k),psi(i+1,2,k),nskip,0)
   endif
   if (i+4*nskip>o_nx) then  
      ! do the [y2,y3] piece, this is the final iteration
      call interpn(psi(i,2,k),psi(i+nskip,2,k),psi(i+2*nskip,2,k),psi(i+3*nskip,2,k),psi(i+1+2*nskip,2,k),nskip,2)
   endif
   ! using [y0,y1,y2,y3], interpolate the [y1,y2] piece
   call interpn(psi(i,2,k),psi(i+nskip,2,k),psi(i+2*nskip,2,k),psi(i+3*nskip,2,k),psi(i+1+nskip,2,k),nskip,1)
enddo
end subroutine









end module

      SUBROUTINE interpn(y0,y1,y2,y3,ynew,n,xbeg)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Interpolates "y" to "n-1" points in the interval 
!c     [xbeg,xbeg+1] of the 3 intervals [0,1],[1,2],[2,3] 
!c     --- y0,y1,y2,y3 are assumed to be evenly spaced
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer i,n,xbeg
      real*8 x0,x1,x2,x3,y0,y1,y2,y3,ynew(*),delx
      real*8 denom0,denom1,denom2,denom3,fact0,fact1,fact2,fact3,newx

      x0 = 0
      x1 = 1
      x2 = 2
      x3 = 3

      denom0 = -6  !(x0-x1)*(x0-x2)*(x0-x3)   ! (-1)(-2)(-3)=-6
      denom1 =  2  !(x1-x0)*(x1-x2)*(x1-x3)   ! ( 1)(-1)(-2)= 2
      denom2 = -2  !(x2-x0)*(x2-x1)*(x2-x3)   ! ( 2)( 1)(-1)=-2
      denom3 =  6  !(x3-x0)*(x3-x1)*(x3-x2)   ! ( 3)( 2)( 1)= 6

      delx = (x1-x0)/n
      do i=1,n-1
        newx = xbeg + i*delx
        
        fact0 = (newx-x1)*(newx-x2)*(newx-x3)/denom0
        fact1 = (newx-x0)*(newx-x2)*(newx-x3)/denom1
        fact2 = (newx-x0)*(newx-x1)*(newx-x3)/denom2
        fact3 = (newx-x0)*(newx-x1)*(newx-x2)/denom3

        ynew(i) = y0*fact0 + y1*fact1 + y2*fact2 + y3*fact3
      enddo

      return
      end subroutine







subroutine ghost_update_x_reshape(psi,n)
use params
use ghost
implicit none
integer :: n
real*8 :: psi(nx,ny,nz,n)
call ghost_update_x(psi,n)
end
subroutine ghost_update_y_reshape(psi,n)
use params
use ghost
implicit none
integer :: n
real*8 :: psi(nx,ny,nz,n)
call ghost_update_y(psi,n)
end


