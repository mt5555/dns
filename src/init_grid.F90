!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! initialization for fft grid
! domain is from 0..1
! first point is at grid point 0, last point is at gridpoint
! 1-delta_x
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_grid
use params
use fft_interface
implicit none

!local variables
real*8 :: one=1
integer i,j,k

! periodic case
o_nx=o_nx+1
o_ny=o_ny+1
if (nz>1) then
   o_nz=o_nz+1
else
   o_nz=1
endif


delx = one/g_nx
dely = one/g_ny
delz = one/g_nz

do i=1,o_nx
   g_xcord(i)=(i-1)*delx	
enddo
do j=1,o_ny
   g_ycord(j)=(j-1)*dely	
enddo
do k=1,o_nz
   g_zcord(k)=(k-1)*delz	
enddo

do i=nx1,nx2
   xcord(i)=(i-1)*delx	
enddo
do j=ny1,ny2
   ycord(j)=(j-1)*dely	
enddo
do k=nz1,nz2
   zcord(k)=(k-1)*delz	
enddo



end subroutine