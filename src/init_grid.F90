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
integer i,j,k,l


delx = one/g_nx
dely = one/g_ny
delz = one/g_nz

do i=1,g_nx
   g_xcord(i)=(i-1)*delx	
enddo
do j=1,g_ny
   g_ycord(j)=(j-1)*dely	
enddo
do k=1,g_nz
   g_zcord(k)=(k-1)*delz	
enddo

call fft_get_mcord(g_imcord,g_nx)
call fft_get_mcord(g_jmcord,g_ny)
call fft_get_mcord(g_kmcord,g_nz)

do i=nx1,nx2
   l = i-nx1+1 + nslabx*mpicoords(1)
   xcord(i)=g_xcord(l)
   imcord(i)=g_imcord(l)
enddo
do j=ny1,ny2
   l = j-ny1+1 + nslaby*mpicoords(2)
   ycord(j)=g_ycord(l)
   jmcord(j)=g_jmcord(l)
enddo
do k=nz1,nz2
   l = k-nz1+1 + nslabz*mpicoords(3)
   zcord(k)=g_zcord(l)
   kmcord(k)=g_kmcord(l)
enddo






end subroutine
