#include "macros.h"







subroutine compute_div(Q,divmx,divi)
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

#ifdef USE_MPI
   g_mx=divmx
   call MPI_allreduce(g_mx,divmx,1,MPI_REAL8,MPI_MAX,comm_3d,ierr)
   g_mx=divi
   call MPI_allreduce(g_mx,divi,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
end subroutine











subroutine compute_spectrum(pin,spectrum,iwave_max,pe)
use params
use mpi
implicit none
integer :: iwave_max,ierr
integer :: pe             ! compute spectrum on this processor
real*8 :: pin(nx,ny,nz)
real*8 :: spectrum(0:iwave_max)

! local variables
real*8 :: p(nx,ny,nz)
real*8 rwave
integer :: iwave
real*8 :: spectrum_in(0:iwave_max)
real*8 :: work(nx,ny,nz)
integer i,j,k

rwave=sqrt(  (g_nx/2.0)**2 + (g_ny/2.0)**2 + (g_nz/2.0)**2 )
if (nint(rwave)>iwave_max) then
   call abort("compute_spectrum: called with insufficient storege for spectrum()")
endif
iwave_max=nint(rwave)

p=pin
call fft3d(p,work)
spectrum=0

do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
    rwave = imcord(i)**2 + jmcord(j)**2 + kmcord(k)**2
    iwave = nint(sqrt(rwave))
    ASSERT("compute spectrum: iwave>iwave_max",iwave<=iwave_max)
    spectrum(iwave)=spectrum(iwave)+p(i,j,k)**2    
enddo
enddo
enddo

#ifdef USE_MPI
spectrum_in=spectrum
call MPI_reduce(spectrum_in,spectrum,iwave_max,MPI_REAL8,MPI_SUM,pe,comm_3d,ierr)
#endif

iwave = max(g_nx,g_ny,g_nz)
iwave = (iwave/2)           ! max wave number in sphere.

! for all waves outside sphere, sum into one wave number:
do i=iwave+2,iwave_max
   spectrum(iwave+1)=spectrum(iwave+1)+spectrum(i)
enddo
iwave_max=iwave+1



end subroutine