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











subroutine compute_spectrum(pin,p,work,spectrum,spectrum_x,spectrum_y,spectrum_z,iwave_max,pe)
!
!  INPUT:  iwave_max:  size of spectrum()
!  OUTPUT: iwave_max:  number of coefficients returned in spectrum()
!          spectrum()  spherical wave number spectrum
!          spectrum_x  spectrum in x
!          spectrum_y  spectrum in y
!          spectrum_z  spectrum in z
!
!
use params
use mpi
implicit none
integer :: iwave_max,ierr
integer :: pe             ! compute spectrum on this processor
real*8 :: pin(nx,ny,nz)
real*8 :: work(nx,ny,nz)
real*8 :: p(nx,ny,nz)
real*8 :: spectrum(0:iwave_max)
real*8 :: spectrum_x(0:g_nx/2)
real*8 :: spectrum_y(0:g_ny/2)
real*8 :: spectrum_z(0:g_nz/2)

! local variables
real*8 rwave
integer :: iwave
real*8 :: spectrum_in(0:iwave_max)
real*8 :: energy
integer i,j,k


rwave=sqrt(  (g_nx/2.0)**2 + (g_ny/2.0)**2 + (g_nz/2.0)**2 )
if (nint(rwave)>iwave_max) then
   call abort("compute_spectrum: called with insufficient storege for spectrum()")
endif
iwave_max=nint(rwave)


p=pin
call fft3d(p,work)
spectrum=0
spectrum_x=0
spectrum_y=0
spectrum_z=0

do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
    rwave = imcord(i)**2 + jmcord(j)**2 + kmcord(k)**2
    iwave = nint(sqrt(rwave))

    energy = 8
    if (kmcord(k)==0) energy=energy/2
    if (jmcord(j)==0) energy=energy/2
    if (imcord(i)==0) energy=energy/2
    energy=energy*p(i,j,k)*p(i,j,k)

    spectrum(iwave)=spectrum(iwave)+energy
    spectrum_x(abs(imcord(i)))=spectrum_x(abs(imcord(i))) + energy
    spectrum_y(abs(jmcord(j)))=spectrum_y(abs(jmcord(j))) + energy
    spectrum_z(abs(kmcord(k)))=spectrum_z(abs(kmcord(k))) + energy

enddo
enddo
enddo


#ifdef USE_MPI
spectrum_in=spectrum
call MPI_reduce(spectrum_in,spectrum,1+iwave_max,MPI_REAL8,MPI_SUM,pe,comm_3d,ierr)
spectrum_in(0:g_nx/2)=spectrum_x
call MPI_reduce(spectrum_in,spectrum_x,1+(g_nx/2),MPI_REAL8,MPI_SUM,pe,comm_3d,ierr)
spectrum_in(0:g_ny/2)=spectrum_y
call MPI_reduce(spectrum_in,spectrum_y,1+(g_ny/2),MPI_REAL8,MPI_SUM,pe,comm_3d,ierr)
spectrum_in(0:g_nz/2)=spectrum_z
call MPI_reduce(spectrum_in,spectrum_z,1+(g_nz/2),MPI_REAL8,MPI_SUM,pe,comm_3d,ierr)
#endif

if (g_nz == 1)  then
   iwave = min(g_nx,g_ny)
else
   iwave = min(g_nx,g_ny,g_nz)
endif
iwave = (iwave/2)           ! max wave number in sphere.

! for all waves outside sphere, sum into one wave number:
do i=iwave+2,iwave_max
   spectrum(iwave+1)=spectrum(iwave+1)+spectrum(i)
enddo
iwave_max=iwave+1



end subroutine

