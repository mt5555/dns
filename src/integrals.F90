#include "macros.h"

subroutine compute_ke(u,pe,ke)
use params
implicit none
integer :: pe             ! compute totals on this processor
real*8 :: u(nx,ny,nz,3)
real*8 :: ke

! local variables
integer i,j,k,n
real*8 ke2
                                     ! vor = ints(2)

ke=0
do n=1,3
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   ke=ke+u(i,j,k,n)**2
enddo
enddo
enddo
enddo
ke=ke/2


#ifdef MPI
ke2=ke
call MPI_reduce(ke2,ke,1,MPI_REAL8,MPI_SUM,pe,comm_3d,ierr)
#endif

end subroutine



#if 0
subroutine compute_ke_diff(Q,delt,pe,ke,ke_diff,ke_diff_vis)
use params
implicit none
integer pe
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: divmx,divi

!local variablesn
real*8 :: Q_tmp(nx,ny,nz,n_var)
real*8 :: work(nx,ny,nz)
integer i,j,k

Q_tmp=Q
diff=0
call ns3D(rhs,Q_old,time,diff)
call divfree(rhs,Q_tmp(1,1,1,1),Q_tmp(1,1,1,2),Q_tmp(1,1,1,3))
Q_tmp=Q+delt*rhs

! estimate of total KE disapation:
! d/dt(.5*Q^2) 
call compute_ke(Q_tmp,io_pe,tmp)
call compute_ke(Q_old,io_pe,ke_diss_tot)
ke_diss_tot = (tmp-ke_diss_tot)/delt

Q_old=.5*(Q_tmp+Q_old)
diff=Q_old*diff
ke_diss_diff=sum(diff(nx1:nx2,ny1:ny2,nz1:nz2,1:3))

print *,'delta ke: ',ke_diss_tot,ke_diss_diff
end subroutine
#endif








subroutine compute_div(Q,pe,divmx,divi)
use params
implicit none
integer pe
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: divmx,divi

!local variablesn
real*8 :: p(nx,ny,nz)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
integer i,j,k

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

#ifdef MPI
   call MPI_REDUCE(mx,MPI_MAX)
   g_mx=divmx
   call MPI_reduce(g_mx,mx,1,MPI_REAL8,MPI_MAX,pe,comm_3d,ierr)
   g_mx=divi
   call MPI_reduce(g_mx,divi,1,MPI_REAL8,MPI_SUM,pe,comm_3d,ierr)
#endif
end subroutine











subroutine compute_spectrum(p,spectrum,iwave_max,pe)
use params
implicit none
integer :: iwave_max
integer :: pe             ! compute spectrum on this processor
real*8 :: p(nx,ny,nz)
real*8 :: spectrum(0:iwave_max)

! local variables
real*8 rwave,iwave
real*8 :: spectrum_in(0:iwave_max)
real*8 :: work(nx,ny,nz)
integer i,j,k

rwave=sqrt(  (g_nx/2.0)**2 + (g_ny/2.0)**2 + (g_nz/2.0)**2 )
if (nint(rwave)>iwave_max) then
   call abort("compute_spectrum: called with insufficient storege for spectrum()")
endif
iwave_max=nint(rwave)


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

#ifdef MPI
spectrum_in=spectrum
call MPI_reduce(spectrum_in,spectrum,iwave_max,MPI_REAL8,MPI_SUM,pe,comm_3d,ierr)
#endif

end subroutine