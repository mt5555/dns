#include "macros.h"
module spectrum
use params
implicit none

#if 0

module to compute spherical and other 1D spectrums and
transfer function spectrums


Example on how to compute spectrum of u dot forcing term:

if (compute_integrals .and. compute_transfer) then
   call compute_spectrum_z_fft(Qhat,io_pe,spec_velocity,iwave_max)
   call compute_spectrum_z_fft(rhs,io_pe,spec_term,iwave_max)
   call transfer_spec_components(time,spec_f,spec_velocity,spec_term,0)
endif

< accumulate forcing into rhs... >

if (compute_integrals .and. compute_transfer) then
   call compute_spectrum_z_fft(rhs,io_pe,spec_term,iwave_max)
   call transfer_spec_components(time,spec_f,spec_velocity,spec_term,1)
endif
   


#endif

logical :: compute_transfer=.false.

real*8,private ::  spec_x(0:g_nx/2)
real*8,private ::  spec_y(0:g_ny/2)
real*8,private ::  spec_z(0:g_nz/2)
real*8,private ::  spec_r(0:max(g_nx,g_ny,g_nz))
real*8,private ::  spec_r_old(0:max(g_nx,g_ny,g_nz))
real*8,private ::  transfer_r(0:max(g_nx,g_ny,g_nz))
real*8,private ::  time_old=-1

integer,private :: iwave=-1



real*8,private :: transfer_comp_time         ! time at which below terms evaluated at:
real*8 ::  spec_diff(0:max(g_nx,g_ny,g_nz))  ! u dot diffusion term
real*8 ::  spec_f(0:max(g_nx,g_ny,g_nz))     ! u dot forcing term
real*8 ::  spec_rhs(0:max(g_nx,g_ny,g_nz))   ! transfer spec of u dot RHS
real*8 ::  spec_velocity(0:max(g_nx,g_ny,g_nz),3)     ! provided for calling program to use
real*8 ::  spec_term(0:max(g_nx,g_ny,g_nz),3)     ! provided for calling program to use

contains






subroutine transfer_spec_components(time,spec_tot,spec_u,spec_t,accum)
!
!  
!  spec_u:  spectrum of 3 components of velocity
!  spec_t:  spectrum of 3 components of a term in the equations
!
! OUTPUT:
!  spec_tot:  spec_u dot spec_t             if accum=0
!  spec_tot:  spec_u dot spec_t - spec_tot    if accum=1
! 
!  usually, use accum=0.  But to compute spectrum of a function which is
!  accumulated directly into the rhs, call this routine with rhs spectrum before 
!  (with accum=0) and then call with rhs after (with accum=1)
! 
use params
implicit none
integer :: operation,accum,decomp
real*8 ::  spec_tot(0:max(g_nx,g_ny,g_nz))
real*8 ::  spec_u(0:max(g_nx,g_ny,g_nz),3)
real*8 ::  spec_t(0:max(g_nx,g_ny,g_nz),3)
real*8 :: time

! local
real*8 ::  spec_new(0:max(g_nx,g_ny,g_nz))
integer i,n

if (accum==0)  then
   spec_tot=0
else
   spec_tot=-spec_tot
endif

do n=1,ndim
do i=1,iwave
   spec_tot(i)=spec_tot(i) + spec_t(i,n)*spec_u(i,n)
enddo
enddo

transfer_comp_time=time

end subroutine





subroutine compute_spec(time,Q,q1,work1,work2)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time

!local
integer :: iwave_max,i
real*8 spec_x2(0:g_nx/2)
real*8 spec_y2(0:g_ny/2)
real*8 spec_z2(0:g_nz/2)
real*8 ::  spec_r2(0:max(g_nx,g_ny,g_nz))

iwave_max=max(g_nx,g_ny,g_nz)
spec_r=0
spec_r2=0
spec_x=0
spec_y=0
spec_z=0

q1=Q


do i=1,ndim
   call compute_spectrum(q1(1,1,1,i),work1,work2,spec_r2,spec_x2,spec_y2,spec_z2,iwave_max,io_pe)
   spec_r=spec_r+.5*spec_r2
   spec_x=spec_x + .5*spec_x2
   spec_y=spec_y + .5*spec_y2
   spec_z=spec_z + .5*spec_z2
enddo

spec_r_old=spec_r
time_old=time

end subroutine






subroutine compute_tran(time,Q,q1,work1,work2)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time

!local
integer :: iwave_max,i
real*8 ::  spec_r2(0:max(g_nx,g_ny,g_nz))
real*8 :: diss_tot
character(len=80) :: message

iwave_max=max(g_nx,g_ny,g_nz)
spec_r=0
spec_r2=0
spec_x=0
spec_y=0
spec_z=0

q1=Q

! compute spectrum in spec_r
do i=1,ndim
   call fft3d(q1(1,1,1,i),work1)
   call compute_spectrum_fft(q1(1,1,1,i),io_pe,spec_r2,iwave_max)
   spec_r=spec_r+.5*spec_r2
enddo

! compute time rate of change in transfer_r()
if (time-time_old > 0) then

   diss_tot=0
   do i=0,iwave
      transfer_r(i)=(spec_r(i)-spec_r_old(i) ) / (time-time_old)
      diss_tot = diss_tot + transfer_r(i)
   enddo
   
   print *,'diss tot: ',diss_tot
   print *,'min: ',minval(transfer_r(0:iwave)) 
   print *,'max: ',maxval(transfer_r(0:iwave)) 
   
endif

end subroutine







subroutine output_spec(time,Q,q1,q2,q3,work1,work2)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: q2(nx,ny,nz,n_var)
real*8 :: q3(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time

! local variables
integer i,j,k,n
integer :: ierr
character,save :: access="0"

real*8 :: x,divx,divi
character(len=80) :: message
CPOINTER fid




! append to output files, unless this is first call.
if (access=="0") then
   access="w"
else
   access="a"
endif



write(message,'(a,f10.4)') " KE spectrum",time
call logplotascii(spec_r,iwave,message(1:25))
!call logplotascii(spec_x,g_nx/2,message)
!call logplotascii(spec_y,g_ny/2,message)
!call logplotascii(spec_z,g_nz/2,message)



! for incompressible equations, print divergence as diagnostic:
if (equations==NS_UVW) then
   call compute_div(Q,q1,work1,work2,divx,divi)
   write(message,'(3(a,e12.5))') 'max(div)=',divx
   call print_message(message)	
endif



if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time_initial
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".spec"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "spec_write(): Error opening file errno=",ierr
      call abort(message)
   endif
   call cwrite8(fid,time,1)
   x=1+iwave; call cwrite8(fid,x,1)
   call cwrite8(fid,spec_r,1+iwave)
   x=1+g_nx/2; call cwrite8(fid,x,1)
   call cwrite8(fid,spec_x,1+g_nx/2)
   x=1+g_ny/2; call cwrite8(fid,x,1)
   call cwrite8(fid,spec_y,1+g_ny/2)
   x=1+g_nz/2; call cwrite8(fid,x,1)
   call cwrite8(fid,spec_z,1+g_nz/2)
   call cclose(fid,ierr)
endif
end subroutine







subroutine output_tran(time,Q,q1,q2,q3,work1,work2)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: q2(nx,ny,nz,n_var)
real*8 :: q3(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time

! local variables
integer i,j,k,n
integer :: ierr
character,save :: access="0"

real*8 :: x,ymax,ymax2
character(len=80) :: message
CPOINTER fid

if (time-time_old <= 0) return


! append to output files, unless this is first call.
if (access=="0") then
   access="w"
else
   access="a"
endif


ymax=maxval(transfer_r(0:iwave))
ymax=max(ymax,-minval(transfer_r(0:iwave)))

ymax2=.001
if (ymax>ymax2) ymax2=.0025
if (ymax>ymax2) ymax2=.0050
if (ymax>ymax2) ymax2=.0100
if (ymax>ymax2) ymax2=.0250
if (ymax>ymax2) ymax2=.0500
if (ymax>ymax2) ymax2=.1000
if (ymax>ymax2) ymax2=.2500
if (ymax>ymax2) ymax2=.5000
if (ymax>ymax2) ymax2=1


write(message,'(a,f10.4)') " KE transfer spectrum",time
call plotascii(transfer_r,iwave,message(1:25),-ymax2,ymax2)



if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time_initial
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".spect"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "transfer spec_write(): Error opening file errno=",ierr
      call abort(message)
   endif
   call cwrite8(fid,.5*(time+time_old),1)
   x=1+iwave; call cwrite8(fid,x,1)
   call cwrite8(fid,transfer_r,1+iwave)
   call cclose(fid,ierr)
endif
end subroutine









subroutine compute_spectrum(pin,p,work,spectrum,spectrum_x,spectrum_y,spectrum_z,iwave_max,pe)
!
!  INPUT:  iwave_max:  size of spectrum()
!  OUTPUT: iwave:      number of coefficients returned in spectrum()
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
iwave=iwave+1



end subroutine








subroutine compute_spectrum_fft(p,pe,spec_r,iwave_max)
use params
use mpi
implicit none
integer :: iwave_max,ierr
integer :: pe             ! compute spectrum on this processor
real*8 :: p(nx,ny,nz)
real*8 :: spec_r(0:iwave_max)

! local variables
real*8 rwave
real*8 :: spec_r_in(0:iwave_max)
real*8 :: energy
integer i,j,k


rwave=sqrt(  (g_nx/2.0)**2 + (g_ny/2.0)**2 + (g_nz/2.0)**2 )
if (nint(rwave)>iwave_max) then
   call abort("compute_spectrum: called with insufficient storege for spectrum()")
endif
iwave_max=nint(rwave)

spec_r=0

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

    spec_r(iwave)=spec_r(iwave)+energy

enddo
enddo
enddo


#ifdef USE_MPI
spec_r_in=spec_r
call MPI_reduce(spec_r_in,spec_r,1+iwave_max,MPI_REAL8,MPI_SUM,pe,comm_3d,ierr)
#endif

if (g_nz == 1)  then
   iwave = min(g_nx,g_ny)
else
   iwave = min(g_nx,g_ny,g_nz)
endif
iwave = (iwave/2)           ! max wave number in sphere.

! for all waves outside sphere, sum into one wave number:
do i=iwave+2,iwave_max
   spec_r(iwave+1)=spec_r(iwave+1)+spec_r(i)
enddo
iwave=iwave+1

end subroutine






subroutine compute_spectrum_z_fft(p,pe,spec_r,iwave_max)
use params
use mpi
implicit none
integer :: iwave_max,ierr
integer :: pe             ! compute spectrum on this processor
real*8 :: p(g_nz2,nslabx,ny_2dz)
real*8 :: spec_r(0:iwave_max)

! local variables
real*8 rwave
real*8 :: spec_r_in(0:iwave_max)
real*8 :: energy
integer i,j,k,jm,km,im


rwave=sqrt(  (g_nx/2.0)**2 + (g_ny/2.0)**2 + (g_nz/2.0)**2 )
if (nint(rwave)>iwave_max) then
   call abort("compute_spectrum: called with insufficient storege for spectrum()")
endif
iwave_max=nint(rwave)

spec_r=0

do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nslabx
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

         rwave = im**2 + jm**2 + km**2
         iwave = nint(sqrt(rwave))
         
         energy = 8
         if (km==0) energy=energy/2
         if (jm==0) energy=energy/2
         if (im==0) energy=energy/2
         energy=energy*p(i,j,k)*p(i,j,k)
         
         spec_r(iwave)=spec_r(iwave)+energy
         
enddo
enddo
enddo


#ifdef USE_MPI
spec_r_in=spec_r
call MPI_reduce(spec_r_in,spec_r,1+iwave_max,MPI_REAL8,MPI_SUM,pe,comm_3d,ierr)
#endif

if (g_nz == 1)  then
   iwave = min(g_nx,g_ny)
else
   iwave = min(g_nx,g_ny,g_nz)
endif
iwave = (iwave/2)           ! max wave number in sphere.

! for all waves outside sphere, sum into one wave number:
do i=iwave+2,iwave_max
   spec_r(iwave+1)=spec_r(iwave+1)+spec_r(i)
enddo
iwave=iwave+1

end subroutine



end module
