#include "macros.h"
module spectrum
use params
use sforcing
implicit none

#if 0

module to compute spherical and other 1D spectrums and
transfer function spectrums

Spectrum routines know about Lz for setting aspect ration.
Spherical shells are stored in integer wave number k,
and are of thickness of thickness 2*pi/Lz

Max wave number:   2*pi*(g_nz/2)
Number of shells:  .5*g_nz*Lz
However, we only allocate space for .5*g_nz shells, so this code will
print an error message if Lz>1



#endif

logical :: compute_transfer=.false.

real*8,private ::  spec_x(0:g_nx/2,n_var)
real*8,private ::  spec_y(0:g_ny/2,n_var)
real*8,private ::  spec_z(0:g_nz/2,n_var)
real*8,private ::  spec_r(0:max(g_nx,g_ny,g_nz),n_var)
real*8,private ::  spec_r_new(0:max(g_nx,g_ny,g_nz))
real*8,private ::  edot_r(0:max(g_nx,g_ny,g_nz))
real*8,private ::  spec_r_2d(0:max(g_nx,g_ny),0:g_nz/2,n_var)
real*8,private ::  q2spec_r_2d(0:max(g_nx,g_ny),0:g_nz/2)
real*8,private ::  time_old=-1

!helicity related spectra (sk)
real*8,private ::  spec_helicity_rp(0:max(g_nx,g_ny,g_nz))
real*8,private ::  spec_helicity_rn(0:max(g_nx,g_ny,g_nz))
real*8,private ::  spec_varhe(0:max(g_nx,g_ny,g_nz))
real*8,private ::  spec_meanhe(0:max(g_nx,g_ny,g_nz))
real*8 ::  spec_E(0:max(g_nx,g_ny,g_nz)) !E(k) from helicity free modes 
real*8 ::  spec_kEk(0:max(g_nx,g_ny,g_nz))  ! k E(k)
real*8 ::  spec_meanE(0:max(g_nx,g_ny,g_nz)) ! mean energy in each mode
real*8 ::  spec_varE(0:max(g_nx,g_ny,g_nz)) ! variance of energy in each mode
real*8 ::  E33(0:max(g_nx,g_ny,g_nz))  ! E_33 component of spec tensor
real*8 ::  II2sq(0:max(g_nx,g_ny,g_nz))  ! I_2^2 square of II component along RR
real*8 ::  RRsq(0:max(g_nx,g_ny,g_nz))  ! RR^2 square of real part only
real*8 ::  I2I3(0:max(g_nx,g_ny,g_nz))  ! I_2*I_3
real*8 ::  R2I3(0:max(g_nx,g_ny,g_nz))  ! mod_rr*I_3 
real*8 ::  par(0:max(g_nx,g_ny,g_nz))   ! abs(R\cross I)/(R^2+I^2)
real*8 ::  cos_tta_spec_p(0:max(g_nx,g_ny,g_nz)) !pos spec of cos_tta betn RR and II
real*8 ::  cos_tta_spec_n(0:max(g_nx,g_ny,g_nz)) !neg spec of cos_tta
real*8 ::  rhel_spec_p(0:max(g_nx,g_ny,g_nz)) !pos spec of rel hel
real*8 ::  rhel_spec_n(0:max(g_nx,g_ny,g_nz)) !neg spec of rel hel
real*8 ::  rhel_rms_spec(0:max(g_nx,g_ny,g_nz)) ! spec of rms of rel hel
real*8 ::  costta_pdf(0:max(g_nx,g_ny,g_nz),200) !pdfs of cos(tta) for each k
real*8 ::  tmp_pdf(200)                          !pdfs of cos(tta) for each k
real*8 ::  rhel_pdf(0:max(g_nx,g_ny,g_nz),200)  !pdfs of rel hel for each k

! cospectrum in x,y,z directions.
! n_var=1,3:  uv, uw, vw
real*8,private ::  cospec_r(0:max(g_nx,g_ny,g_nz),n_var)
real*8,private ::  cospec_x(0:g_nx/2,n_var)   
real*8,private ::  cospec_y(0:g_ny/2,n_var)
real*8,private ::  cospec_z(0:g_nz/2,n_var)
real*8,private ::  q2spec_x(0:g_nx/2)
real*8,private ::  q2spec_y(0:g_ny/2)
real*8,private ::  q2spec_z(0:g_nz/2)
real*8,private ::  q2spec_r(0:max(g_nx,g_ny,g_nz))
real*8,private ::  bspec_x(0:g_nx/2,n_var)
real*8,private ::  bspec_y(0:g_ny/2,n_var)
real*8,private ::  bspec_z(0:g_nz/2,n_var)
real*8,private ::  bspec_r(0:max(g_nx,g_ny,g_nz))


integer,private :: iwave=-1
integer,private :: iwave_2d=-1
integer,private :: nbin=200


real*8 ::  transfer_comp_time         ! time at which below terms evaluated at:

real*8 ::  spec_diff(0:max(g_nx,g_ny,g_nz))  ! u dot diffusion term
real*8 ::  spec_diff_new(0:max(g_nx,g_ny,g_nz)) 
real*8 ::  spec_f(0:max(g_nx,g_ny,g_nz))     ! u dot forcing term
real*8 ::  spec_model(0:max(g_nx,g_ny,g_nz)) ! spectrum of div(tau) or smagorinsky term
real*8 ::  spec_model_new(0:max(g_nx,g_ny,g_nz)) 
real*8 ::  spec_rhs(0:max(g_nx,g_ny,g_nz))   ! transfer spec of u dot RHS
real*8 ::  spec_ens(0:max(g_nx,g_ny,g_nz))   ! enstropy spectrum
real*8 ::  spec_ens2(0:max(g_nx,g_ny,g_nz))   ! enstropy spectrum

real*8 ::  spec_tmp(0:max(g_nx,g_ny,g_nz))   ! storage, for calling program convienience

contains



subroutine compute_spec(time,Q,q1,work1,work2)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)   ! (u,v,w)
real*8 :: q1(nx,ny,nz,n_var) 
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time

!local
integer :: iwave_max,i,n,j
real*8 ::  spec_r2(0:max(g_nx,g_ny,g_nz))
real*8 ::  spec_d2(0:max(g_nx,g_ny,g_nz))

iwave_max=max(g_nx,g_ny,g_nz)
spec_r=0
spec_diff=0
spec_r2=0
spec_d2=0
spec_x=0
spec_y=0
spec_z=0

! compute FFT in q1  of velocity:
do i=1,ndim
   q1(:,:,:,i)=Q(:,:,:,i)
   call fft3d(q1(1,1,1,i),work1)
enddo
! compute FFT of scalars:
do i=np1,np2
   q1(:,:,:,i)=Q(:,:,:,i)
   call fft3d(q1(1,1,1,i),work1)
enddo
!
!   pe=pe+.5*grav*Q(i,j,3)**2 - .5*grav*H0**2
!   ke = ke + .5*Q(i,j,3)*Q(i,j,n)**2   n=1,2
!
!
! passive scalars:
do n=np1,np2
   call compute_spectrum(q1(1,1,1,n),work1,work2,spec_r(0,n),spec_r2,&
       spec_x(0,n),spec_y(0,n),spec_z(0,n),iwave_max,1)
   spec_r(:,n)=.5*spec_r(:,n)
enddo
do i=1,ndim
   call compute_spectrum(q1(1,1,1,i),work1,work2,spec_r2,spec_d2,&
       spec_x(0,i),spec_y(0,i),spec_z(0,i),iwave_max,1)
   spec_r(:,1)=spec_r(:,1)+.5*spec_r2
   ! for now, use the value computed in RHS
   ! spec_diff=spec_diff + spec_d2
enddo
spec_x=.5*spec_x
spec_y=.5*spec_y
spec_z=.5*spec_z


time_old=time

#undef SPEC_DIFF_CHECK
#ifdef SPEC_DIFF_CHECK
spec_diff_new=spec_diff  ! make a copy of spec_diff for check below
#endif

! q1 already contains the FFT of Q, so set skip_fft (last arg) to 1:
if (ndim>=3) then
   call compute_helicity_spectrum(Q,q1,work1,1)
endif
end subroutine




subroutine compute_pv2_spec(time,Q,q1,q2,q3,work1,work2)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)   ! (u,v,w)
real*8 :: q1(nx,ny,nz,n_var) 
real*8 :: q2(nx,ny,nz,n_var) 
real*8 :: q3(nx,ny,nz,n_var) 
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time, pv2energy
real*8 pv(nx,ny,nz)  ! used for potential vorticity
real*8 pv2(nx,ny,nz)  ! used for potential enstrophy

!local
integer :: iwave_max,n,pv_type, i,j,k,iw
real*8 ::  spec_r2(0:max(g_nx,g_ny,g_nz))
real*8 :: rwave,xfac
!
! use the full pv
!
pv_type = 1
iwave_max=max(g_nx,g_ny,g_nz)
!bw
!bw For potential enstrophy Q = q^2/2
!bw
q2spec_r=0
q2spec_x=0
q2spec_y=0
q2spec_z=0
q2spec_r_2d=0

! 
!   pv =  (vorticity + f) dot grad rho
!   Q = pv^2/2 
!
!bw   Compute the vorticity and the potential vorticity
!bw
!bw
!bw I probably don't have the work arrays correct in these calls.
!bw
! compute pv in work1, vorticity in q1
call potential_vorticity(work1,q1,Q,q2,q3,pv_type)

call fft3d(work1,q1)

call compute_spectrum(work1,q2,q3,q2spec_r,spec_r2,&
     q2spec_x,q2spec_y,q2spec_z,iwave_max,1)

!bw computing pv is complicated so set all those to zero for the moment
q2spec_r=.5*q2spec_r
q2spec_x=.5*q2spec_x
q2spec_y=.5*q2spec_y
q2spec_z=.5*q2spec_z


call compute_spectrum_2d(work1,q1,q2,q2spec_r_2d,iwave_max,1)


end subroutine




subroutine compute_spec_2d(time,Q,q1,work1,work2)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)   ! (u,v,w)
real*8 :: q1(nx,ny,nz,n_var) 
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time

!local
integer :: iwave_max,i

iwave_max=max(g_nx,g_ny)
spec_r_2d=0


do i=1,n_var
   call compute_spectrum_2d(Q(1,1,1,i),work1,work2,spec_r_2d(0,0,i),iwave_max,0)
enddo

spec_r_2d=spec_r_2d/2



end subroutine






subroutine compute_spec_shallow(time,Q,q1,work1,work2)
use params
implicit none
real*8 :: Q(nx,ny,nz,3)   ! (u,v,h)
real*8 :: q1(nx,ny,nz,3)  ! (u,v,h)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time

!local
integer :: iwave_max,i
real*8 ::  spec_r2(0:max(g_nx,g_ny,g_nz))
real*8 ::  spec_d2(0:max(g_nx,g_ny,g_nz))

iwave_max=max(g_nx,g_ny,g_nz)
spec_r=0
spec_r2=0
spec_d2=0
spec_x=0
spec_y=0
spec_z=0


q1(:,:,:,1)=sqrt(Q(:,:,:,3))*Q(:,:,:,1)
q1(:,:,:,2)=sqrt(Q(:,:,:,3))*Q(:,:,:,2)
q1(:,:,:,3)=(Q(:,:,:,3))

!
!   pe=pe+.5*grav*Q(i,j,3)**2 - .5*grav*H0**2
!   ke = ke + .5*Q(i,j,3)*Q(i,j,n)**2   n=1,2
!
!   ke = .5 * spec(q1(:,:,:,1)) + .5*spec(q1(:,:,:,2))  
!   pe = .5*grav* ( spec(q1(:,:,:,3)) - H0**2(0,0,0 mode only)
!
!
!
!  PE part
!
call compute_spectrum(q1(1,1,1,3),work1,work2,spec_r2,spec_d2,&
       spec_x(0,1),spec_y(0,1),spec_z(0,1),iwave_max,0)
spec_r(:,1)=spec_r(:,1)+.5*grav*spec_r2
spec_r(0,1)=spec_r(0,1) - .5*grav*H0**2

!
!  KE part
!                 spec_[x,y,z] = KE only
!
do i=1,ndim
   call compute_spectrum(q1(1,1,1,i),work1,work2,spec_r2,spec_d2,&
       spec_x(0,i),spec_y(0,i),spec_z(0,i),iwave_max,0)
   spec_r(:,1)=spec_r(:,1)+.5*spec_r2
enddo
spec_x=.5*spec_x
spec_y=.5*spec_y
spec_z=.5*spec_z



time_old=time


end subroutine





subroutine compute_Edotspec_shallow(time,Q,q1,work1,work2)
use params
implicit none
real*8 :: Q(nx,ny,nz,3)   ! (u,v,h)
real*8 :: q1(nx,ny,nz,3)  ! work array
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time

!local
integer :: iwave_max,i
real*8 ::  spec_r2(0:max(g_nx,g_ny,g_nz))
real*8 ::  spec_d2(0:max(g_nx,g_ny,g_nz))

iwave_max=max(g_nx,g_ny,g_nz)
spec_r_new=0




spec_r2=0
spec_d2=0
spec_x=0
spec_y=0
spec_z=0


q1(:,:,:,1)=sqrt(Q(:,:,:,3))*Q(:,:,:,1)
q1(:,:,:,2)=sqrt(Q(:,:,:,3))*Q(:,:,:,2)
q1(:,:,:,3)=(Q(:,:,:,3))

!
!   pe=pe+.5*grav*Q(i,j,3)**2 - .5*grav*H0**2
!   ke = ke + .5*Q(i,j,3)*Q(i,j,n)**2   n=1,2
!
!   ke = .5 * spec(q1(:,:,:,1)) + .5*spec(q1(:,:,:,2))  
!   pe = .5*grav* ( spec(q1(:,:,:,3)) - H0**2(0,0,0 mode only)
!
!
!
!  PE part
!
call compute_spectrum(q1(1,1,1,3),work1,work2,spec_r2,spec_d2,&
       spec_x(0,1),spec_y(0,1),spec_z(0,1),iwave_max,0)
spec_r_new=spec_r_new+.5*grav*spec_r2
spec_r_new(0)=spec_r_new(0) - .5*grav*H0**2

!
!  KE part
!                 spec_[x,y,z] = KE only
!
do i=1,ndim
   call compute_spectrum(q1(1,1,1,i),work1,work2,spec_r2,spec_d2,&
       spec_x(0,i),spec_y(0,i),spec_z(0,i),iwave_max,0)
   spec_r_new=spec_r_new+.5*spec_r2
enddo

! compute time rate of change in edot_r()
if (time-time_old > 0) then
   do i=0,iwave
      edot_r(i)=(spec_r_new(i)-spec_r(i,1) ) / (time-time_old)
   enddo
endif



end subroutine






subroutine compute_Edotspec(time,Q,q1,work1,work2)
!
! This routine assumes compute_spec was called at the end
! of the last time step.  so we have this data at time_old:
! spec_r, spec_diff
! This routine will save this data, and compute new values
! of spec_r, spec_diff at the current time.  
!
! It will then compute the total KE rate of change spectrum, 
! edot_r
!
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)   ! 
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time

!local
integer :: iwave_max,i
real*8 ::  spec_r2(0:max(g_nx,g_ny,g_nz))
real*8 ::  spec_d2(0:max(g_nx,g_ny,g_nz))
character(len=100) :: message


#ifdef SPEC_DIFF_CHECK
!spec_diff_new:  spec_diff computed at t=0 call to compute_spec().
!spec_diff:      computed in ns.F90 because flag was set
!check that they are the same:
print *,'**************************************'
print *,'maxval spec_diff check: ',maxval(abs(spec_diff-spec_diff_new))
print *,'norms: ',maxval(abs(spec_diff)),maxval(abs(spec_diff_new))
print *,'**************************************'
#endif





!
! save data computed with last call to compute_spec  (time_old)
!


iwave_max=max(g_nx,g_ny,g_nz)
spec_r_new=0
spec_diff_new=0

spec_r2=0
spec_d2=0
spec_x=0
spec_y=0
spec_z=0


q1=Q

! compute spectrum in spec_r
do i=1,ndim
   call compute_spectrum(q1(1,1,1,i),work1,work2,spec_r2,spec_d2,&
       spec_x(0,i),spec_y(0,i),spec_z(0,i),iwave_max,0)
   spec_r_new=spec_r_new+.5*spec_r2
   ! for now, use the value computed in RHS
   ! spec_diff_new=spec_diff_new + spec_d2  
enddo


! compute time rate of change in edot_r()
if (time-time_old > 0) then
   do i=0,iwave
      edot_r(i)=(spec_r_new(i)-spec_r(i,1) ) / (time-time_old)
   enddo
endif

end subroutine







subroutine output_spec(time,time_file)
use params
implicit none
real*8 :: time,time_file

! local variables
integer i,j,k,n
integer :: ierr
character,save :: access="0"

real*8 :: x
character(len=80) :: message
CPOINTER fid




! append to output files, unless this is first call.
if (access=="0" .or. time==time_file) then
   access="w"
else
   access="a"
endif

! If we are also computing 2d spectrum, dont bother to plot
! 3D spectrum 
if (iwave_2d <= iwave ) then
   write(message,'(a,f10.4)') " Energy Spectrum t=",time
   call logplotascii(spec_r(0,1),iwave,message(1:25))
   !call logplotascii(spec_x,g_nx/2,message)
   !call logplotascii(spec_y,g_ny/2,message)
   !call logplotascii(spec_z,g_nz/2,message)
endif





if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time_file
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".spec"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "spec_write(): Error opening file errno=",ierr
      call abortdns(message)
   endif
   call cwrite8(fid,time,1)
   x=1+iwave; call cwrite8(fid,x,1)
   call cwrite8(fid,spec_r(0,1),1+iwave)
   x=1+g_nx/2; call cwrite8(fid,x,1)
   do i=1,3
      call cwrite8(fid,spec_x(0,i),1+g_nx/2)
   enddo
   x=1+g_ny/2; call cwrite8(fid,x,1)
   do i=1,3
      call cwrite8(fid,spec_y(0,i),1+g_ny/2)
   enddo
   x=1+g_nz/2; call cwrite8(fid,x,1)
   do i=1,3
      call cwrite8(fid,spec_z(0,i),1+g_nz/2)
   enddo
   call cclose(fid,ierr)

   if (npassive>0) then
   write(message,'(f10.4)') 10000.0000 + time_file
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".pspec"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "spec_write(): Error opening file errno=",ierr
      call abortdns(message)
   endif
   x=npassive; call cwrite8(fid,x,1)
   call cwrite8(fid,time,1)

   x=1+iwave; call cwrite8(fid,x,1)
   do i=np1,np2
      call cwrite8(fid,spec_r(0,i),1+iwave) 
   enddo

   x=1+g_nx/2; call cwrite8(fid,x,1)
   do i=np1,np2
      call cwrite8(fid,spec_x(0,i),1+g_nx/2)
   enddo
   x=1+g_ny/2; call cwrite8(fid,x,1)
   do i=np1,np2
      call cwrite8(fid,spec_y(0,i),1+g_ny/2)
   enddo
   x=1+g_nz/2; call cwrite8(fid,x,1)
   do i=np1,np2	
      call cwrite8(fid,spec_z(0,i),1+g_nz/2)
   enddo
   call cclose(fid,ierr)
   endif
endif
end subroutine



subroutine output_helicity_spec(time,time_file)
use params
implicit none
real*8 :: time,time_file

! local variables
integer i,j,k,n
integer :: ierr
character,save :: access="0"

real*8 :: x
character(len=80) :: message
CPOINTER fid




! append to output files, unless this is first call.
if (access=="0" .or. time==time_file) then
   access="w"
else
   access="a"
endif


if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time_file
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".hspec"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "output_helicity_spec(): Error opening file errno=",ierr
      call abortdns(message)
   endif
   call cwrite8(fid,time,1)
   x=1+iwave; call cwrite8(fid,x,1)
   call cwrite8(fid,spec_helicity_rn,1+iwave)
   call cwrite8(fid,spec_helicity_rp,1+iwave)
   call cclose(fid,ierr)


   write(message,'(f10.4)') 10000.0000 + time_file
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".kspec"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "output_helicity_spec(): Error opening file errno=",ierr
      call abortdns(message)
   endif
   call cwrite8(fid,time,1)
   x=1+iwave; call cwrite8(fid,x,1)
   call cwrite8(fid,spec_kEk,1+iwave)
   call cclose(fid,ierr)


   write(message,'(f10.4)') 10000.0000 + time_file
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".cospec"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "output_cospec(): Error opening file errno=",ierr
      call abortdns(message)
   endif
   x=g_nx/2; call cwrite8(fid,x,1)
   x=g_ny/2; call cwrite8(fid,x,1)
   x=g_nz/2; call cwrite8(fid,x,1)
   x=1+iwave; call cwrite8(fid,x,1)
   call cwrite8(fid,time,1)
   do i=1,3
      call cwrite8(fid,cospec_x(0,i),g_nx/2+1)
   enddo
   do i=1,3
      call cwrite8(fid,cospec_y(0,i),g_ny/2+1)
   enddo
   do i=1,3
      call cwrite8(fid,cospec_z(0,i),g_nz/2+1)
   enddo
   do i=1,3
      call cwrite8(fid,cospec_r(0,i),1+iwave)
   enddo


   call cclose(fid,ierr)
endif
end subroutine




subroutine output_2d_spec(time,time_file)
use params
implicit none
real*8 :: time,time_file

! local variables
integer i,j,k,n
integer :: ierr
character,save :: access="0"

real*8 :: x
character(len=80) :: message
CPOINTER fid

if (iwave_2d<0) then
   ! no spectrum was computed
   call print_message("Warning: output_2d_spec() called, but no 2d spectra was computed")
   return
endif

write(message,'(a,f10.4)') " 2D Spectrum t=",time
call logplotascii(spec_r_2d(0,:,1),iwave_2d,message(1:25))


! append to output files, unless this is first call.
if (access=="0" .or. time==time_file) then
   access="w"
else
   access="a"
endif


if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time_file
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".spec2d"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "output_2d_spec(): Error opening file errno=",ierr
      call abortdns(message)
   endif
   call cwrite8(fid,time,1)
   x=n_var; call cwrite8(fid,x,1)  
   x=1+iwave_2d; call cwrite8(fid,x,1)  
   x=1+g_nz/2; call cwrite8(fid,x,1)
   do i=1,n_var
   do k=0,g_nz/2
      call cwrite8(fid,spec_r_2d(0,k,i),1+iwave_2d)
   enddo
   enddo
   call cclose(fid,ierr)

!   if (npassive>0) then
!      call abortdns("Error: ouput of 2d spectrum for scalers not coded")
!      ! add spectrum of scalars
!   endif

endif
end subroutine





subroutine output_hfree_spec(time,time_file)
use params
implicit none
real*8 :: time,time_file

! local variables
integer i,j,k,n
integer :: ierr
character,save :: access="0"

real*8 :: x
character(len=100) :: message
CPOINTER fid




! append to output files, unless this is first call.
if (access=="0" .or. time==time_file) then
   access="w"
else
   access="a"
endif

if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time_file
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname))// message(2:10) // ".hf_spec"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)')"output_hfree_spec(): Error opening file errno=",ierr
      call abortdns(message)
   endif
   call cwrite8(fid,time,1)
   x=1+iwave; call cwrite8(fid,x,1)   
   call cwrite8(fid,spec_helicity_rn,1+iwave)
   call cwrite8(fid,spec_helicity_rp,1+iwave)
   call cwrite8(fid,spec_meanhe,1+iwave)
   call cwrite8(fid,spec_varhe,1+iwave)
   call cwrite8(fid,spec_E,1+iwave)
   call cwrite8(fid,spec_kEk,1+iwave)
   call cwrite8(fid,spec_meanE,1+iwave)
   call cwrite8(fid,spec_varE,1+iwave)

   ! these quantities are disabled; see corresponding code in compute_hfree_spec
#if 0
   call cwrite8(fid,E33,1+iwave)
   call cwrite8(fid,II2sq,1+iwave)
   call cwrite8(fid,RRsq,1+iwave)
   call cwrite8(fid,I2I3,1+iwave)
   call cwrite8(fid,R2I3,1+iwave)
   call cwrite8(fid,par,1+iwave)
#endif

   call cwrite8(fid,rhel_spec_n,1+iwave)
   call cwrite8(fid,rhel_spec_p,1+iwave)
   call cwrite8(fid,rhel_rms_spec,1+iwave)
   call cwrite8(fid,cos_tta_spec_n,1+iwave)
   call cwrite8(fid,cos_tta_spec_p,1+iwave)

   x=nbin
   call cwrite8(fid,x,1)
   do i=1,1+iwave
      tmp_pdf=costta_pdf(i-1,:)
      call cwrite8(fid,tmp_pdf,nbin)
   enddo
   do i=1,1+iwave
      tmp_pdf=rhel_pdf(i-1,:)
      call cwrite8(fid,tmp_pdf,nbin)
   enddo   
   call cclose(fid,ierr)
   
endif
end subroutine





subroutine output_pv2_spec(time,time_file)
!bw
!bw  This subroutine outputs the potential enstrophy spectrum
!bw  The notation here is Q = q^2/2 = q2
!bw  
!bw  
!bw
use params
implicit none
real*8 :: time,time_file

! local variables
integer i,j,k,n
integer :: ierr
character,save :: access="0"

real*8 :: x
character(len=80) :: message
CPOINTER fid




! append to output files, unless this is first call.
if (access=="0" .or. time==time_file) then
   access="w"
else
   access="a"
endif

!
! Don't plot the potential enstrophy spectrum at all
! If we are also computing 2d spectrum, dont bother to plot
! 3D spectrum 
!if (iwave_2d <= iwave ) then
!   write(message,'(a,f10.4)') " Energy Spectrum t=",time
!   call logplotascii(spec_r(0,1),iwave,message(1:25))
!   !call logplotascii(spec_x,g_nx/2,message)
!   !call logplotascii(spec_y,g_ny/2,message)
!   !call logplotascii(spec_z,g_nz/2,message)
!endif

if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time_file
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".pv2spec"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "spec_write(): Error opening file errno=",ierr
      call abortdns(message)
   endif
   call cwrite8(fid,time,1)
   x=1+iwave; call cwrite8(fid,x,1)
   x=1+g_nx/2; call cwrite8(fid,x,1)
   x=1+g_ny/2; call cwrite8(fid,x,1)	 
   x=1+g_nz/2; call cwrite8(fid,x,1)	 
   x=1+iwave_2d; call cwrite8(fid,x,1)	
   call cwrite8(fid,q2spec_r,1+iwave)
   call cwrite8(fid,q2spec_x,1+g_nx/2)
   call cwrite8(fid,q2spec_y,1+g_ny/2)
   call cwrite8(fid,q2spec_z,1+g_nz/2)
   do k=0,g_nz/2
      call cwrite8(fid,q2spec_r_2d(0,k),1+iwave_2d)
   enddo
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
real*8 ::  spec_r2(0:max(g_nx,g_ny,g_nz))

real*8 :: x,ymax,ymax2
real*8 :: sum_tot, sum_tot2, sum_diss, sum_f

character(len=80) :: message
CPOINTER fid

if (time-time_old <= 0) return


! append to output files, unless this is first call.
if (access=="0") then
   access="w"
else
   access="a"
endif

if (transfer_comp_time/=time_old) then
   call print_message("WARNING: transfer spectrum not calculated")	
   return
endif

if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time_initial
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
      // message(2:10) // ".spect"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') &
      "transfer spec_write(): Error opening file errno=",ierr
      call abortdns(message)
      endif
      
      if (equations==NS_UVW .or. equations==CNS) then
         x=4                    ! number of spectrums in file for each time.  
         call cwrite8(fid,x,1)
         
                                ! d/dt total KE spectrum
         call cwrite8(fid,.5*(time+time_old),1) 
         x=1+iwave; call cwrite8(fid,x,1)
      call cwrite8(fid,edot_r,1+iwave)
      
      ! transfer function (only in NS case.  other cases spec_rhs is
      ! not computed and this is garbage.
      spec_r2 = spec_rhs - (spec_diff+spec_f)
      call cwrite8(fid,transfer_comp_time,1)  
      x=1+iwave; call cwrite8(fid,x,1)
      call cwrite8(fid,spec_r2,1+iwave)
      
                                ! diffusion spectrum
      call cwrite8(fid,transfer_comp_time,1) 
      x=1+iwave; call cwrite8(fid,x,1)
      call cwrite8(fid,spec_diff,1+iwave)
      
      ! forcing spectrum
      call cwrite8(fid,transfer_comp_time,1) 
      x=1+iwave; call cwrite8(fid,x,1)
      call cwrite8(fid,spec_f,1+iwave)
      endif
      if (equations==SHALLOW) then
      x=3   ! number of spectrums in file for each time.  
      call cwrite8(fid,x,1)
      
      ! e-dot
      call cwrite8(fid,.5*(time+time_old),1) 
      x=1+iwave; call cwrite8(fid,x,1)
      call cwrite8(fid,edot_r,1+iwave)

      ! diffusion spectrum
      call cwrite8(fid,.5*(time+time_old),1) 
      x=1+iwave; call cwrite8(fid,x,1)
      call cwrite8(fid,spec_diff,1+iwave)

      ! div(tau) (or smag) spectrum
      call cwrite8(fid,.5*(time+time_old),1) 
      x=1+iwave; call cwrite8(fid,x,1)
      call cwrite8(fid,spec_model,1+iwave)
   endif


   call cclose(fid,ierr)
#if 0
   write(*,'(a)') '    i      d(E_k)/dt       T      D'
   do i=0,min(iwave,60)
      write(*,'(i4,5f12.4)') i,edot_r(i),&
           edot_r(i)-spec_diff(i),spec_diff(i)
   enddo
#endif


#if 0
   sum_tot=0
   sum_tot2=0
   sum_diss=0
   sum_f=0
   write(*,'(a)') '    i      d(E_k)/dt      d(E_k)/dt      T      D      F'
   do i=0,min(iwave,20)
      write(*,'(i4,5f12.4)') i,edot_r(i),&
           spec_r2(i)+spec_diff(i)+spec_f(i), &
           spec_r2(i),spec_diff(i),spec_f(i)
      sum_tot=sum_tot+edot_r(i)
      sum_tot2=sum_tot2+spec_rhs(i)
      sum_f=sum_f + spec_f(i)
      sum_diss=sum_diss + spec_diff(i)
   enddo
   print *,'diss tot:  ',sum_tot
   print *,'diss tot2: ',sum_tot2
   print *,'sum_f:     ',sum_f
   print *,'sum_diff:  ',sum_diss
   print *,'min: ',minval(edot_r(0:iwave)) 
   print *,'max: ',maxval(edot_r(0:iwave)) 
   

   ymax=maxval(spec_r2)
   ymax=max(ymax,-minval(spec_r2))
   
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
   call plotascii(spec_r2,iwave,message(1:25),-ymax2,ymax2)
#endif   

endif
end subroutine









subroutine compute_spectrum(pin,p,work,spectrum,spec_d,spectrum_x,spectrum_y,&
   spectrum_z,iwave_max,skip_fft)
!
!  INPUT:  iwave_max:  size of spectrum()
!  OUTPUT: iwave:      number of coefficients returned in spectrum()
!          spectrum()  spherical wave number spectrum
!          spec_d()    spherical wave number spectrum of diffusion term
!          spectrum_x  spectrum in x
!          spectrum_y  spectrum in y
!          spectrum_z  spectrum in z
!
!  skip_fft=0    pin = grid point data - take FFT
!  skip_fft=1    pin = FFT data, skip the FFT
!
use params
use mpi
implicit none
integer :: iwave_max,ierr,skip_fft
real*8 :: pin(nx,ny,nz)
real*8 :: work(nx,ny,nz)
real*8 :: p(nx,ny,nz)
real*8 :: spectrum(0:iwave_max)
real*8 :: spec_d(0:iwave_max)
real*8 :: spectrum_x(0:g_nx/2)
real*8 :: spectrum_y(0:g_ny/2)
real*8 :: spectrum_z(0:g_nz/2)

! local variables
real*8 rwave
real*8 :: spectrum_in(0:max(g_nx,g_ny,g_nz))
real*8 :: energy,denergy,xfac,xw
integer i,j,k,iw


rwave=Lz*sqrt(  (g_nx/2.0)**2 + (g_ny/2.0)**2 + (g_nz/(2.0*Lz))**2 )
if (nint(rwave)>iwave_max) then
   call abortdns("compute_spectrum: called with insufficient storege for spectrum()")
endif
iwave_max=nint(rwave)


p=pin
if (skip_fft==0) call fft3d(p,work)
spectrum=0
spec_d=0
spectrum_x=0
spectrum_y=0
spectrum_z=0

do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
    rwave = imcord(i)**2 + jmcord(j)**2 + (kmcord(k)/Lz)**2
    iw = nint(Lz*sqrt(rwave))

    xfac = 8
    if (kmcord(k)==0) xfac=xfac/2
    if (jmcord(j)==0) xfac=xfac/2
    if (imcord(i)==0) xfac=xfac/2
    energy=xfac*p(i,j,k)*p(i,j,k)

    spectrum(iw)=spectrum(iw)+energy
    spectrum_x(abs(imcord(i)))=spectrum_x(abs(imcord(i))) + energy
    spectrum_y(abs(jmcord(j)))=spectrum_y(abs(jmcord(j))) + energy
    spectrum_z(abs(kmcord(k)))=spectrum_z(abs(kmcord(k))) + energy


    xw=rwave*pi2_squared
    if (mu_hyper==4) then
       xw=xw**4  ! viscosity = (del**2)**mu_hyper
    endif
    spec_d(iw)=spec_d(iw)  -mu*xw*energy


enddo
enddo
enddo


#ifdef USE_MPI
spectrum_in=spectrum
call mpi_reduce(spectrum_in,spectrum,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spectrum_in(0:g_nx/2)=spectrum_x
call mpi_reduce(spectrum_in,spectrum_x,1+(g_nx/2),MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spectrum_in(0:g_ny/2)=spectrum_y
call mpi_reduce(spectrum_in,spectrum_y,1+(g_ny/2),MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spectrum_in(0:g_nz/2)=spectrum_z
call mpi_reduce(spectrum_in,spectrum_z,1+(g_nz/2),MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)

spectrum_in=spec_d
call mpi_reduce(spectrum_in,spec_d,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)


#endif

! compute maximum wave number of complete spherical shell
! max wave number in shell integer k = (Lz*im, Lz*jm, km)
!
! (k/Lz) = im < g_nx/2          k < Lz*g_nx/2
! k = km < g_nz/2      
!
! 
if (g_nz == 1)  then
   iwave = min(g_nx/2,g_ny/2)
else
   iwave = floor(min(Lz*g_nx/2d0,Lz*g_ny/2d0,g_nz/2d0))
endif

! for all waves outside sphere, sum into one wave number:
do i=iwave+2,iwave_max
   spectrum(iwave+1)=spectrum(iwave+1)+spectrum(i)
   spec_d(iwave+1)=spec_d(iwave+1)+spec_d(i)
enddo
iwave=iwave+1



end subroutine








subroutine compute_spectrum_2d(pin,p,work,spectrum,iwave_max,skip_fft)
!
!  INPUT:  iwave_max:  size of spectrum()
!  OUTPUT: iwave:      number of coefficients returned in spectrum()
!          spectrum()  spherical wave number spectrum
!
!
!  skip_fft=0    pin = grid point data - take FFT
!  skip_fft=1    pin = FFT data, skip the FFT
!
use params
use mpi
implicit none
integer :: iwave_max,ierr,skip_fft
real*8 :: pin(nx,ny,nz)
real*8 :: work(nx,ny,nz)
real*8 :: p(nx,ny,nz)
real*8 :: spectrum(0:max(g_nx,g_ny),0:g_nz/2)


! local variables
real*8 rwave
real*8 :: spectrum_in(0:max(g_nx,g_ny),0:g_nz/2)
real*8 :: energy,denergy,xfac,xw
integer ::  i,j,k,n,km

rwave=sqrt(  (g_nx/2.0)**2 + (g_ny/2.0)**2 )
if (nint(rwave)>iwave_max) then
   call abortdns("compute_spectrum_2d: called with insufficient storege for spectrum()")
endif
iwave_max=nint(rwave)


p=pin
if (skip_fft==0) call fft3d(p,work)
spectrum=0


do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2

   !jm=jmcord(j)
   !im=imcord(i)
   km=abs(kmcord(k))
   
   rwave = imcord(i)**2 + jmcord(j)**2 
   iwave_2d = nint(sqrt(rwave))
   
   xfac = 8
   if (kmcord(k)==0) xfac=xfac/2
   if (jmcord(j)==0) xfac=xfac/2
   if (imcord(i)==0) xfac=xfac/2
   energy=xfac*p(i,j,k)*p(i,j,k)
   
   spectrum(iwave_2d,km)=spectrum(iwave_2d,km)+energy

enddo
enddo
enddo


#ifdef USE_MPI
spectrum_in=spectrum
n=(1+max(g_nx,g_ny))*(1+g_nz/2)
call mpi_reduce(spectrum_in,spectrum,n,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
#endif

iwave_2d = min(g_nx/2,g_ny/2)
! for all waves outside sphere, sum into one wave number:
do i=iwave_2d+2,iwave_max
   spectrum(iwave_2d+1,:)=spectrum(iwave_2d+1,:)+spectrum(i,:)
enddo
iwave_2d=iwave_2d+1



end subroutine













subroutine compute_spectrum_z_fft(p1,p2,spec)
use params
use mpi
implicit none
integer :: ierr
real*8 :: p1(g_nz2,nx_2dz,ny_2dz)
real*8 :: p2(g_nz2,nx_2dz,ny_2dz)
real*8 :: spec(0:max(g_nx,g_ny,g_nz))

! local variables
real*8 rwave
real*8 :: spec_r_in(0:max(g_nx,g_ny,g_nz))
real*8 :: energy
integer i,j,k,jm,km,im,iwave_max,iw


rwave=Lz*sqrt(  (g_nx/2.0)**2 + (g_ny/2.0)**2 + (g_nz/(2.0*Lz))**2 )
iwave_max=nint(rwave)
if (iwave_max>max(g_nx,g_ny,g_nz)) then
   call abortdns("compute_spectrum_z_fft: called with insufficient storege for spectrum()")
endif

spec=0

do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nx_2dz
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

         rwave = im**2 + jm**2 + (km/Lz)**2
         iw = nint(Lz*sqrt(rwave))
         
         energy = 8
         if (km==0) energy=energy/2
         if (jm==0) energy=energy/2
         if (im==0) energy=energy/2
         energy=energy*p1(k,i,j)*p2(k,i,j)
         
         spec(iw)=spec(iw)+energy
         
enddo
enddo
enddo


#ifdef USE_MPI
spec_r_in=spec
call mpi_reduce(spec_r_in,spec,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
#endif

! compute maximum wave number of complete spherical shell
! max wave number in shell integer k = (Lz*im, Lz*jm, km)
!
! (k/Lz) = im < g_nx/2          k < Lz*g_nx/2
! k = km < g_nz/2      
!
! 
if (g_nz == 1)  then
   iwave = min(g_nx/2,g_ny/2)
else
   iwave = floor(min(Lz*g_nx/2d0,Lz*g_ny/2d0,g_nz/2d0))
endif


! for all waves outside sphere, sum into one wave number:
do i=iwave+2,iwave_max
   spec(iwave+1)=spec(iwave+1)+spec(i)
enddo
iwave=iwave+1

end subroutine




subroutine compute_spectrum_fft(p1,p2,spec)
use params
use mpi
implicit none
integer :: ierr
real*8 :: p1(nx,ny,nz)
real*8 :: p2(nx,ny,nz)
real*8 :: spec(0:max(g_nx,g_ny,g_nz))

! local variables
real*8 rwave
real*8 :: spec_r_in(0:max(g_nx,g_ny,g_nz))
real*8 :: energy
integer i,j,k,jm,km,im,iwave_max,iw


if (Lz/=1) call abortdns("Error: compute_spetrum_fft cant handle Lz<>1")
rwave=sqrt(  (g_nx/2.0)**2 + (g_ny/2.0)**2 + (g_nz/2.0)**2 )
iwave_max=nint(rwave)
if (iwave_max>max(g_nx,g_ny,g_nz)) then
   call abortdns("compute_spectrum_z_fft: called with insufficient storege for spectrum()")
endif


spec=0

do j=ny1,ny2
   jm=jmcord(j)
   do i=nx1,nx2
      im=imcord(i)
      do k=nz1,nz2
         km=kmcord(k)

         rwave = im**2 + jm**2 + (km/Lz)**2
         iw = nint(Lz*sqrt(rwave))
         
         energy = 8
         if (km==0) energy=energy/2
         if (jm==0) energy=energy/2
         if (im==0) energy=energy/2
         energy=energy*p1(i,j,k)*p2(i,j,k)
         
         spec(iw)=spec(iw)+energy
         
enddo
enddo
enddo


#ifdef USE_MPI
spec_r_in=spec
call mpi_reduce(spec_r_in,spec,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
#endif

!max wave number in sphere.
if (g_nz == 1)  then
   iwave = min(g_nx/2,g_ny/2)
else
   iwave = floor(min(Lz*g_nx/2d0,Lz*g_ny/2d0,g_nz/2d0))
endif


! for all waves outside sphere, sum into one wave number:
do i=iwave+2,iwave_max
   spec(iwave+1)=spec(iwave+1)+spec(i)
enddo
iwave=iwave+1

end subroutine



subroutine compute_helicity_spectrum(Q,p1,work,skip_fft)
!
! skip_fft=1:
!       input: Q    p1, work are work arrays
! skip_fft=0:
!      input: p1 (which should be Qhat).  Q is not used.  
!
!
use params
use mpi
implicit none
integer :: ierr,skip_fft
real*8 :: Q(nx,ny,nz,*)
real*8 :: p1(nx,ny,nz,*)
real*8 :: work(nx,ny,nz)

! local variables
real*8 rwave
real*8 :: spec_r_in(0:max(g_nx,g_ny,g_nz))
real*8 ::  spec_x_in(0:g_nx/2,n_var)   
real*8 ::  spec_y_in(0:g_ny/2,n_var)
real*8 ::  spec_z_in(0:g_nz/2,n_var)

real*8 :: energy,vx,wx,uy,wy,uz,vz,heltot
real*8 :: diss1,diss2,hetot,co_energy(3),xw,xfac
integer i,j,k,jm,km,im,iwave_max,n,iw

if (ndim<3) then
   call abortdns("compute_helicity_specturm: can only be used in 3D")
endif

if (skip_fft==0) then
p1(:,:,:,1:3)=Q(:,:,:,1:3)
do n=1,3
   call fft3d(p1(1,1,1,n),work)
enddo
endif

rwave=Lz*sqrt(  (g_nx/2.0)**2 + (g_ny/2.0)**2 + (g_nz/(2.0*Lz))**2 )
iwave_max=nint(rwave)
if (iwave_max > max(g_nx,g_ny,g_nz)) then
   call abortdns("compute_helicity_spec(): called with insufficient storege for spectrum()")
endif


hetot=0
diss1=0
diss2=0
spec_helicity_rp=0
spec_helicity_rn=0
cospec_x=0
cospec_y=0
cospec_z=0
cospec_r=0
spec_kEk=0


do j=ny1,ny2
   jm=jmcord(j)
   do i=nx1,nx2
      im=imcord(i)
      do k=nz1,nz2
         km=kmcord(k)

         rwave = im**2 + jm**2 + (km/Lz)**2
         iw = nint(Lz*sqrt(rwave))
         
         !ux = - pi2*im*p1(i+imsign(i),j,k,1)
         vx = - pi2*im*p1(i+imsign(i),j,k,2)
         wx = - pi2*im*p1(i+imsign(i),j,k,3)
         
         uy = - pi2*jm*p1(i,j+jmsign(j),k,1)
         !vy = - pi2*jm*p1(i,j+jmsign(j),k,2)
         wy = - pi2*jm*p1(i,j+jmsign(j),k,3)
         
         uz =  - pi2*km*p1(i,j,k+kmsign(k),1)/Lz
         vz =  - pi2*km*p1(i,j,k+kmsign(k),2)/Lz
         !wz =  - pi2*km*p1(i,j,k+kmsign(k),3)/Lz
      
         ! vorcity: ( (wy - vz), (uz - wx), (vx - uy) )

         xfac = 8
         if (km==0) xfac=xfac/2
         if (jm==0) xfac=xfac/2
         if (im==0) xfac=xfac/2

         ! compute k E(k)
         xw=sqrt(rwave*pi2_squared)
         spec_kEk(iw)=spec_kEk(iw)+xw*xfac* &
            (p1(i,j,k,1)**2 + p1(i,j,k,2)**2 + p1(i,j,k,3)**2)

         co_energy(1) = xfac*p1(i,j,k,1)*p1(i,j,k,2)
         co_energy(2) = xfac*p1(i,j,k,1)*p1(i,j,k,3)
         co_energy(3) = xfac*p1(i,j,k,2)*p1(i,j,k,3)

         cospec_x(abs(im),1:ndim)=cospec_x(abs(im),1:ndim)+co_energy(1:ndim)  ! uv
         cospec_y(abs(jm),1:ndim)=cospec_y(abs(jm),1:ndim)+co_energy(1:ndim)  ! uw
         cospec_z(abs(km),1:ndim)=cospec_z(abs(km),1:ndim)+co_energy(1:ndim)  ! vw
         cospec_r(iw,1:ndim)=cospec_r(iw,1:ndim)+co_energy(1:ndim)

         energy = xfac*(p1(i,j,k,1)*(wy-vz) + &
                          p1(i,j,k,2)*(uz-wx) + &
                          p1(i,j,k,3)*(vx-uy)) 
         if (energy>0) spec_helicity_rp(iw)=spec_helicity_rp(iw)+energy
         if (energy<0) spec_helicity_rn(iw)=spec_helicity_rn(iw)+energy



         hetot=hetot+energy
         diss1=diss1 -2*energy*iw**2*pi2_squared
         diss2=diss2 -2*energy*rwave*pi2_squared
enddo
enddo
enddo


#ifdef USE_MPI
spec_r_in=spec_helicity_rp
call mpi_reduce(spec_r_in,spec_helicity_rp,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=spec_helicity_rn
call mpi_reduce(spec_r_in,spec_helicity_rn,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=spec_kEk
call mpi_reduce(spec_r_in,spec_kEk,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)

do n=1,ndim
   spec_r_in=cospec_r(:,n)
   call mpi_reduce(spec_r_in,cospec_r(0,n),(1+iwave_max),MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
enddo

spec_x_in=cospec_x
call mpi_reduce(spec_x_in,cospec_x,(1+g_nx/2)*ndim,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_y_in=cospec_y
call mpi_reduce(spec_y_in,cospec_y,(1+g_ny/2)*ndim,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_z_in=cospec_z
call mpi_reduce(spec_z_in,cospec_z,(1+g_nz/2)*ndim,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)


rwave=diss1
call mpi_reduce(rwave,diss1,1,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
rwave=diss2
call mpi_reduce(rwave,diss2,1,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
rwave=hetot
call mpi_reduce(rwave,hetot,1,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
#endif

#if 0
if (my_pe==io_pe) then
   print *,'helicity data from compute_helicity_spec()'
   print *,'helicity: ',hetot
   print *,'helicity dissipation (spectrum): ',diss1*mu
   print *,'helicity dissipation (exact):    ',diss2*mu
endif
#endif

! compute maximum wave number of complete spherical shell
! max wave number in shell integer k = (Lz*im, Lz*jm, km)
!
! (k/Lz) = im < g_nx/2          k < Lz*g_nx/2
! k = km < g_nz/2      
!
! 
if (g_nz == 1)  then
   iwave = min(g_nx/2,g_ny/2)
else
   iwave = floor(min(Lz*g_nx/2d0,Lz*g_ny/2d0,g_nz/2d0))
endif


! for all waves outside sphere, sum into one wave number:
do i=iwave+2,iwave_max
   spec_helicity_rp(iwave+1)=spec_helicity_rp(iwave+1)+spec_helicity_rp(i)
   spec_helicity_rn(iwave+1)=spec_helicity_rn(iwave+1)+spec_helicity_rn(i)
   cospec_r(iwave+1,1:ndim)=cospec_r(iwave+1,1:ndim)+cospec_r(i,1:ndim)
   spec_kEk(iwave+1)=spec_kEk(iwave+1)+spec_kEk(i)

enddo
iwave=iwave+1

if (my_pe==io_pe) then
   heltot=0
   do i=0,iwave
      heltot=heltot+spec_helicity_rp(i)+spec_helicity_rn(i)
   enddo
   print *,'total helicity: ',heltot
endif
   
end subroutine




subroutine compute_hfree_spec(pgrid,cmodes_r,cmodes_i,p1)
!
!various manipulations of helicity (sk)

use params
use mpi
implicit none
integer :: ierr,skip_fft
real*8 :: pgrid(nx,ny,nz,3)
real*8 :: cmodes_r(nx,ny,nz,3)
real*8 :: cmodes_i(nx,ny,nz,3)
real*8 :: p1(nx,ny,nz,3)


! local variables
real*8 rwave
real*8 :: spec_r_in(0:max(g_nx,g_ny,g_nz))
real*8 :: countn(0:max(g_nx,g_ny,g_nz)), countp(0:max(g_nx,g_ny,g_nz)),pcount(0:max(g_nx,g_ny,g_nz))
real*8 ::  spec_x_in(0:g_nx/2,n_var)   
real*8 ::  spec_y_in(0:g_ny/2,n_var)
real*8 ::  spec_z_in(0:g_nz/2,n_var)
real*8 :: energy,vx,wx,uy,wy,uz,vz,heltot
real*8 :: diss1,diss2,hetot,co_energy(3),xw,RR(3),II(3),mod_rr,mod_ii
real*8 :: WR(3),WI(3)
real*8 :: cos_tta,rhel,delta,rp,ip,e1,e2,xfac,maxct,minct,rhmin,rhmax
integer i,j,k,jm,km,im,iwave_max,n,ibin,a,b,ra,rb,ind,iw

complex*16 tmp


if (Lz/=1) call abortdns("Error: compute_hfree_spec cant handle Lz<>1")

if (ndim<3) then
   call abortdns("compute_hfree_spec: can only be used in 3D")
endif

#undef TESTEXP
#ifdef TESTEXP
! 5,-3,1   xfac=8   efac=64    N=64
! 0,-3,1   xfac=4   efac=32
! 0,0,1   xfac=2   efac=16
! 0,0,0   xfac=1   efac=8
im=5
jm=-3
km=1
print *,'initial mode: im,jm,km: ',im,jm,km

pgrid=0
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   tmp=im*pi2*(0,1)*xcord(i) + &
       jm*pi2*(0,1)*ycord(j) + &
       km*pi2*(0,1)*zcord(k) 
   tmp=exp(tmp)
   pgrid(i,j,k,1)=real(tmp)
enddo
enddo
enddo
#endif

p1=pgrid
! compute fft in p1. (use cmodes_r as work array)
do n = 1,3
   call fft3d(p1(1,1,1,n),cmodes_r)
enddo
do n = 1,3
   call sincos_to_complex_field(p1(1,1,1,n),cmodes_r(1,1,1,n),cmodes_i(1,1,1,n))
enddo

#ifdef TESTEXP
rp=0
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   rp=rp+pgrid(i,j,k,1)**2
enddo
enddo
enddo
print *,'energy (grid space) = ',rp/g_nx/g_ny/g_nz
e1=0
e2=0
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   ! wave number (im,jm,km)   im positive or negative
   im=imcord(i)
   jm=jmcord(j)
   km=kmcord(k)
   if (abs(p1(i,j,k,1))>1e-14 ) then
      write(*,'(3i4,f10.5)') im,jm,km,p1(i,j,k,1)
   endif

   xfac=8
   if (km==0) xfac=xfac/2
   if (jm==0) xfac=xfac/2
   if (im==0) xfac=xfac/2
   e1=e1+xfac*p1(i,j,k,1)**2
enddo
enddo
enddo

do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   ! wave number (im,jm,km)   im positive or negative
   im=imcord_exp(i)
   jm=jmcord_exp(j)
   km=kmcord_exp(k)
   rp = cmodes_r(i,j,k,1)
   ip = cmodes_i(i,j,k,1)
   if (abs(rp)>1e-14 .or. abs(ip)>1e-14  ) then
      write(*,'(3i4,a,f10.5,a,f10.5,a,f10.5)') im,jm,km,&
      '  (',rp,',',ip,')  '
   endif

   xfac=64
   if (km==0) xfac=xfac/2
   if (jm==0) xfac=xfac/2
   if (im==0) xfac=xfac/2
   e2=e2+xfac*(rp**2+ip**2)
enddo
enddo
enddo
print *,'energy (sin/cos) = ',e1
print *,'energy (complex) = ',e2
!return
stop
#endif


rwave=sqrt(  (g_nx/2.0)**2 + (g_ny/2.0)**2 + (g_nz/2.0)**2 )
iwave_max=nint(rwave)
e1=0
e2=0
hetot=0
diss1=0
diss2=0
spec_helicity_rp=0
spec_helicity_rn=0
spec_meanhe=0
spec_varhe=0
spec_E = 0
spec_kEk=0
spec_varE=0

#if 0

E33 = 0
II2sq = 0
RRsq = 0
I2I3  = 0
R2I3 = 0
par = 0

#endif

cos_tta_spec_n = 0
cos_tta_spec_p = 0
rhel_spec_n = 0
rhel_spec_p = 0
costta_pdf = 0
rhel_pdf = 0
countp = 0
countn = 0


do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   xfac=8
   jm=jmcord(j)
   im=imcord(i)
   km=kmcord(k)
   if (km==0) xfac=xfac/2
   if (jm==0) xfac=xfac/2
   if (im==0) xfac=xfac/2
   e1=e1+.5*xfac*(p1(i,j,k,1)**2 + p1(i,j,k,2)**2 + p1(i,j,k,3)**2)
enddo
enddo
enddo


! for cos_tta pdf, bin index = a + b*cos_tta; calculate a and b
minct = -1
maxct = 1
a = (nbin*minct - maxct)/(minct - maxct)
b = (1 - nbin)/(minct - maxct)

! for rhel pdf, bin index = ra + rb*rhel; calculate ra and rb
rhmin = -1 
rhmax = 1
ra = (nbin*rhmin - rhmax)/(rhmin - rhmax)
rb = (1-nbin)/(rhmin - rhmax)


do j=ny1,ny2
   jm=jmcord_exp(j)
   do i=nx1,nx2
      im=imcord_exp(i)
      do k=nz1,nz2
         km=kmcord_exp(k)
         
         rwave = im**2 + jm**2 + (km/Lz)**2
         iw = nint(Lz*sqrt(rwave))

         !     compute angle between Re and Im parts of u(k)
         RR = cmodes_r(i,j,k,:)
         II = cmodes_i(i,j,k,:)
         
         mod_rr = sqrt(RR(1)*RR(1) + RR(2)*RR(2) + RR(3)*RR(3))
         mod_ii = sqrt(II(1)*II(1) + II(2)*II(2) + II(3)*II(3))
         
         if (iw > 0) then
            cos_tta = (RR(1)*II(1) + RR(2)*II(2) + RR(3)*II(3))/&
                 (mod_rr*mod_ii)
            
            ! histogram of cosine of angle between RR and II 
!            ind = nint(a + b*abs(cos_tta))	
            ind = nint(a + b*(cos_tta))
            
            if (ind>nbin) then 
               write(6,*)"iw,cos(tta),ind = ",iw,cos_tta,ind
               call abortdns("Error: spectrum.F90: ind>nbin")
               
            endif
            if (ind<1) then
               write(6,*)"iw,cos(tta),ind = ",iw,cos_tta,ind
               call abortdns("Error: spectrum.F90: ind<1")            
            endif

            costta_pdf(iw,ind) = costta_pdf(iw,ind) + 1.0        
            
            ! mean value of angles in each wavenumber        
!            cos_tta_spec(iw) = cos_tta_spec(iw) + abs(cos_tta)
            if (cos_tta >= 0) then
               cos_tta_spec_p(iw) = cos_tta_spec_p(iw) + (cos_tta)
               countp(iw) = countp(iw)+1	
            else 
               cos_tta_spec_n(iw) = cos_tta_spec_n(iw) + (cos_tta)
               countn(iw) = countn(iw)+1

            endif
           endif
           
         
#if 0
         ! compute vorticity           
         ! sqrt(-1) * 2pi * (im,jm,km) cross (RR-sqrt(-1)II)
         WR(1) = pi2*(-jm*II(3)+km*II(2)/Lz)  
         WR(2) = pi2*(im*II(3) - km*II(1)/Lz)
         WR(3) = pi2*(-im*II(2) + jm*II(1))
         WI(1) = pi2*(-jm*RR(3) + km*RR(2)/Lz)  
         WI(2) = pi2*(im*RR(3) - km*RR(1)/Lz)
         WI(3) = pi2*(-im*RR(2) + jm*RR(1))	
#endif            
         xfac = 64
         if (km==0) xfac=xfac/2
         if (jm==0) xfac=xfac/2
         if (im==0) xfac=xfac/2
         
         !     compute E(k) and kE(k)
         xw=sqrt(rwave*pi2_squared)
         e2 = e2 + .5*xfac*(sum(RR*RR)+ sum(II*II))
         
         ! these are turned off because they don't seem to give anything useful
#if 0
         ! In coordinate system with x_1 || k, x_2 || RR and x_3 ||(k \cross RR):
         
         ! E_33 component of energy tensor (= I_3^2 = square of component of II orthogonal to RR)
         E33(iw)=E33(iw)+ 0.5*xfac*(mod_ii**2 - (sum(II*RR)/(mod_rr))**2)
         
         ! I_2^2 (square of the component of II along RR) (presents zero contribution to helicity?)
         II2sq(iw) = II2sq(iw) + 0.5*xfac*(sum(II*RR)/(mod_rr))**2

         ! RR^2 (square of RR) 
         RRsq(iw) = RRsq(iw)+0.5*xfac*sum(RR*RR)

         ! I_2*I_3 (the part of the energy tensor that doesn't contribute to 
         ! either energy or helicity
         I2I3(iw) = I2I3(iw) + 0.5*xfac*((sum(II*RR)/(mod_rr))*&
              sqrt(mod_ii**2 -(sum(II*RR)/(mod_rr))**2))

         ! R_2*I_3 = RR*I_3 (the part of the energy tensor that contributes 
         ! to the helicity
         R2I3(iw) = R2I3(iw) + 0.5*xfac*(mod_rr*sqrt(mod_ii**2 -(sum(II*RR)/(mod_rr))**2))
         
         !par = abs(RR\cross II)/(RR^2+II^2)
         if (iw > 0) then
         par(iw) = par(iw) + sqrt((RR(2)*II(3) - II(2)*RR(3))**2 + &
              (II(1)*RR(3) - RR(1)*II(3))**2 + &
              (RR(1)*II(2) - II(1)*RR(2))**2)/(mod_rr**2 + mod_ii**2)
         endif
#endif

         !	helicity(k) = k\cdot RR(k) cross II(k)            
         energy = xfac * 2 * pi2 * (im*(RR(2)*II(3) - II(2)*RR(3)) + &
              jm*(II(1)*RR(3) - RR(1)*II(3)) + &
              (km/Lz)*(RR(1)*II(2) - II(1)*RR(2)))
         	
         ! relative helicity and its distribution in current wavenumber
         if (iw > 0) then
            rhel = energy/(xw*xfac*(sum(RR*RR)+sum(II*II)))
!            write(6,*)"rhel = ",rhel
!            if (rhel < rhmin) rhmin = rhel
!            if (rhel > rhmax) rhmax = rhel
            
            !  histogram of relative helicity in each shell
            ind = nint(ra + rb*(rhel))	   
	
            
            if (ind>nbin) then 
               write(6,*)"iw,rhel,ind = ",iw,rhel,ind
               call abortdns("Error: spectrum.F90: ind>nbin")
               
            endif
            if (ind<1) then 
               write(6,*)"iw,rhel,ind = ",iw,rhel,ind
               call abortdns("Error: spectrum.F90: ind<1")
               
            endif
            
            rhel_pdf(iw,ind) = rhel_pdf(iw,ind) + 1        
            
            !distribution of total relative helicity in each wavenumber
            if (rhel >= 0) then
               rhel_spec_p(iw) = rhel_spec_p(iw) + rhel
            else 
               rhel_spec_n(iw) = rhel_spec_n(iw) + rhel
	       
            endif

            !spectrum of variance of relative helicity
            rhel_rms_spec(iw) = rhel_rms_spec(iw) + rhel**2	            
            
         endif
         
         
         !     cutoff for recalculating the spectra
         delta = 1      !this value can be changed by hand
         
         !     omit modes where cos_tta is less than cutoff delta 
         !         (we are looking for 'non-helical' modes)
         
         !	if (abs(cos_tta) > delta) then
         
         !if (abs(rhel) <= delta) then
         if (.true.) then
            ! store E(k),kE(k), varE(k)
            spec_E(iw)=spec_E(iw) + 0.5*xfac*(sum(RR*RR)+ sum(II*II))
            spec_kEk(iw)=spec_kEk(iw) + xw*xfac*(sum(RR*RR)+ sum(II*II))
            spec_varE(iw) = spec_varE(iw) + (0.5*xfac*(sum(RR*RR)+ sum(II*II)))**2

            ! store helicity(k), mean helicity spec_meanhe(k) 
	    ! and spec_varh(k) (variance of the helicity in each wavenumber)
            if (energy>0) spec_helicity_rp(iw)= & 
                 spec_helicity_rp(iw)+energy
            if (energy<0) spec_helicity_rn(iw)= &
                 spec_helicity_rn(iw) + energy
            spec_meanhe(iw) = spec_meanhe(iw) + energy
            spec_varhe(iw) = spec_varhe(iw) + energy**2
            pcount(iw) = pcount(iw) + 1
       
            hetot=hetot+energy
            diss1=diss1 -2*energy*iw**2*pi2_squared
            diss2=diss2 -2*energy*rwave*pi2_squared  
         endif
      enddo
   enddo
enddo




#ifdef USE_MPI
spec_r_in=spec_helicity_rp
call mpi_reduce(spec_r_in,spec_helicity_rp,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=spec_helicity_rn
call mpi_reduce(spec_r_in,spec_helicity_rn,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=spec_meanhe
call mpi_reduce(spec_r_in,spec_meanhe,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=spec_varhe
call mpi_reduce(spec_r_in,spec_varhe,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=spec_E
call mpi_reduce(spec_r_in,spec_E,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=spec_kEk
call mpi_reduce(spec_r_in,spec_kEk,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=spec_varE
call mpi_reduce(spec_r_in,spec_varE,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)


! these quantities correspond to those disabled above. Make sure to turn those on before turning this on.
#if 0
spec_r_in=E33
call mpi_reduce(spec_r_in,E33,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=II2sq
call mpi_reduce(spec_r_in,II2sq,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=RRsq
call mpi_reduce(spec_r_in,RRsq,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=I2I3
call mpi_reduce(spec_r_in,I2I3,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=R2I3
call mpi_reduce(spec_r_in,R2I3,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=par
call mpi_reduce(spec_r_in,par,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
#endif

spec_r_in = cos_tta_spec_n
call mpi_reduce(spec_r_in,cos_tta_spec_n,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)

spec_r_in = cos_tta_spec_p
call mpi_reduce(spec_r_in,cos_tta_spec_p,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)

spec_r_in = countn
call mpi_reduce(spec_r_in,countn,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)

spec_r_in = countp
call mpi_reduce(spec_r_in,countp,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)

spec_r_in = rhel_spec_n
call mpi_reduce(spec_r_in,rhel_spec_n,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)

spec_r_in = rhel_spec_p
call mpi_reduce(spec_r_in,rhel_spec_p,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)

spec_r_in = rhel_rms_spec
call mpi_reduce(spec_r_in,rhel_rms_spec,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)


do n = 1,nbin
   spec_r_in = costta_pdf(:,n)
   call mpi_reduce(spec_r_in,costta_pdf(0,n),(1+iwave_max),MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
enddo

do n = 1,nbin
   spec_r_in = rhel_pdf(:,n)
   call mpi_reduce(spec_r_in,rhel_pdf(0,n),(1+iwave_max),MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
enddo


do n=1,ndim
   spec_r_in=cospec_r(:,n)
   call mpi_reduce(spec_r_in,cospec_r(0,n),(1+iwave_max),MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
enddo


spec_x_in=cospec_x
call mpi_reduce(spec_x_in,cospec_x,(1+g_nx/2)*ndim,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_y_in=cospec_y
call mpi_reduce(spec_y_in,cospec_y,(1+g_ny/2)*ndim,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_z_in=cospec_z
call mpi_reduce(spec_z_in,cospec_z,(1+g_nz/2)*ndim,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)


rwave=diss1
call mpi_reduce(rwave,diss1,1,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
rwave=diss2
call mpi_reduce(rwave,diss2,1,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
rwave=hetot
call mpi_reduce(rwave,hetot,1,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
rwave=e1
call mpi_reduce(rwave,e1,1,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
rwave=e2
call mpi_reduce(rwave,e2,1,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
#endif



if (my_pe==io_pe) then
   print *,'KE: (from sin/cos modes)',e1
   print *,'KE: (from complex modes)',e2
   print *,'helicity: ',hetot
   print *,'helicity dissipation (spectrum): ',diss1*mu
   print *,'helicity dissipation (exact):    ',diss2*mu
endif

! maximum wave number for spherical shells:
if (g_nz == 1)  then
   iwave = min(g_nx/2,g_ny/2)
else
   iwave = floor(min(Lz*g_nx/2d0,Lz*g_ny/2d0,g_nz/2d0))
endif

! for all waves outside sphere, sum into one wave number, iwave+1:
do i=iwave+2,iwave_max
   spec_helicity_rp(iwave+1)=spec_helicity_rp(iwave+1)+spec_helicity_rp(i)
   spec_helicity_rn(iwave+1)=spec_helicity_rn(iwave+1)+spec_helicity_rn(i)
   spec_meanhe(iwave+1)=spec_meanhe(iwave+1)+spec_meanhe(i)
   spec_varhe(iwave+1)=spec_varhe(iwave+1)+spec_varhe(i)
   spec_E(iwave+1)=spec_E(iwave+1)+spec_E(i)  
   spec_kEk(iwave+1)=spec_kEk(iwave+1)+spec_kEk(i)
   spec_varE(iwave+1)=spec_varE(iwave+1)+spec_varE(i)
   E33(iwave+1)=E33(iwave+1)+E33(i)
   II2sq(iwave+1)=II2sq(iwave+1)+II2sq(i)
   RRsq(iwave+1)=RRsq(iwave+1)+ RRsq(i)
   I2I3(iwave+1)=I2I3(iwave+1)+I2I3(i)
   R2I3(iwave+1)=R2I3(iwave+1)+R2I3(i)
   par(iwave+1)=par(iwave+1)+par(i)		   
   cos_tta_spec_n(iwave+1) = cos_tta_spec_n(iwave+1) + cos_tta_spec_n(i)
   cos_tta_spec_p(iwave+1) = cos_tta_spec_p(iwave+1) + cos_tta_spec_p(i)
   rhel_spec_n(iwave+1) = rhel_spec_n(iwave+1) + rhel_spec_n(i)
   rhel_spec_p(iwave+1) = rhel_spec_p(iwave+1) + rhel_spec_p(i)
   rhel_rms_spec(iwave+1) = rhel_rms_spec(iwave+1) + rhel_rms_spec(i)
   costta_pdf(iwave+1,:) = costta_pdf(iwave+1,:) + costta_pdf(i,:)
   rhel_pdf(iwave+1,:) = rhel_pdf(iwave+1,:) + rhel_pdf(i,:)

enddo

!average costta and other mean values over sphere (exclude iwave=0)
do i = 1,iwave+1
cos_tta_spec_n(i) = cos_tta_spec_n(i)/countn(i)
cos_tta_spec_p(i) = cos_tta_spec_p(i)/countp(i)
spec_meanhe(i) = (spec_meanhe(i)/pcount(i))
spec_varhe(i) = (spec_varhe(i)/pcount(i))
spec_meanE(i) = (spec_E(i)/pcount(i))
spec_varE(i) = (spec_varE(i)/pcount(i))

!below do not need to be averaged over shells!
!rhel_spec_n(i) = rhel_spec_n(i)/pcountn(i)
!rhel_spec_p(i) = rhel_spec_p(i)/pcountp(i)
!rhel_rms_spec(i) = (rhel_rms_spec(i)/rcount(i))


enddo

!compute pdfs from histograms (exclude iwave = 0)
do i = 1,iwave+1
if (sum(rhel_pdf(i,:)) > 0.0) rhel_pdf(i,:) = rhel_pdf(i,:)/sum(rhel_pdf(i,:))

if (sum(costta_pdf(i,:)) > 0.0) costta_pdf(i,:) = costta_pdf(i,:)/sum(costta_pdf(i,:))

enddo

iwave=iwave+1

if (my_pe==io_pe) then
   heltot=0
   do i=0,iwave
      heltot=heltot+spec_helicity_rp(i)+spec_helicity_rn(i)
   enddo
   print *,'total helicity: ',heltot
endif
   
end subroutine



subroutine compute_shear_cospectrum(Q,p1,work,skip_fft)
!
! skip_fft=1:
!       input: Q    p1, work are work arrays
! skip_fft=0:
!      input: p1 (which should be Qhat).  Q is not used.  
!
!
use params
use mpi
implicit none
integer :: ierr,skip_fft
real*8 :: Q(nx,ny,nz,*)
real*8 :: p1(nx,ny,nz,*)
real*8 :: work(nx,ny,nz)

! local variables
real*8 rwave
real*8 :: spec_r_in(0:max(g_nx,g_ny,g_nz))
real*8 ::  spec_x_in(0:g_nx/2,n_var)   
real*8 ::  spec_y_in(0:g_ny/2,n_var)
real*8 ::  spec_z_in(0:g_nz/2,n_var)

real*8 :: energy,vx,wx,uy,wy,uz,vz,heltot
real*8 :: diss1,diss2,hetot,co_energy(3),xw,xfac
integer i,j,k,jm,km,im,iwave_max,n



end subroutine

end module
