#include "macros.h"
module spectrum
use params
use sforcing
implicit none

#if 0

module to compute spherical and other 1D spectrums and
transfer function spectrums

#endif

logical :: compute_transfer=.false.

real*8,private ::  spec_x(0:g_nx/2,n_var)
real*8,private ::  spec_y(0:g_ny/2,n_var)
real*8,private ::  spec_z(0:g_nz/2,n_var)
real*8,private ::  spec_r(0:max(g_nx,g_ny,g_nz),n_var)
real*8,private ::  spec_r_new(0:max(g_nx,g_ny,g_nz))
real*8,private ::  edot_r(0:max(g_nx,g_ny,g_nz))
real*8,private ::  time_old=-1

real*8,private ::  spec_helicity_rp(0:max(g_nx,g_ny,g_nz))
real*8,private ::  spec_helicity_rn(0:max(g_nx,g_ny,g_nz))

! cospectrum in x,y,z directions.
! n_var=1,3:  uv, uw, vw
real*8,private ::  cospec_r(0:max(g_nx,g_ny,g_nz),n_var)
real*8,private ::  cospec_x(0:g_nx/2,n_var)   
real*8,private ::  cospec_y(0:g_ny/2,n_var)
real*8,private ::  cospec_z(0:g_nz/2,n_var)


integer,private :: iwave=-1



real*8 ::  transfer_comp_time         ! time at which below terms evaluated at:
real*8 ::  spec_E(0:max(g_nx,g_ny,g_nz)) !E(k) from helicity free modes 
real*8 ::  spec_kEk(0:max(g_nx,g_ny,g_nz))  ! k E(k)
real*8 ::  cos_tta_spec(0:max(g_nx,g_ny,g_nz)) !spec of cos_tta betn RR and II 
real*8 ::  spec_diff(0:max(g_nx,g_ny,g_nz))  ! u dot diffusion term
real*8 ::  spec_diff_new(0:max(g_nx,g_ny,g_nz)) 
real*8 ::  spec_f(0:max(g_nx,g_ny,g_nz))     ! u dot forcing term
real*8 ::  spec_model(0:max(g_nx,g_ny,g_nz)) ! spectrum of div(tau) or smagorinsky term
real*8 ::  spec_model_new(0:max(g_nx,g_ny,g_nz)) 
real*8 ::  spec_rhs(0:max(g_nx,g_ny,g_nz))   ! transfer spec of u dot RHS

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
integer :: iwave_max,i,n
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

q1=Q


!
!   pe=pe+.5*grav*Q(i,j,3)**2 - .5*grav*H0**2
!   ke = ke + .5*Q(i,j,3)*Q(i,j,n)**2   n=1,2
!
!

! passive scalars:
do n=np1,np2
   call compute_spectrum(q1(1,1,1,n),work1,work2,spec_r(0,n),spec_r2,&
       spec_x(0,n),spec_y(0,n),spec_z(0,n),iwave_max,0)
enddo


do i=1,ndim
   call fft3d(q1(1,1,1,i),work1)
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
character(len=80) :: message


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


write(message,'(a,f10.4)') " Energy t=",time
call logplotascii(spec_r(0,1),iwave,message(1:25))
!call logplotascii(spec_x,g_nx/2,message)
!call logplotascii(spec_y,g_ny/2,message)
!call logplotascii(spec_z,g_nz/2,message)






if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time_file
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".spec"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "spec_write(): Error opening file errno=",ierr
      call abort(message)
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
      call abort(message)
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
      call abort(message)
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
      call abort(message)
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
      call abort(message)
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


subroutine output_hfree_spec(time,time_file)
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
   message = rundir(1:len_trim(rundir)) &
        // runname(1:len_trim(runname))&
        // message(2:10) // ".hf_spec"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)')"output_hfree_spec(): Error opening file errno=",ierr
      call abort(message)
   endif
   x=1+iwave; call cwrite8(fid,x,1)
   call cwrite8(fid,time,1)
   call cwrite8(fid,spec_helicity_rn,1+iwave)
   call cwrite8(fid,spec_helicity_rp,1+iwave)
   call cwrite8(fid,spec_E,1+iwave)
   call cwrite8(fid,spec_kEk,1+iwave)
   call cwrite8(fid,cos_tta_spec,1+iwave)
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
      call abort(message)
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
real*8 :: spectrum_in(0:iwave_max)
real*8 :: energy,denergy,xfac,xw
integer i,j,k


rwave=sqrt(  (g_nx/2.0)**2 + (g_ny/2.0)**2 + (g_nz/2.0)**2 )
if (nint(rwave)>iwave_max) then
   call abort("compute_spectrum: called with insufficient storege for spectrum()")
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
    rwave = imcord(i)**2 + jmcord(j)**2 + kmcord(k)**2
    iwave = nint(sqrt(rwave))

    xfac = 8
    if (kmcord(k)==0) xfac=xfac/2
    if (jmcord(j)==0) xfac=xfac/2
    if (imcord(i)==0) xfac=xfac/2
    energy=xfac*p(i,j,k)*p(i,j,k)

    spectrum(iwave)=spectrum(iwave)+energy
    spectrum_x(abs(imcord(i)))=spectrum_x(abs(imcord(i))) + energy
    spectrum_y(abs(jmcord(j)))=spectrum_y(abs(jmcord(j))) + energy
    spectrum_z(abs(kmcord(k)))=spectrum_z(abs(kmcord(k))) + energy


    xw=rwave*pi2_squared
    if (mu_hyper==4) then
       xw=xw**4  ! viscosity = (del**2)**mu_hyper
    endif
    spec_d(iwave)=spec_d(iwave)  -mu*xw*energy


enddo
enddo
enddo


#ifdef USE_MPI
spectrum_in=spectrum
call MPI_reduce(spectrum_in,spectrum,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spectrum_in(0:g_nx/2)=spectrum_x
call MPI_reduce(spectrum_in,spectrum_x,1+(g_nx/2),MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spectrum_in(0:g_ny/2)=spectrum_y
call MPI_reduce(spectrum_in,spectrum_y,1+(g_ny/2),MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spectrum_in(0:g_nz/2)=spectrum_z
call MPI_reduce(spectrum_in,spectrum_z,1+(g_nz/2),MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)

spectrum_in=spec_d
call MPI_reduce(spectrum_in,spec_d,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)

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
   spec_d(iwave+1)=spec_d(iwave+1)+spec_d(i)
enddo
iwave=iwave+1



end subroutine













subroutine compute_spectrum_z_fft(p1,p2,spec)
use params
use mpi
implicit none
integer :: ierr
real*8 :: p1(g_nz2,nslabx,ny_2dz)
real*8 :: p2(g_nz2,nslabx,ny_2dz)
real*8 :: spec(0:max(g_nx,g_ny,g_nz))

! local variables
real*8 rwave
real*8 :: spec_r_in(0:max(g_nx,g_ny,g_nz))
real*8 :: energy
integer i,j,k,jm,km,im,iwave_max


rwave=sqrt(  (g_nx/2.0)**2 + (g_ny/2.0)**2 + (g_nz/2.0)**2 )
iwave_max=nint(rwave)

spec=0

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
         energy=energy*p1(k,i,j)*p2(k,i,j)
         
         spec(iwave)=spec(iwave)+energy
         
enddo
enddo
enddo


#ifdef USE_MPI
spec_r_in=spec
call MPI_reduce(spec_r_in,spec,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
#endif

if (g_nz == 1)  then
   iwave = min(g_nx,g_ny)
else
   iwave = min(g_nx,g_ny,g_nz)
endif
iwave = (iwave/2)           ! max wave number in sphere.

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
integer i,j,k,jm,km,im,iwave_max


rwave=sqrt(  (g_nx/2.0)**2 + (g_ny/2.0)**2 + (g_nz/2.0)**2 )
iwave_max=nint(rwave)

spec=0

do j=ny1,ny2
   jm=jmcord(j)
   do i=nx1,nx2
      im=imcord(i)
      do k=nz1,nz2
         km=kmcord(k)

         rwave = im**2 + jm**2 + km**2
         iwave = nint(sqrt(rwave))
         
         energy = 8
         if (km==0) energy=energy/2
         if (jm==0) energy=energy/2
         if (im==0) energy=energy/2
         energy=energy*p1(i,j,k)*p2(i,j,k)
         
         spec(iwave)=spec(iwave)+energy
         
enddo
enddo
enddo


#ifdef USE_MPI
spec_r_in=spec
call MPI_reduce(spec_r_in,spec,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
#endif

if (g_nz == 1)  then
   iwave = min(g_nx,g_ny)
else
   iwave = min(g_nx,g_ny,g_nz)
endif
iwave = (iwave/2)           ! max wave number in sphere.

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
real*8 :: diss1,diss2,hetot,co_energy(3),xw
integer i,j,k,jm,km,im,iwave_max,n

if (ndim<3) then
   call abort("compute_helicity_specturm: can only be used in 3D")
endif

if (skip_fft==0) then
p1(:,:,:,1:3)=Q(:,:,:,1:3)
do n=1,3
   call fft3d(p1(1,1,1,n),work)
enddo
endif

rwave=sqrt(  (g_nx/2.0)**2 + (g_ny/2.0)**2 + (g_nz/2.0)**2 )
iwave_max=nint(rwave)

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

         rwave = im**2 + jm**2 + km**2
         iwave = nint(sqrt(rwave))
         
         !ux = - pi2*im*p1(i+imsign(i),j,k,1)
         vx = - pi2*im*p1(i+imsign(i),j,k,2)
         wx = - pi2*im*p1(i+imsign(i),j,k,3)
         
         uy = - pi2*jm*p1(i,j+jmsign(j),k,1)
         !vy = - pi2*jm*p1(i,j+jmsign(j),k,2)
         wy = - pi2*jm*p1(i,j+jmsign(j),k,3)
         
         uz =  - pi2*km*p1(i,j,k+kmsign(k),1)
         vz =  - pi2*km*p1(i,j,k+kmsign(k),2)
         !wz =  - pi2*km*p1(i,j,k+kmsign(k),3)
      
         ! vorcity: ( (wy - vz), (uz - wx), (vx - uy) )

         energy = 8
         if (km==0) energy=energy/2
         if (jm==0) energy=energy/2
         if (im==0) energy=energy/2

         ! compute k E(k)
         xw=sqrt(rwave*pi2_squared)
         spec_kEk(iwave)=spec_kEk(iwave)+xw*energy* &
            (p1(i,j,k,1)**2 + p1(i,j,k,2)**2 + p1(i,j,k,3)**2)

         co_energy(1) = energy*p1(i,j,k,1)*p1(i,j,k,2)
         co_energy(2) = energy*p1(i,j,k,1)*p1(i,j,k,3)
         co_energy(3) = energy*p1(i,j,k,2)*p1(i,j,k,3)

         cospec_x(abs(im),1:ndim)=cospec_x(abs(im),1:ndim)+co_energy(1:ndim)  ! uv
         cospec_y(abs(jm),1:ndim)=cospec_y(abs(jm),1:ndim)+co_energy(1:ndim)  ! uw
         cospec_z(abs(km),1:ndim)=cospec_z(abs(km),1:ndim)+co_energy(1:ndim)  ! vw
         cospec_r(iwave,1:ndim)=cospec_r(iwave,1:ndim)+co_energy(1:ndim)

         energy = energy*(p1(i,j,k,1)*(wy-vz) + &
                          p1(i,j,k,2)*(uz-wx) + &
                          p1(i,j,k,3)*(vx-uy)) 
         if (energy>0) spec_helicity_rp(iwave)=spec_helicity_rp(iwave)+energy
         if (energy<0) spec_helicity_rn(iwave)=spec_helicity_rn(iwave)+energy



         hetot=hetot+energy
         diss1=diss1 -2*energy*iwave**2*pi2_squared
         diss2=diss2 -2*energy*rwave*pi2_squared
enddo
enddo
enddo


#ifdef USE_MPI
spec_r_in=spec_helicity_rp
call MPI_reduce(spec_r_in,spec_helicity_rp,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=spec_helicity_rn
call MPI_reduce(spec_r_in,spec_helicity_rn,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=spec_kEk
call MPI_reduce(spec_r_in,spec_kEk,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)

do n=1,ndim
   spec_r_in=cospec_r(:,n)
   call MPI_reduce(spec_r_in,cospec_r(0,n),(1+iwave_max),MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
enddo

spec_x_in=cospec_x
call MPI_reduce(spec_x_in,cospec_x,(1+g_nx/2)*ndim,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_y_in=cospec_y
call MPI_reduce(spec_y_in,cospec_y,(1+g_ny/2)*ndim,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_z_in=cospec_z
call MPI_reduce(spec_z_in,cospec_z,(1+g_nz/2)*ndim,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)


rwave=diss1
call MPI_reduce(rwave,diss1,1,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
rwave=diss2
call MPI_reduce(rwave,diss2,1,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
rwave=hetot
call MPI_reduce(rwave,hetot,1,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
#endif


if (my_pe==io_pe) then
   print *,'helicity: ',hetot
   print *,'helicity dissipation (spectrum): ',diss1*mu
   print *,'helicity dissipation (exact):    ',diss2*mu
endif

if (g_nz == 1)  then
   iwave = min(g_nx,g_ny)
else
   iwave = min(g_nx,g_ny,g_nz)
endif
iwave = (iwave/2)           ! max wave number in sphere.

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




subroutine compute_hfree_spec(Q,p1,work,skip_fft)
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
real*8 :: cmodes(2,nx,ny,nz,3)
real*8 :: energy,vx,wx,uy,wy,uz,vz,heltot
real*8 :: diss1,diss2,hetot,co_energy(3),xw,RR(3),II(3),mod_rr,mod_ii,&
      cos_tta, delta
integer i,j,k,jm,km,im,iwave_max,n

      if (ndim<3) then
         call abort("compute_helicity_specturm: can only be used in 3D")
      endif
      
      if (skip_fft==0) then
         p1(:,:,:,1:3)=Q(:,:,:,1:3)
         do n=1,3
         call fft3d(p1(1,1,1,n),work)
      enddo
      endif
      
      rwave=sqrt(  (g_nx/2.0)**2 + (g_ny/2.0)**2 + (g_nz/2.0)**2 )
      iwave_max=nint(rwave)
      
      do n = 1,3
         call sincos_to_complex(p1(1,1,1,n),cmodes(1,1,1,1,n),g_nmax)
      enddo
      

      hetot=0
      diss1=0
      diss2=0
      spec_helicity_rp=0
      spec_helicity_rn=0
      cos_tta_spec = 0

      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)
            do k=nz1,nz2
               km=kmcord(k)
               
               rwave = im**2 + jm**2 + km**2
               iwave = nint(sqrt(rwave))
               
!     compute angle between Re and Im parts of u(k)
               RR = cmodes(1,i,j,k,:)
               II = cmodes(2,i,j,k,:)
               
               mod_rr = sqrt(RR(1)*RR(1) + RR(2)*RR(2) + RR(3)*RR(3))
               mod_ii = sqrt(II(1)*II(1) + II(2)*II(2) + II(3)*II(3))
               
               cos_tta = (RR(1)*II(1) + RR(2)*II(2) + RR(3)*II(3))/&
               (mod_rr*mod_ii)
               
!     spectrum of angles
               cos_tta_spec(iwave) = cos_tta_spec(iwave) + cos_tta
               
!     cutoff for recalculating the spectra
               delta = 0.1      !this value can be changed by hand
               
!     omit modes where cos_tta is less than cutoff delta
               if (cos_tta > delta) then
            
!     ux = - pi2*im*p1(i+imsign(i),j,k,1)
                  vx = - pi2*im*p1(i+imsign(i),j,k,2)
                  wx = - pi2*im*p1(i+imsign(i),j,k,3)
                  
                  uy = - pi2*jm*p1(i,j+jmsign(j),k,1)
!     vy = - pi2*jm*p1(i,j+jmsign(j),k,2)
                  wy = - pi2*jm*p1(i,j+jmsign(j),k,3)
                  
                  uz =  - pi2*km*p1(i,j,k+kmsign(k),1)
                  vz =  - pi2*km*p1(i,j,k+kmsign(k),2)
!     wz =  - pi2*km*p1(i,j,k+kmsign(k),3)
                  
!     vorcity: ( (wy - vz), (uz - wx), (vx - uy) )

                  energy = 8
                  if (km==0) energy=energy/2
                  if (jm==0) energy=energy/2
                  if (im==0) energy=energy/2
                  
!     compute E(k) and kE(k)
                  xw=sqrt(rwave*pi2_squared)
                  spec_E(iwave)=spec_E(iwave) + 	energy* &
                  (p1(i,j,k,1)**2 + p1(i,j,k,2)**2 + p1(i,j,k,3)**2)
                  spec_kEk(iwave)=spec_kEk(iwave) + xw*energy* &
                  (p1(i,j,k,1)**2 + p1(i,j,k,2)**2 + p1(i,j,k,3)**2)
                  
                  energy = energy*(p1(i,j,k,1)*(wy-vz) + &
                  p1(i,j,k,2)*(uz-wx) + p1(i,j,k,3)*(vx-uy)) 
                  if (energy>0) spec_helicity_rp(iwave)= & 
                  spec_helicity_rp(iwave)+energy
                  if (energy<0) spec_helicity_rn(iwave)= &
                  spec_helicity_rn(iwave) + energy
                  
                  hetot=hetot+energy
                  diss1=diss1 -2*energy*iwave**2*pi2_squared
                  diss2=diss2 -2*energy*rwave*pi2_squared
               endif
            enddo
         enddo
      enddo


#ifdef USE_MPI
spec_r_in=spec_helicity_rp
call MPI_reduce(spec_r_in,spec_helicity_rp,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=spec_helicity_rn
call MPI_reduce(spec_r_in,spec_helicity_rn,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=spec_E
call MPI_reduce(spec_r_in,spec_E,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in=spec_kEk
call MPI_reduce(spec_r_in,spec_kEk,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_r_in = cos_tta_spec
call MPI_reduce(spec_r_in,cos_tta_spec,1+iwave_max,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)

do n=1,ndim
   spec_r_in=cospec_r(:,n)
   call MPI_reduce(spec_r_in,cospec_r(0,n),(1+iwave_max),MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
enddo

spec_x_in=cospec_x
call MPI_reduce(spec_x_in,cospec_x,(1+g_nx/2)*ndim,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_y_in=cospec_y
call MPI_reduce(spec_y_in,cospec_y,(1+g_ny/2)*ndim,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
spec_z_in=cospec_z
call MPI_reduce(spec_z_in,cospec_z,(1+g_nz/2)*ndim,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)


rwave=diss1
call MPI_reduce(rwave,diss1,1,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
rwave=diss2
call MPI_reduce(rwave,diss2,1,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
rwave=hetot
call MPI_reduce(rwave,hetot,1,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
#endif


if (my_pe==io_pe) then
   print *,'helicity: ',hetot
   print *,'helicity dissipation (spectrum): ',diss1*mu
   print *,'helicity dissipation (exact):    ',diss2*mu
endif

if (g_nz == 1)  then
   iwave = min(g_nx,g_ny)
else
   iwave = min(g_nx,g_ny,g_nz)
endif
iwave = (iwave/2)           ! max wave number in sphere.

! for all waves outside sphere, sum into one wave number:
do i=iwave+2,iwave_max
   spec_helicity_rp(iwave+1)=spec_helicity_rp(iwave+1)+spec_helicity_rp(i)
   spec_helicity_rn(iwave+1)=spec_helicity_rn(iwave+1)+spec_helicity_rn(i)
   spec_E(iwave+1)=spec_E(iwave+1)+spec_E(i)  
   spec_kEk(iwave+1)=spec_kEk(iwave+1)+spec_kEk(i)
   cos_tta_spec(iwave+1) = cos_tta_spec(iwave+1) + cos_tta_spec(i)   
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


end module
