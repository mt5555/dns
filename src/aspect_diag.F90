#include "macros.h"
subroutine output_model(doit_model,doit_diag,time,Q,Qhat,q1,q2,q3,work1,work2)
use params
use spectrum
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(g_nz2,nslabx,ny_2dz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: q2(nx,ny,nz,n_var)
real*8 :: q3(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time
logical :: doit_model,doit_diag

! local variables
integer,parameter :: nints_e=49,npints_e=51
real*8 :: ints_e(nints_e)
real*8 :: pints_e(npints_e,n_var)
real*8 :: x,zero_len
real*8 :: divx,divi
integer i,j,k,n,ierr,csig
integer :: n1,n1d,n2,n2d,n3,n3d
character(len=80) :: message
CPOINTER fid,fidj,fidS,fidcore




if (compute_transfer) then
   compute_transfer=.false.
   ! spec_r computed last time step
   ! spec_diff, spec_f, spec_rhs were computed in RHS computation at the
   ! beginning of this flag (becuase compute_transfer flag was set)
   ! So they are all known at time_old. Now compute spec_r_new 
   ! (used to compute edot_r)
   call compute_Edotspec(time,Q,q1,work1,work2)
   ! output all the spectrum:
   call output_tran(time,Q,q1,q2,q3,work1,work2)
endif


! compute spectrum
! always compute at first timestep because transfer cannot be computed
! on last timestep.   
if (doit_model .or. time==time_initial) then
if ( g_bdy_x1==PERIODIC .and. &
     g_bdy_y1==PERIODIC .and. &
     g_bdy_z1==PERIODIC) then


   call compute_spec(time,Q,q1,work1,work2)
   call compute_spec_2d(time,Q,q1,work1,work2)
   call output_spec(time,time_initial)
   call output_helicity_spec(time,time_initial)  ! put all hel spec in same file
   call output_2d_spec(time,time_initial)  

   !set this flag so that for next timestep, we will compute and save
   !spectral transfer functions:
   compute_transfer=.true.


   ! for incompressible equations, print divergence as diagnostic:
   if (equations==NS_UVW) then
      call compute_div(Q,q1,work1,work2,divx,divi)
      write(message,'(3(a,e12.5))') 'max(div)=',divx
      call print_message(message)	
   endif


endif
endif

! do PDF's and scalars if doit_model=.true., OR if this is a restart
! but we have computed new passive scalars.
if ((compute_passive_on_restart .and. time==time_initial) .or. &
    doit_model) then
   ! do the rest of this suburoutine
else
   return
endif



!
! the "expensive" scalars
!
! call compute_expensive_scalars(  )
! output these to a file...


end subroutine




subroutine compute_expensive_scalars
end subroutine




subroutine compute_pv_dissipation(q_diss,Q,Qhat,work1,work2,d1)
use params
use fft_interface
!bw
!bw  This routine computes the dissipation of Q = q^2/2
!bw  
!bw  enstr_diss = - nu < q laplacian (omega dot 
!bw                                   grad (rho_o - b z + \tilde{rho}) >
!bw  
!bw               - kappa < q laplacian (grad \tilde{rho}
!bw                                      dot  omega_a) >
!bw  
!bw  
!bw  
!bw  
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(g_nz2,nslabx,ny_2dz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: d1(nx,ny,nz)
real*8 :: vor(nx,ny,nz,3)
real*8 :: potvor(nx,ny,nz)
real*8 :: dummy(1)
real*8 :: enstr_diss, omegadotrho_nu(nx,ny,nz), omegadotrho_kappa(nx,ny,nz)
integer :: pv_type, i, j, k, im, jm, km


!
! For this subroutine pv_type = 1 (the full potential vorticity)
!

!bw
!bw I probably don't have the work arrays correct in these calls.
!bw
call vorticity(Q,vor,work1,work2)
call potential_vorticity(potvor,vor,Q,d1,work,pv_type)
!bw
!bw
!bw Compute this quantity in 3 stages in physical space
!bw   1. omega dot grad \tilde{rho}
!bw   2. omega_3 b
!bw   3. omega_3 f
!bw Then fft, then
!bw

omegadotrho_nu = 0

   call der(Q(1,1,1,4),d1,dummy,work1,DX_ONLY,1)
   do k=nz1,nz2
      do j=ny1,ny2
         do i=nx1,nx2
            omegadotrho_nu(i,j,k) = omegadotrho_nu(i,j,k) &
                              + d1(i,j,k)*vor(i,j,k,1)
         enddo
      enddo
   enddo
  ! compute theta_y
   call der(u(1,1,1,4),d1,dummy,work,DX_ONLY,2)
   do k=nz1,nz2
      do j=ny1,ny2
         do i=nx1,nx2
            omegadotrho_nu(i,j,k) = omegadotrho_nu(i,j,k) + d1(i,j,k)*vor(i,j,k,2)
         enddo
      enddo
   enddo
 
omegadotrho_kappa(i,j,k) = omegadotrho_nu(i,j,k)
  
   ! compute theta_z -- do not forget to add in coriolis here
   call der(u(1,1,1,4),d1,dummy,work,DX_ONLY,3)
   do k=nz1,nz2
      do j=ny1,ny2
         do i=nx1,nx2
            omegadotrho_nu(i,j,k) = omegadotrho_nu(i,j,k) + &
                            d1(i,j,k)/Lz*(vor(i,j,k,3)-bous)
            omegadotrho_kappa(i,j,k) = omegadotrho_kappa(i,j,k) + &
                                d1(i,j,k)/Lz*(vor(i,j,k,3)+fcor)
         enddo
      enddo
   enddo
!bw 
!bw Now laplacian both omegadotrho_nu  and omegadotrho_kappa
!bw
   call fft3d(omegadotrho_nu,work1)
   call fft3d(omegadotrho_kappa,work1)
!bw
   do k=nz1,nz2
      km=(kmcord(k))
      do j=ny1,ny2
         jm=(jmcord(j))
         do i=nx1,nx2
            im=(imcord(i))
            rwave = im**2 + jm**2 + (km/Lz)**2
            omegadotrho_nu(i,j,k) = -omegadotrho_nu(i,j,k)*rwave
            omegadotrho_kappa(i,j,k) = -omegadotrho_kappa(i,j,k)*rwave
         enddo
      enddo
   enddo
!bw Transform back -- I hope this is the right way to do this.
   call ifft3d(omegadotrho_nu,work1)
   call ifft3d(omegadotrho_kappa,work1)
!bw Now take the area average
   enstr_diss = 0
   do k=nz1,nz2
      km=(kmcord(k))
      do j=ny1,ny2
         jm=(jmcord(j))
         do i=nx1,nx2
            im=(imcord(i))
            xfac = 8
            if (km==0) xfac=xfac/2
            if (jm==0) xfac=xfac/2
            if (im==0) xfac=xfac/2
            enstr_diss = enstr_diss &
                         + nu*omegadotrho_nu(i,j,k) &
                         + kamma*omegadotrho_kappa(i,j,k)
         enddo
      enddo
   enddo


end subroutine







