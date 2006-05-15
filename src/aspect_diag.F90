#include "macros.h"
subroutine output_model(doit_model,doit_diag,time,Q,Qhat,q1,q2,q3,work1,work2)
use params
use spectrum
use isoave
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
integer :: nints_e
integer,parameter :: nints_e_max=100
real*8 :: ints_e(nints_e_max)
real*8 :: x,zero_len
real*8 :: divx,divi
integer i,j,k,n,ierr,csig
integer :: n1,n1d,n2,n2d,n3,n3d,pv_type,stype


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
   call compute_pv_spec(time,Q,q1,q2,q3,work1,work2)
!bw   call compute_bous
   call output_spec(time,time_initial)
!bw   call output_helicity_spec(time,time_initial)  ! put all hel spec in same file
!bw   call output_pv_spec(time,time_initial) ! Too complicated for now
   call output_pv2_spec(time,time_initial)
   call output_2d_spec(time,time_initial)  
!bw    call output_bous

   !set this flag so that for next timestep, we will compute and save
   !spectral transfer functions:
   compute_transfer=.false.


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
!bw compute scalars as a function of time, and output to a file
call compute_expensive_scalars(Q,Qhat,q1,q2,q3,work1,work2,nints_e,ints_e)
if (nints_e > nints_e_max) call abort("Error: aspec_diag.F90: nints_e_max too small")


! output turb scalars
if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".scalars-bous"
   call copen(message,"w",fid,ierr)
   write(6,*) "Opening -bous file"
   if (ierr/=0) then
      write(message,'(a,i5)') "diag_output(): Error opening .scalars-bous file errno=",ierr
      call abort(message)
   endif
   x=nints_e; call cwrite8(fid,x,1)
   call cwrite8(fid,time,1)
   call cwrite8(fid,ints_e,nints_e)
   call cclose(fid,ierr)
endif



!
! two-point correlations
!
if (diag_struct==1) then
   pv_type=1
   stype=3                   ! structure functions of u,v,w
   if (npassive==1) stype=4; ! structure functions of u,v,w and PV
   
   ! compute pv in work2, vorticity in q1
   call potential_vorticity(work1,q1,Q,q2,q3,pv_type)
   work2 = Q(:,:,:,4)  ! make a copy
   Q(:,:,:,4) = work1  ! overwrite 4'th component with PV
    
   ! angle averaged functions:
   call isoavep(Q,q1,q2,q3,stype,csig)
   ! if csig>0, isoavep did not complete - interrupted by SIGURG
   if (my_pe==io_pe .and. csig==0) then
      write(message,'(f10.4)') 10000.0000 + time
      message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".bisostr"
      call copen(message,"w",fid,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "output_model(): Error opening .isostr file errno=",ierr
         call abort(message)
      endif
      ! add the enstrophy dissipation (stored in position 8) 
      ! to the structure function file
      call writeisoave2(fid,time,ints_e(8),1)
      call cclose(fid,ierr)
   endif

   Q(:,:,:,4) = work2  ! restore
endif



end subroutine




subroutine compute_expensive_scalars(Q,Qhat,vor,potvor,theta_z,omegadotrho_nu,omegadotrho_kappa,nints_e,ints_e)
use params
use fft_interface
use mpi
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(g_nz2,nslabx,ny_2dz,n_var)
real*8 :: vor(nx,ny,nz,3)
real*8 :: potvor(nx,ny,nz)
real*8 :: theta_z(nx,ny,nz)
real*8 :: omegadotrho_nu(nx,ny,nz)
real*8 :: omegadotrho_kappa(nx,ny,nz)
integer :: nints_e
real*8  :: ints_e(*)

! local variables
real*8,allocatable :: ints2(:)
real*8 :: enstr_diss
real*8 :: dummy(1)
integer :: pv_type, i, j, k, im, jm, km, ierr
real*8 :: xw,xw2,u2,xw_viss,xfac, vx,wx,uy,wy,uz,vz,pe,totale,potens,pv
real*8 :: ke,pe_diss,ke_diss,h_diss,ens_diss,rwave,pv_diss,kappa


!
!  Compute scalars from the spectral data of prognostic variables 
!  stored in Qhat.  We need to add pe and pe_diss here.  
!
!  note: Qhat is stored in a more efficient FFT format, so
!  these loops do not look like the loops over FFT coefficients
!  that appear later. 
!
if (npassive==0) call abort("Error: aspect_diag.F90: not coded for npassive==0")
!if (mu_hyper_value/=0) call abort("Error: aspect_diag.F90: not coded for hyper viscosity")
if (mu_hypo_value/=0) call abort("Error: aspect_diag.F90: not coded for hypo viscosity")

! viscosity of theta:
kappa = mu/schmidt(np1)

pe = 0
ke = 0
ke_diss=0
pe_diss=0
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nslabx
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

            xw=(im*im + jm*jm + km*km/Lz/Lz)*pi2_squared
            xw_viss=mu*xw

            xfac = 2*2*2
            if (km==0) xfac=xfac/2
            if (jm==0) xfac=xfac/2
            if (im==0) xfac=xfac/2
            
            u2=Qhat(k,i,j,1)*Qhat(k,i,j,1) + &
                 Qhat(k,i,j,2)*Qhat(k,i,j,2) + &
                 Qhat(k,i,j,3)*Qhat(k,i,j,3)
            
            ke = ke + .5*xfac*u2
            pe = pe + .5*Qhat(k,i,j,np1)**2

            ke_diss = ke_diss + .5*xfac*xw_viss*u2
            pe_diss = pe_diss + &
                 .5*xfac*(xw*kappa)*Qhat(k,i,j,np1)**2

#if 0               
            vx = - pi2*im*Qhat(k,i+z_imsign(i),j,2)
            wx = - pi2*im*Qhat(k,i+z_imsign(i),j,3)
            uy = - pi2*jm*Qhat(k,i,j+z_jmsign(j),1)
            wy = - pi2*jm*Qhat(k,i,j+z_jmsign(j),3)
            uz =  - pi2*km*Qhat(k+z_kmsign(k),i,j,1)/Lz
            vz =  - pi2*km*Qhat(k+z_kmsign(k),i,j,2)/Lz
            ! vorcity: ( (wy - vz), (uz - wx), (vx - uy) )
            ! compute 2*k^2 u vor:
            h_diss = h_diss + 2*xfac*mu*xw*&
                 (Qhat(k,i,j,1)*(wy-vz) + &
                 Qhat(k,i,j,2)*(uz-wx) + &
                 Qhat(k,i,j,3)*(vx-uy)) 
            ens_diss = ens_diss + 2*xfac*mu*xw_viss*  &
                 ((wy-vz)**2 + (uz-wx)**2 + (vx-uy)**2) 
#endif            
      enddo
   enddo
enddo

ke = ke * Lz
pe = pe * Lz
ke = ke * Lz
pe = pe * Lz




!
!bw  Now computes the dissipation of Q = q^2/2
!bw  
!bw  enstr_diss = - nu < q laplacian (omega dot 
!bw                                   grad (rho_o - b z + \tilde{rho}))>
!bw  
!bw               - kappa < q laplacian (grad \tilde{rho}
!bw                                      dot  omega_a) >
!bw  
!bw  
!bw  
!bw  
!bw
!bw
!bw Compute this quantity in 3 stages in physical space
!bw   1. omega dot grad \tilde{rho}
!bw   2. omega_3 b
!bw   3. omega_3 f
!bw Then fft, then
!bw


! potvor = grad(Q(:,:,:,4)) dot vorticity  
! (use omegadotrho_* arrays as work arrays)
pv_type=2
call potential_vorticity(potvor,vor,Q,omegadotrho_nu,omegadotrho_kappa,pv_type)

! compute d/dz of theta,
call der(Q(1,1,1,np1),theta_z,dummy,omegadotrho_nu,DX_ONLY,3)

omegadotrho_nu = potvor - bous*theta_z/Lz
omegadotrho_kappa = potvor + fcor*theta_z/Lz

! now overwrite potvor with pvtype=1 potvor:
potvor = potvor + fcor*theta_z/Lz - bous*theta_z/Lz


!bw 
!bw Now laplacian both omegadotrho_nu  and omegadotrho_kappa
!bw
call fft3d(omegadotrho_nu,vor)
call fft3d(omegadotrho_kappa,vor)

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
!bw
!bw Fourier transform back to physical space because we need these quantites
!bw times the total potential vorticity integrated in the volume.
!bw
call fft3d(omegadotrho_nu,vor)
call fft3d(omegadotrho_kappa,vor)


pv = 0
pv_diss = 0
potens = 0
enstr_diss = 0
enstr_diss = 0
do k=nz1,nz2
   do j=ny1,ny2
      do i=nx1,nx2
         pv = pv + potvor(i,j,k)
         potens = potens + (potvor(i,j,k)*potvor(i,j,k)*.5)
         pv_diss = pv_diss + (mu*omegadotrho_nu(i,j,k) &
              +  kappa*omegadotrho_kappa(i,j,k)) 
         enstr_diss = enstr_diss + potvor(i,j,k)* &
              ( mu*omegadotrho_nu(i,j,k) +  kappa*omegadotrho_kappa(i,j,k))
         
      enddo
   enddo
enddo
! grid space calculations need to be normalized:
pv=pv/g_nx/g_ny/g_nz
pv_diss=pv_diss/g_nx/g_ny/g_nz
potens=potens/g_nx/g_ny/g_nz
enstr_diss=enstr_diss/g_nx/g_ny/g_nz


! store the integrals for output:
ints_e(1)=ke 
ints_e(2)=pe
ints_e(3)=ke_diss 
ints_e(4)=pe_diss
ints_e(5)=pv
ints_e(6)=pv_diss
ints_e(7)=potens
ints_e(8)=enstr_diss
nints_e=8


! global sum over all processors:
#ifdef USE_MPI
   allocate(ints2(nints_e))
   ints2(1:nints_e)=ints_e(1:nints_e)
   call mpi_allreduce(ints2,ints_e,nints_e,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   deallocate(ints2)
#endif



print *,'aspect_diag ke=',ke
print *,'aspect_diag pe=',pe



end subroutine compute_expensive_scalars

