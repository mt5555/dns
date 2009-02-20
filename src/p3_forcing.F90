!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Copyright 2007.  Los Alamos National Security, LLC. This material was
!produced under U.S. Government contract DE-AC52-06NA25396 for Los
!Alamos National Laboratory (LANL), which is operated by Los Alamos
!National Security, LLC for the U.S. Department of Energy. The
!U.S. Government has rights to use, reproduce, and distribute this
!software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
!LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
!FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
!derivative works, such modified software should be clearly marked, so
!as not to confuse it with the version available from LANL.
!
!Additionally, this program is free software; you can redistribute it
!and/or modify it under the terms of the GNU General Public License as
!published by the Free Software Foundation; either version 2 of the
!License, or (at your option) any later version. Accordingly, this
!program is distributed in the hope that it will be useful, but WITHOUT
!ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
!for more details.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "macros.h"

module p3_forcing
implicit none

type wnforcing_d
   integer :: n
   ! note: SGI f90 does not allow allocatable arrays in a struct
   integer, pointer :: index(:,:)
end type

integer           :: numb1,numb
integer,parameter :: numb_max=512
integer,parameter :: numbs=3         ! max stochastic forcing bands
real*8            :: FM(numbs)       ! stochastic forcing normalization
type(wnforcing_d),allocatable :: wnforcing(:)

integer,private :: ntot=0
integer :: init_sforcing=0
integer :: comm_sforcing  ! MPI communicator for all pe's involved in forcing

real*8 :: ener_target(numb_max)

contains




subroutine p3_sforce(rhs,Qhat,f_diss,fxx_diss)
!
! Add a forcing term to rhs, in spectral space, assuming P3DFFT decomposition.
!
use params
implicit none
complex*16 :: Qhat(p3_nx,p3_ny,p3_nz,3)
complex*16 :: rhs(p3_nx,p3_ny,p3_nz,3)
real*8 :: f_diss,param,fxx_diss

! determinisitic with E1=E2=.5
if (forcing_type==1) call p3_forcing12(rhs,Qhat,f_diss,fxx_diss,0)

end subroutine



subroutine p3_forcing12(rhs,Qhat,f_diss,fxx_diss,model_spec)
!
! Add a forcing term to rhs.
! Force 3D wave numbers 1 back to the sphere E=1**(-5/3)
! Force 3D wave numbers 2 back to the sphere E=2**(-5/3)
!
! model_spec==0    E(1)=E(2)=.5            'iso12' option
! model_spec==1    Overholt & Pope
! model_spec==2    Balu
! model_spec==3    high wave number for rotation case
! model_spec==4    E(1) = 1.0              'iso1' option
!
use params
use mpi
implicit none
complex*16 :: Qhat(p3_nx,p3_ny,p3_nz,3)
complex*16 :: rhs(p3_nx,p3_ny,p3_nz,3)
integer :: model_spec
integer km,jm,im,i,j,k,n,wn,ierr,kfmax
real*8 xw,xfac,f_diss,tauf,tau_inv,fxx_diss
real*8 ener(numb_max),temp(numb_max),Q2
character(len=80) :: message


if (0==init_sforcing) then
   numb1=1
   if (model_spec==0) then                        ! 'iso12' option
      numb=2   ! apply forcing in bands 1,2
      call sforcing_init()
      do wn=numb1,numb
         ener_target(wn)=.5
      enddo
   endif
   if (model_spec==1) then                         ! 'iso' option   Overhold&Pope
      numb=8
      call sforcing_init()
      do wn=numb1,numb
         ener_target(wn)=(real(wn)/numb)**4
      enddo
   endif
   if (model_spec==2) then     ! balu forcing      ' bal  
      numb1=dealias_23_kmax-1
      numb=dealias_23_kmax
      call sforcing_init()
      do wn=numb1,numb
         ener_target(wn)=wn**(-5./3.)
      enddo
   endif
   if (model_spec==3) then     ! For the rotation case  
      numb1=max(forcing_peak_waveno-8,1)
      numb=forcing_peak_waveno+8
      call sforcing_init()
      do wn=numb1,numb
         ener_target(wn)=exp(-.5*(wn-forcing_peak_waveno)**2)/sqrt(2*pi)
      enddo
   endif
   if (model_spec==4) then                        ! 'iso1' option
      numb=1   ! apply forcing in bands 1 only
      call sforcing_init()
      do wn=numb1,numb
         ener_target(wn)=1.0
      enddo
   endif
   if (numb>numb_max) call abortdns("sforcing12: numb_max too small")
   if (io_pe==my_pe) then
      print *,'Using Deterministic Low Wave Number Forcing. Target spectrum:'
      do wn=numb1,numb
         print *,'shell k=',wn,'  energy=',ener_target(wn)
      enddo
   endif
endif


if (g_u2xave==0) then
   ! on first call, this will not have been computed, so lets
   ! compute it here:  ints(2)
   g_u2xave=0
! note: we assume P3DFFT was compiled with STRIDE1 option.
! that means 1st dimension is really Z direction
!            2nd dimension is really X direction
!            3nd dimension is really Y direction
! X wave number im = p3_jmcord(j)
! Y wave number jm = p3_kmcord(k)
! Z wave number km = p3_imcord(i)
do k=1,p3_nz
   jm=p3_kmcord(k)
   do j=1,p3_ny
      im=p3_jmcord(j)
      do i=1,p3_nx
         km=p3_imcord(i)
         
            
            xw=(im*im + jm*jm + km*km)*pi2_squared
            
            xfac = 2
            if (im==0) xfac=xfac/2
            
            g_u2xave = g_u2xave + xfac*xw*&
                 real(Qhat(i,j,k,1)*conjg(Qhat(i,j,k,1)) + &
                 Qhat(i,j,k,2)*conjg(Qhat(i,j,k,2)) + &
                 Qhat(i,j,k,3)*conjg(Qhat(i,j,k,3)) )
            
         enddo
      enddo
   enddo
#ifdef USE_MPI
   xw=g_u2xave
   call mpi_allreduce(xw,g_u2xave,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif

endif


if (ntot==0) return
! only CPUS which belong to "comm_sforcing" beyond this point!


! relaxation coefficient (below) and kfmax (above) 
! comes from Overhold and Pope 1998.    .5 = tau/tau_kolmogorov
! tau = .5 tau_kolmogorov = .5 eta^2 / mu = .5 sqrt(mu/epsilon) = 
!                                           .5/sqrt( <ux,ux> ) 
!  
tau_inv=sqrt(g_u2xave)/.5



f_diss=0
fxx_diss=0
ener=0
do wn=numb1,numb

   do n=1,wnforcing(wn)%n
      i=wnforcing(wn)%index(n,1)
      j=wnforcing(wn)%index(n,2)
      k=wnforcing(wn)%index(n,3)
      xfac=2
      if (p3_imcord(i)==0) xfac=xfac/2
      ener(wn)=ener(wn)+.5*xfac*real( &
           Qhat(i,j,k,1)*conjg(Qhat(i,j,k,1)) +&
           Qhat(i,j,k,2)*conjg(Qhat(i,j,k,2)) +&
           Qhat(i,j,k,3)*conjg(Qhat(i,j,k,3))  )
   enddo
enddo
#ifdef USE_MPI
   temp=ener
   call mpi_allreduce(temp,ener,numb,MPI_REAL8,MPI_SUM,comm_sforcing,ierr)
#endif


do wn=numb1,numb
   ! Qf = Q*sqrt(ener_target/ener)
   ! forcing = 1/tau (Qf-Q) = 1/tau * (sqrt(ener_target/ener)-1) Q
   tauf=tau_inv*(sqrt(ener_target(wn)/ener(wn))-1)

   ! make sure relaxation is stable:
   if (delt*tauf>.5) tauf=.5/delt

!   if (io_pe==my_pe) &
!   write(*,'(a,3i4,2f17.10,e17.10)') 'FORCING:',my_pe,wn,wnforcing(wn)%n,ener(wn),ener_target(wn),tauf

   if (tauf>0) then ! only apply forcing if net input is positive
   do n=1,wnforcing(wn)%n
      i=wnforcing(wn)%index(n,1)
      j=wnforcing(wn)%index(n,2)
      k=wnforcing(wn)%index(n,3)
      rhs(i,j,k,1) = rhs(i,j,k,1) + tauf*Qhat(i,j,k,1)
      rhs(i,j,k,2) = rhs(i,j,k,2) + tauf*Qhat(i,j,k,2)
      rhs(i,j,k,3) = rhs(i,j,k,3) + tauf*Qhat(i,j,k,3)

      Q2 = real(Qhat(i,j,k,1)*conjg(Qhat(i,j,k,1)) + &
           Qhat(i,j,k,2)*conjg(Qhat(i,j,k,2)) + &
           Qhat(i,j,k,3)*conjg(Qhat(i,j,k,3)) )

      xfac=2
      if (p3_imcord(i)==0) xfac=xfac/2
      f_diss = f_diss + xfac*tauf*Q2

      xw=-(p3_imcord(i)**2 + p3_jmcord(j)**2 + p3_kmcord(k)**2)*pi2_squared
      fxx_diss = fxx_diss + xfac*tauf*xw*Q2

   enddo
   endif

enddo
end subroutine 
   
   







   
subroutine sforcing_init
!
! build data structure of modes with wave number 0<k<2.5
!
use params
use mpi
implicit none
real*8 :: xw
integer km,jm,im,i,j,k,n,wn,ierr
integer :: color,key
character(len=80) :: message

init_sforcing=1
allocate(wnforcing(numb1:numb))


   do n=numb1,numb
      wnforcing(n)%n=0
   enddo
   
   ! count the number of wavenumbers in each band on this CPU.  
! note: we assume P3DFFT was compiled with STRIDE1 option.
! that means 1st dimension is really Z direction
!            2nd dimension is really X direction
!            3nd dimension is really Y direction
! X wave number im = p3_jmcord(j)
! Y wave number jm = p3_kmcord(k)
! Z wave number km = p3_imcord(i)
   do k=1,p3_nz
      jm=p3_kmcord(k)
      do j=1,p3_ny
         im=p3_jmcord(j)
         do i=1,p3_nx
            km=p3_imcord(i)
            
            xw=sqrt(real(km**2+jm**2+im**2))
            do n=numb1,numb
               if (xw>=n-.5 .and. xw<n+.5) then
                  wnforcing(n)%n=wnforcing(n)%n+1
               endif
            enddo

         enddo
      enddo
   enddo
   
   ! allocate storage
   do n=numb1,numb
      i=wnforcing(n)%n
      if (i>0) allocate(wnforcing(n)%index(i,3))
      wnforcing(n)%n=0  ! reset counter to use again below
   enddo
   
   ! store all the indexes
   do k=1,p3_nz
      jm=p3_kmcord(k)
      do j=1,p3_ny
         im=p3_jmcord(j)
         do i=1,p3_nx
            km=p3_imcord(i)

            xw=sqrt(real(km**2+jm**2+im**2))
            do n=numb1,numb
            if (xw>=n-.5 .and. xw<n+.5) then
               wnforcing(n)%n=wnforcing(n)%n+1
               wnforcing(n)%index(wnforcing(n)%n,1)=i
               wnforcing(n)%index(wnforcing(n)%n,2)=j
               wnforcing(n)%index(wnforcing(n)%n,3)=k
            endif
            enddo

         enddo
      enddo
   enddo



ntot=0
do wn=numb1,numb
   ntot=ntot+wnforcing(wn)%n
enddo

color = 0
if (ntot>0) color=1
key=0

! everyone with ntot>0 joins a new group, comm_sforcing
#ifdef USE_MPI
call mpi_comm_split(comm_3d,color,key,comm_sforcing,ierr);
if (ntot==0) call mpi_comm_free(comm_sforcing,ierr)
#endif



end subroutine



end module
