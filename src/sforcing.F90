#include "macros.h"

module sforcing
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




subroutine sforce(rhs,Qhat,f_diss,fxx_diss)
!
! Add a forcing term to rhs, in spectral space.
!
use params
implicit none
real*8 :: Qhat(*)
real*8 :: rhs(*)
real*8 :: f_diss,param,fxx_diss

! determinisitic with E1=E2=.5
if (forcing_type==1) call sforcing12(rhs,Qhat,f_diss,fxx_diss,0)

! stochastic, wave numbers 1,2
if (forcing_type==2) call sforcing_random12(rhs,Qhat,f_diss,fxx_diss,0)

! determinisitic with Overholt and Pope spectrum
if (forcing_type==3) call sforcing12(rhs,Qhat,f_diss,fxx_diss,1)

! stochastic, wave numbers 2,3
if (forcing_type==4) call sforcing_random12(rhs,Qhat,f_diss,fxx_diss,0)

! determinisitic - Balu: single high wave number 
! do nothing - this forcing is handled in ns.F90
if (forcing_type==5) call sforcing12(rhs,Qhat,f_diss,fxx_diss,2)

! determinisitic with E1=E2=.5, also impose helicity
if (forcing_type==6) call sforcing12_helicity(rhs,Qhat,f_diss,fxx_diss,0)




end subroutine




subroutine gforce(Q,rhs,rhsz,q4,q4z,work,f_diss,fxx_diss)
!
! store forcing term to rhs, in gird space
!
use params
use fft_interface
implicit none
real*8 :: work(nx,ny,nz)
real*8 :: Q(nx,ny,nz,3)
real*8 :: q4(nx,ny,nz,3)
real*8 :: rhs(nx,ny,nz,3)
real*8 :: q4z(g_nz2,nslabx,ny_2dz,3) 
real*8 :: rhsz(g_nz2,nslabx,ny_2dz,3) 
real*8 :: f_diss,fxx_diss

! local
integer :: n
real*8 :: fdiss,fxxdiss,xfac,xw
integer im,jm,km,i,j,k


q4=Q
do n=1,3
   call z_fft3d_trashinput(q4(1,1,1,n),rhsz(1,1,1,n),work)
enddo
q4z=0

call sforce(q4z,rhsz,fdiss,fxxdiss)
do n=1,3
   call z_ifft3d(q4z(1,1,1,n),rhs(1,1,1,n),work)
enddo
f_diss=fdiss  ! strange bug on SGI, dnsghost, f_diss is non zero here, but
              ! return value in ns_ghost.F90 is 0.  if we introduce 'fdiss',
              ! everthing is ok
fxx_diss=fxxdiss
return
end subroutine






subroutine sforcing12(rhs,Qhat,f_diss,fxx_diss,model_spec)
!
! Add a forcing term to rhs.
! Force 3D wave numbers 1 back to the sphere E=1**(-5/3)
! Force 3D wave numbers 2 back to the sphere E=2**(-5/3)
!
! model_spec==0    E(1)=E(2)=.5
! model_spec==1    Overholt & Pope
! model_spec==2    Balu
!
use params
use mpi
implicit none
real*8 :: Qhat(g_nz2,nslabx,ny_2dz,3) 
real*8 :: rhs(g_nz2,nslabx,ny_2dz,3) 
integer :: model_spec
integer km,jm,im,i,j,k,n,wn,ierr,kfmax
real*8 xw,xfac,f_diss,tauf,tau_inv,fxx_diss
real*8 ener(numb_max),temp(numb_max)
character(len=80) :: message


if (0==init_sforcing) then
   numb1=1
   if (model_spec==0) then
      numb=2   ! apply forcing in bands 1,2
      call sforcing_init()
      do wn=numb1,numb
         ener_target(wn)=.5
      enddo
   endif
   if (model_spec==1) then
      numb=8
      call sforcing_init()
      do wn=numb1,numb
         ener_target(wn)=(real(wn)/numb)**4
      enddo
   endif
   ! balu forcing
   if (model_spec==2) then
      numb1=10
      numb=10
      call sforcing_init()
      do wn=numb1,numb
         ener_target(wn)=.063
      enddo
   endif
   if (numb>numb_max) call abort("sforcing12: numb_max too small")
endif


if (g_u2xave==0) then
   ! on first call, this will not have been computed, so lets
   ! compute it here:  ints(2)
   g_u2xave=0
   do j=1,ny_2dz
      jm=z_jmcord(j)
      do i=1,nslabx
         im=z_imcord(i)
         do k=1,g_nz
            km=z_kmcord(k)
            
            xw=(im*im + jm*jm + km*km)*pi2_squared
            
            xfac = 2*2*2
            if (km==0) xfac=xfac/2
            if (jm==0) xfac=xfac/2
            if (im==0) xfac=xfac/2
            
            
            g_u2xave = g_u2xave + xfac*xw*(Qhat(k,i,j,1)**2 + &
                 Qhat(k,i,j,2)**2 + &
                 Qhat(k,i,j,3)**2) 
            
         enddo
      enddo
   enddo
#ifdef USE_MPI
   xw=g_u2xave
   call MPI_allreduce(xw,g_u2xave,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
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
      xfac=8
      if (z_kmcord(k)==0) xfac=xfac/2
      if (z_jmcord(j)==0) xfac=xfac/2
      if (z_imcord(i)==0) xfac=xfac/2
      ener(wn)=ener(wn)+.5*xfac*(Qhat(k,i,j,1)**2+Qhat(k,i,j,2)**2+Qhat(k,i,j,3)**2)
   enddo
enddo
#ifdef USE_MPI
   temp=ener
   call MPI_allreduce(temp,ener,numb,MPI_REAL8,MPI_SUM,comm_sforcing,ierr)
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
      rhs(k,i,j,1) = rhs(k,i,j,1) + tauf*Qhat(k,i,j,1)
      rhs(k,i,j,2) = rhs(k,i,j,2) + tauf*Qhat(k,i,j,2)
      rhs(k,i,j,3) = rhs(k,i,j,3) + tauf*Qhat(k,i,j,3)

      xfac=8
      if (z_kmcord(k)==0) xfac=xfac/2
      if (z_jmcord(j)==0) xfac=xfac/2
      if (z_imcord(i)==0) xfac=xfac/2
      f_diss = f_diss + xfac*tauf*(Qhat(k,i,j,1)**2 + &
           Qhat(k,i,j,2)**2 + &
           Qhat(k,i,j,3)**2) 

      xw=-(z_imcord(i)**2 + z_jmcord(j)**2 + z_kmcord(k)**2)*pi2_squared
      fxx_diss = fxx_diss + xfac*tauf*xw*(Qhat(k,i,j,1)**2 + &
                                Qhat(k,i,j,2)**2 + &
                                Qhat(k,i,j,3)**2) 



   enddo
   endif

enddo
end subroutine 
   
   








subroutine sforcing12_helicity(rhs,Qhat,f_diss,fxx_diss,model_spec)
!
! Add a forcing term to rhs.
! Force 3D wave numbers 1 back to the sphere E=1**(-5/3)
! Force 3D wave numbers 2 back to the sphere E=2**(-5/3)
!
! and in addition, force to a prescribed angle to impose some helicity
!
! model_spec==0    E(1)=E(2)=.5
! model_spec==1    Overholt & Pope
!
use params
use mpi
implicit none
real*8 :: Qhat(g_nz2,nslabx,ny_2dz,3) 
real*8 :: rhs(g_nz2,nslabx,ny_2dz,3) 
integer :: model_spec
integer km,jm,im,i,j,k,k2,n,wn,ierr,kfmax
real*8 xw,xfac,f_diss,tauf,tau_inv,fxx_diss
real*8 ener(numb_max),temp(numb_max),Qdel(3)
real*8,allocatable,save :: rmodes(:,:,:,:)
real*8,allocatable,save :: rmodes2(:,:,:,:)
real*8,allocatable,save :: cmodes(:,:,:,:,:)
real*8 RR(3),II(3)

character(len=80) :: message


if (0==init_sforcing) then
   numb1=1
   if (model_spec==0) then
      numb=2   ! apply forcing in bands 1,2
      call sforcing_init()
      do wn=numb1,numb
         ener_target(wn)=.5
      enddo
   endif
   if (model_spec==1) then
      numb=8
      call sforcing_init()
      do wn=numb1,numb
         ener_target(wn)=(real(wn)/numb)**4
      enddo
   endif
   if (numb>numb_max) call abort("sforcing12_helicity: numb_max too small")
endif



if (g_u2xave==0) then
   ! on first call, this will not have been computed, so lets
   ! compute it here:  ints(2)
   g_u2xave=0
   do j=1,ny_2dz
      jm=z_jmcord(j)
      do i=1,nslabx
         im=z_imcord(i)
         do k=1,g_nz
            km=z_kmcord(k)
            
            xw=(im*im + jm*jm + km*km)*pi2_squared
            
            xfac = 2*2*2
            if (km==0) xfac=xfac/2
            if (jm==0) xfac=xfac/2
            if (im==0) xfac=xfac/2
            
            
            g_u2xave = g_u2xave + xfac*xw*(Qhat(k,i,j,1)**2 + &
                 Qhat(k,i,j,2)**2 + &
                 Qhat(k,i,j,3)**2) 
            
         enddo
      enddo
   enddo
#ifdef USE_MPI
   xw=g_u2xave
   call MPI_allreduce(xw,g_u2xave,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
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

if (.not. allocated(rmodes)) then
   allocate(rmodes(-numb:numb,-numb:numb,-numb:numb,3))
   allocate(rmodes2(-numb:numb,-numb:numb,-numb:numb,3))
   allocate(cmodes(2,-numb:numb,-numb:numb,-numb:numb,3))
endif


f_diss=0
fxx_diss=0
ener=0
rmodes=0
do wn=numb1,numb

   do n=1,wnforcing(wn)%n
      i=wnforcing(wn)%index(n,1)
      j=wnforcing(wn)%index(n,2)
      k=wnforcing(wn)%index(n,3)
      xfac=8
      if (z_kmcord(k)==0) xfac=xfac/2
      if (z_jmcord(j)==0) xfac=xfac/2
      if (z_imcord(i)==0) xfac=xfac/2
      ener(wn)=ener(wn)+.5*xfac*(Qhat(k,i,j,1)**2+Qhat(k,i,j,2)**2+Qhat(k,i,j,3)**2)

      rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),1)=Qhat(k,i,j,1)
      rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),2)=Qhat(k,i,j,2)
      rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),3)=Qhat(k,i,j,3)
   enddo
enddo
#ifdef USE_MPI
   temp=ener
   call MPI_allreduce(temp,ener,numb,MPI_REAL8,MPI_SUM,comm_sforcing,ierr)
   rmodes2=rmodes
   i=2*numb+1
   i=i*i*i
   call MPI_allreduce(rmodes2,rmodes,i,MPI_REAL8,MPI_SUM,comm_sforcing,ierr)
#endif


! Qf = Q*sqrt(ener_target/ener)
do wn=numb1,numb
   do n=1,wnforcing(wn)%n
      i=wnforcing(wn)%index(n,1)
      j=wnforcing(wn)%index(n,2)
      k=wnforcing(wn)%index(n,3)
      rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),:)= &
         sqrt(ener_target(wn)/ener(wn))*&
         rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),:)
   enddo
enddo


! convert to complex FFT coefficients (for testing)
!call random_number(rmodes); rmodes2=rmodes

do n=1,3
   ! note: using rmodes(:,:,:,n) fails under ifc/linux.  why?
   call sincos_to_complex(rmodes(-numb,-numb,-numb,n),&
      cmodes(1,-numb,-numb,-numb,n),numb)
enddo


! apply helicity fix:
do i=-numb,numb
do j=-numb,numb
do k=-numb,numb
   k2=i**2 + j**2 + k**2
   if (k2 < (.5+numb)**2 ) then
      RR = cmodes(1,i,j,k,:)
      II = cmodes(2,i,j,k,:)
   endif
enddo
enddo
enddo

! convert back:
do n=1,3
   call complex_to_sincos(rmodes(-numb,-numb,-numb,n),&
      cmodes(1,-numb,-numb,-numb,n),numb)
enddo

! for testing...
#if 0
do i=-numb,numb
do j=-numb,numb
do k=-numb,numb
   if (abs(rmodes(i,j,k,1)-rmodes2(i,j,k,1))>1e-10) then
   write(*,'(i3,i3,i3,f20.10,f20.10)') i,j,k,rmodes(i,j,k,1),rmodes2(i,j,k,1)
   endif
enddo
enddo
enddo
stop
#endif


do wn=numb1,numb
   ! Qf = Q*sqrt(ener_target/ener)
   ! forcing = 1/tau (Qf-Q) = 1/tau * (sqrt(ener_target/ener)-1) Q

   tauf=tau_inv
   ! make sure tauf is not too large that forcing is unstable:
   if (delt*tau_inv*(sqrt(ener_target(wn)/ener(wn))-1) > .5) then
         tauf = .5/(delt*(sqrt(ener_target(wn)/ener(wn))-1))
   endif


   if (ener(wn)<ener_target(wn)) then ! only apply forcing if net input is positive
   do n=1,wnforcing(wn)%n
      i=wnforcing(wn)%index(n,1)
      j=wnforcing(wn)%index(n,2)
      k=wnforcing(wn)%index(n,3)

      Qdel(1)=tauf*(rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),1)-Qhat(k,i,j,1))
      Qdel(2)=tauf*(rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),2)-Qhat(k,i,j,2))
      Qdel(3)=tauf*(rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),3)-Qhat(k,i,j,3))


      rhs(k,i,j,1) = rhs(k,i,j,1) + Qdel(1)
      rhs(k,i,j,2) = rhs(k,i,j,2) + Qdel(2)
      rhs(k,i,j,3) = rhs(k,i,j,3) + Qdel(3)

      xfac=8
      if (z_kmcord(k)==0) xfac=xfac/2
      if (z_jmcord(j)==0) xfac=xfac/2
      if (z_imcord(i)==0) xfac=xfac/2
      f_diss = f_diss + xfac*(Qdel(1)*Qhat(k,i,j,1) +&
                              Qdel(2)*Qhat(k,i,j,2) + &
                              Qdel(3)*Qhat(k,i,j,3) ) 

      xw=-(z_imcord(i)**2 + z_jmcord(j)**2 + z_kmcord(k)**2)*pi2_squared
      fxx_diss = fxx_diss + xfac*xw*(Qdel(1)**2 + Qdel(2)**2 + Qdel(3)**2) 

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
   do j=1,ny_2dz
      jm=z_jmcord(j)
      do i=1,nslabx
         im=z_imcord(i)
         do k=1,g_nz
            km=z_kmcord(k)

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
   do j=1,ny_2dz
      jm=z_jmcord(j)
      do i=1,nslabx
         im=z_imcord(i)
         do k=1,g_nz
            km=z_kmcord(k)

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
call MPI_Comm_split(comm_3d,color,key,comm_sforcing,ierr);
if (ntot==0) call MPI_Comm_free(comm_sforcing,ierr)
#endif



end subroutine








subroutine sforcing_random12(rhs,Qhat,f_diss,fxx_diss,new_f)
!
! Add a forcing term to rhs.
! Random, isotropic, homogenious in first 2 wave nubmers
!
! if new_f=0, add previously computed forcing stored in rmodes() into RHS. 
! if new_f=1, compute a new forcing ONLY, store in rmodes()
!
! numb1,numb2:   force bands [numb1,numb2]
!
use params
use mpi
implicit none
integer :: new_f
real*8 :: Qhat(g_nz2,nslabx,ny_2dz,3) 
real*8 :: rhs(g_nz2,nslabx,ny_2dz,3) 
integer km,jm,im,i,j,k,n,wn,ierr
real*8 xw,xfac,f_diss,fsum,fxx_diss

real*8,save :: rmodes(-numbs:numbs,-numbs:numbs,-numbs:numbs,3)       ! value at time tmod
real*8,save :: tmod,tmod_old
real*8,save :: tscale=.01


if (0==init_sforcing) then
   if (forcing_type==2) then
      numb1=1
      numb=2
      FM(1) = 18/(2*pi)**2    ! 15,15,0 yields: mu <ux,ux> = .5
      FM(2) = 18/(2*pi)**2
      FM(3) = 0
   endif
   if (forcing_type==4) then
      numb1=2
      numb=3
      FM(1) = 0    ! choose these so that mu <ux,ux> = .5
      FM(2) = 30    ! 0,5,5: epsilon: .09     R_l=35  E=.08
      FM(3) = 30    ! 0,15,15:        .27     R_l=45  E=.19
                    ! 0,30,30:        .56     R_l=49  E=.30
   endif

   call sforcing_init()
   rmodes=0
   tmod=0

   ! check that we are not including any wave numbers > numbs
   do wn=numb1,numb
   do n=1,wnforcing(wn)%n
      i=wnforcing(wn)%index(n,1)
      j=wnforcing(wn)%index(n,2)
      k=wnforcing(wn)%index(n,3)
      if (abs(z_imcord(i))>numbs .or. abs(z_jmcord(j))>numbs .or. abs(z_kmcord(k))>numbs)  then 
           call abort("Index error in sforcing_random12")
      endif
   enddo
   enddo
endif

if (ntot==0) return
!
! only CPUS which belong to "comm_sforcing" beyond this point!
!

!
! rmodes is at time  tmod
!

!
! Compute a new forcing function?  
!
if (new_f==1) then
   if (my_pe==io_pe) call random12(rmodes)
#ifdef USE_MPI
   call MPI_bcast(rmodes,3*(2*numbs+1)**3,MPI_REAL8,io_pe,comm_sforcing ,ierr)
#endif
   return
endif



f_diss=0
fxx_diss=0
do wn=numb1,numb
   do n=1,wnforcing(wn)%n
      i=wnforcing(wn)%index(n,1)
      j=wnforcing(wn)%index(n,2)
      k=wnforcing(wn)%index(n,3)

      rhs(k,i,j,1)=rhs(k,i,j,1) + rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),1)
      rhs(k,i,j,2)=rhs(k,i,j,2) + rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),2)
      rhs(k,i,j,3)=rhs(k,i,j,3) + rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),3)

      xfac=8
      if (z_kmcord(k)==0) xfac=xfac/2
      if (z_jmcord(j)==0) xfac=xfac/2
      if (z_imcord(i)==0) xfac=xfac/2
      fsum= ( &
         Qhat(k,i,j,1)*rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),1) +&
         Qhat(k,i,j,2)*rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),2) +&
         Qhat(k,i,j,3)*rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),3) )
      
      f_diss = f_diss + xfac*fsum
      xw=-(z_imcord(i)**2 + z_jmcord(j)**2 + z_kmcord(k)**2)*pi2_squared
      fxx_diss = fxx_diss + xfac*xw*fsum


   enddo

enddo



end subroutine 
   
   





subroutine random12(rmodes)
use params
implicit none
real*8 :: rmodes(-numbs:numbs,-numbs:numbs,-numbs:numbs,3)       

integer km,jm,im,i,j,k,n,wn,ierr,k2,count,countmax
real*8 xw,xfac,f_diss
real*8 :: R(-numbs:numbs,-numbs:numbs,-numbs:numbs,3,2),Rr,Ri
real*8 :: psix_r(3),psix_i(3)
real*8 :: psiy_r(3),psiy_i(3)
real*8 :: psiz_r(3),psiz_i(3)
real*8 :: ff(numbs,3)

#undef GATHER_STATS
#ifdef GATHER_STATS
ff=0
Ri=3
countmax=1000000
do count=1,countmax

call gaussian(rmodes,((2*numbs+1)**3) *3)
do wn=numb1,numb

   do n=1,wnforcing(wn)%n
      i=wnforcing(wn)%index(n,1)
      j=wnforcing(wn)%index(n,2)
      k=wnforcing(wn)%index(n,3)

      xfac=8
      if (z_kmcord(k)==0) xfac=xfac/2
      if (z_jmcord(j)==0) xfac=xfac/2
      if (z_imcord(i)==0) xfac=xfac/2


      ff(wn,1)= ff(wn,1) + xfac*( &
         rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),1)**Ri )
      ff(wn,2)= ff(wn,2) + xfac*( &
         rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),2)**Ri )
      ff(wn,3)= ff(wn,3) + xfac*( &
         rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),3)**Ri )


   enddo

enddo
enddo
ff=abs(ff/countmax)**(1/Ri)
print *,'count=',countmax
print *,'<R,R> ='
do wn=numb1,numb
   print *,ff(wn,:)
enddo
call abort("end stats")
#endif




rmodes=0
call gaussian(R,((2*numbs+1)**3) *3*2)





do km=-numb,numb 
do jm=-numb,numb
do im=-numb,numb
   !
   ! Choose R gaussian, theta uniform from [0..1]
   ! vorticty = (R1 + i R2)  exp(im*2pi*x) * exp(jm*2pi*y) * exp(km*2pi*z) 
   !
   ! R = vorticity
   ! PSI = laplacian^-1 R
   ! f = curl R 
   !
   ! code above (#define COMPUTE_STATS) shows that
   !  k=1  <R,R> = 180
   !  k=2  <R,R> = 1090
   !  k=3  <R,R> = 1809
   !
   k2 = (im**2 + jm**2 + km**2)
   ! ignore wave numbers outside of our forcing band 
   if ((numb1-.5)**2 <=k2 .and.   k2 < (numb+.5)**2) then

      if (k2>0) then
         xfac=1/(-2*pi*k2)
      else
         xfac=0
      endif
      if (delt>0) xfac=xfac/sqrt(delt)
      ! normalize so that < R**2 >   = F1  in wave number 1
      ! normalize so that < R**2 >   = F2  in wave number 2
      if (k2 < 1.5**2) then
         xfac=xfac*sqrt(FM(1)/180)
      else if (k2 < 2.5**2)  then
         xfac=xfac*sqrt(FM(2)/1090)
      else if (k2 < 3.5**2)  then
         xfac=xfac*sqrt(FM(3)/1809)
      endif
      


      do n=1,3
         psix_r(n) = -im*R(im,jm,km,n,2)*xfac
         psix_i(n) =  im*R(im,jm,km,n,1)*xfac
         psiy_r(n) = -jm*R(im,jm,km,n,2)*xfac
         psiy_i(n) =  jm*R(im,jm,km,n,1)*xfac
         psiz_r(n) = -km*R(im,jm,km,n,2)*xfac
         psiz_i(n) =  km*R(im,jm,km,n,1)*xfac
      enddo


      ! vor(1) = w_y - v_z
      ! vor(2) = u_z - w_x 
      ! vor(3) = v_x - u_y
      do n=1,3
         if (n==1) then
            Rr = psiy_r(3) - psiz_r(2)
            Ri = psiy_i(3) - psiz_i(2)
         else if (n==2) then
            Rr = psiz_r(1) - psix_r(3)
            Ri = psiz_i(1) - psix_i(3)
         else if (n==3) then
            Rr = psix_r(2) - psiy_r(1)
            Ri = psix_i(2) - psiy_i(1)
         endif

         
         i=abs(im)
         j=abs(jm)
         k=abs(km)
         ! code taken from complex_to_sincos().  Any bug fixes here,
         ! be sure to apply to that function also.
         rmodes( i, j, k,n) = rmodes( i, j, k,n) + Rr 
         rmodes( i, j,-k,n) = rmodes( i, j,-k,n) - Ri*zerosign(km) 
         rmodes( i,-j, k,n) = rmodes( i,-j, k,n) - Ri*zerosign(jm) 
         rmodes( i,-j,-k,n) = rmodes( i,-j,-k,n) - Rr*zerosign(jm*km) 
         rmodes(-i, j, k,n) = rmodes(-i, j, k,n) - Ri*zerosign(im)
         rmodes(-i, j,-k,n) = rmodes(-i, j,-k,n) - Rr*zerosign(im*km) 
         rmodes(-i,-j, k,n) = rmodes(-i,-j, k,n) - Rr*zerosign(im*jm) 
         rmodes(-i,-j,-k,n) = rmodes(-i,-j,-k,n) + Ri*zerosign(im*jm*km) 
      enddo
   endif
enddo
enddo
enddo
end subroutine





integer function zerosign(i)
integer i
if (i==0) then
   zerosign=0
else if (i<0) then
   zerosign=-1
else
   zerosign=1
endif


end function





subroutine zdecomp_to_rmodes(p,rmodes,nmax)
!
! extract sine and cosine modes from p, store in rmodes
!
use params
use mpi
implicit none
real*8 :: rmodes(-nmax:nmax,-nmax:nmax,-nmax:nmax)
real*8 :: rmodes2(-nmax:nmax,-nmax:nmax,-nmax:nmax)
real*8 :: p(g_nz2,nslabx,ny_2dz)
integer :: nmax,i,j,k,ierr

rmodes=0
do j=1,ny_2dz
do i=1,nslabx
do k=1,g_nz
   if (abs(z_imcord(i))<=nmax .and. abs(z_jmcord(j))<=nmax .and. &
       abs(z_kmcord(k))<=nmax ) then
      rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k)) = p(k,i,j) 
   endif
enddo
enddo
enddo
#ifdef USE_MPI
i=2*nmax+1
i=i*i*i
rmodes2=rmodes
call MPI_allreduce(rmodes2,rmodes,i,MPI_REAL8,MPI_SUM,comm_3d,ierr)

#endif

end subroutine




subroutine rmodes_to_zdecomp(p,rmodes,nmax)
!
! extract sine and cosine modes from rmodes, store in p
!
use params
implicit none
real*8 :: rmodes(-nmax:nmax,-nmax:nmax,-nmax:nmax)
real*8 :: p(g_nz2,nslabx,ny_2dz)
integer :: nmax,i,j,k

do j=1,ny_2dz
do i=1,nslabx
do k=1,g_nz
   if (abs(z_imcord(i))<=nmax .and. abs(z_jmcord(j))<=nmax .and. &
       abs(z_kmcord(k))<=nmax ) then
       p(k,i,j) = rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k)) 
   endif
enddo
enddo
enddo



end subroutine






subroutine complex_to_sincos(rmodes,cmodes,nmax)
!
!  convert set of complex FFT coefficients to sine and cosines
!
use params
implicit none
real*8 :: Rr,Ri
real*8 :: cmodes(2,-nmax:nmax,-nmax:nmax,-nmax:nmax)
real*8 :: rmodes(-nmax:nmax,-nmax:nmax,-nmax:nmax)
integer :: i,j,k,im,jm,km,nmax,imax,ip,jp,kp
!
! Note: this code also used in random12() above. 
! Any bugfixes here, also apply to random12().
!   
! convert to sine & cosine modes:
!
! (R1 + i R2) (cosx + i sinx)  (cosy + i siny)  (cosz + i sinz)  
! = (real parts only:)
!   R1  cosx cosy cosz                      R1 (1,1,1)
! i R2  cosx cosy sinz  i                  -R2 (1,1,-1) zerosign(km) 
! i R2  cosx siny cosz  i                  -R2 (1,-1,1) zerosign(jm)
!   R1  cosx siny sinz  i**2               -R1 (1,-1,-1) zerosign(km*jm)
! i R2  sinx cosy cosz  i                  -R2 (-1,1,1) zerosign(im) 
!   R1  sinx cosy sinz  i**2               -R1 (-1,1,-1) zerosign(im*km)
!   R1  sinx siny cosz  i**2               -R1 (-1,-1,1) zerosign(im*jm)
! i R2  sinx siny sinz  i**3                R2 (-1,-1,-1) zerosign(im*jm*km)
!  
! 

rmodes=0

do im=-nmax,nmax
do jm=-nmax,nmax
do km=-nmax,nmax
   Rr = cmodes(1,im,jm,km)
   Ri = cmodes(2,im,jm,km)

   i=abs(im)
   j=abs(jm)
   k=abs(km)
   
   rmodes( i, j, k) = rmodes( i, j, k) + Rr 
   rmodes( i, j,-k) = rmodes( i, j,-k) - Ri*zerosign(km) 
   rmodes( i,-j, k) = rmodes( i,-j, k) - Ri*zerosign(jm) 
   rmodes( i,-j,-k) = rmodes( i,-j,-k) - Rr*zerosign(jm*km) 
   rmodes(-i, j, k) = rmodes(-i, j, k) - Ri*zerosign(im)
   rmodes(-i, j,-k) = rmodes(-i, j,-k) - Rr*zerosign(im*km) 
   rmodes(-i,-j, k) = rmodes(-i,-j, k) - Rr*zerosign(im*jm) 
   rmodes(-i,-j,-k) = rmodes(-i,-j,-k) + Ri*zerosign(im*jm*km) 

enddo
enddo
enddo



end subroutine








subroutine sincos_to_complex(rmodes,cmodes,nmax)
!
!  convert set of sine and cosine FFT coefficients to complex coefficients
!
#if 0
   conversion to complex coefficients:
        a cos(lx) cos(my) cos(nz)
      
        sign(>=0)=1
        sign(<0)=-1

        p = number of negative values in sign(l), sign(m), sign(n):   
        1/i**0    1
        1/i**1   -i
        1/i**2   -1
        1/i**3    i


        a/(8 i**p) (exp(ilx)+sign(l)exp(-ilx))
                   (exp(imy)+sign(m)exp(-imy))
                   (exp(inz)+sign(n)exp(-inz))

        ( l, m, n)
        ( l, m,-n) sign(n)
        ( l,-m, n) sign(m)
        ( l,-m,-n) sign(n)*sign(m)
        (-l, m, n) sign(l)
        (-l, m,-n) sign(l)*sign(n)
        (-l,-m, n) sign(l)*sign(m)
        (-l,-m,-n) sign(l)*sign(m)*sign(n)

#endif      
use params
implicit none
real*8 :: rmodes(-nmax:nmax,-nmax:nmax,-nmax:nmax)
real*8 :: cmodes(2,-nmax:nmax,-nmax:nmax,-nmax:nmax)
real*8 :: a,b
integer :: i,j,k,im,jm,km,imax,nmax,sm,ip,jp,kp

imax=2*nmax+2
cmodes=0

do im=-nmax,nmax
do jm=-nmax,nmax
do km=-nmax,nmax

   ip=abs(im)
   jp=abs(jm)
   kp=abs(km)


   a=0; b=0

   ! count the number if sin() terms:
   sm=0; if (im<0) sm=sm+1;  if (jm<0) sm=sm+1;  if (km<0) sm=sm+1

   if (sm==0) then
      a=rmodes(im,jm,km)/8
   else if (sm==1) then
      b=-rmodes(im,jm,km)/8
   else if (sm==2) then
      a=-rmodes(im,jm,km)/8
   else if (sm==3) then
      b=rmodes(im,jm,km)/8
   else
      call abort("this cant happen")
   endif

   cmodes(1,ip,jp,kp)=cmodes(1,ip,jp,kp) + a;    
   cmodes(2,ip,jp,kp)=cmodes(2,ip,jp,kp) + b

   cmodes(1,ip,jp,-kp)=cmodes(1,ip,jp,-kp) + a*sign(1,km)   
   cmodes(2,ip,jp,-kp)=cmodes(2,ip,jp,-kp) + b*sign(1,km)

   cmodes(1,ip,-jp,kp)=cmodes(1,ip,-jp,kp) + a*sign(1,jm)
   cmodes(2,ip,-jp,kp)=cmodes(2,ip,-jp,kp) + b*sign(1,jm)

   cmodes(1,ip,-jp,-kp)=cmodes(1,ip,-jp,-kp) + a*sign(1,jm)*sign(1,km)  
   cmodes(2,ip,-jp,-kp)=cmodes(2,ip,-jp,-kp) + b*sign(1,jm)*sign(1,km)  

   cmodes(1,-ip,jp,kp)=cmodes(1,-ip,jp,kp) + a*sign(1,im)
   cmodes(2,-ip,jp,kp)=cmodes(2,-ip,jp,kp) + b*sign(1,im)

   cmodes(1,-ip,jp,-kp)=cmodes(1,-ip,jp,-kp) + a*sign(1,im)*sign(1,km)
   cmodes(2,-ip,jp,-kp)=cmodes(2,-ip,jp,-kp) + b*sign(1,im)*sign(1,km)

   cmodes(1,-ip,-jp,kp)=cmodes(1,-ip,-jp,kp) + a*sign(1,im)*sign(1,jm)
   cmodes(2,-ip,-jp,kp)=cmodes(2,-ip,-jp,kp) + b*sign(1,im)*sign(1,jm)

   cmodes(1,-ip,-jp,-kp)=cmodes(1,-ip,-jp,-kp) + a*sign(1,im)*sign(1,jm)*sign(1,km)
   cmodes(2,-ip,-jp,-kp)=cmodes(2,-ip,-jp,-kp) + b*sign(1,im)*sign(1,jm)*sign(1,km)

enddo
enddo
enddo
end subroutine





end module
