#include "macros.h"

module sforcing
implicit none

type wnforcing_d
   integer :: n
   ! note: SGI f90 does not allow allocatable arrays in a struct
   integer, pointer :: index(:,:)
end type
#define NUMBANDS 2
type(wnforcing_d) :: wnforcing(NUMBANDS)

integer,private :: ntot=0
integer :: init_sforcing=0
integer :: comm_sforcing  ! MPI communicator for all pe's involved in forcing

real*8 :: tau
contains


subroutine sforce(rhs,Qhat,f_diss)
!
! Add a forcing term to rhs, in spectral space.
!
use params
implicit none
real*8 :: Qhat(*)
real*8 :: rhs(*)
real*8 :: f_diss
if (forcing_type==1) call sforcing12(rhs,Qhat,f_diss)

end subroutine


subroutine gforce(Q,rhs,rhsz,q4,q4z,work,f_diss)
!
! store forcing term to rhs, in gird space
!
use params
use fft_interface
implicit none
real*8 :: work(nx,ny,nz)
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: q4(nx,ny,nz,n_var)
real*8 :: rhs(nx,ny,nz,n_var)
real*8 :: q4z(g_nz2,nslabx,ny_2dz,n_var) 
real*8 :: rhsz(g_nz2,nslabx,ny_2dz,n_var) 
real*8 :: f_diss
integer :: n

q4=Q
do n=1,3
   call z_fft3d_trashinput(q4(1,1,1,n),rhsz(1,1,1,n),work)
enddo
q4z=0
call sforce(q4z,rhsz,f_diss)
do n=1,3
   call z_ifft3d(q4z(1,1,1,n),rhs(1,1,1,n),work)
enddo
end subroutine






subroutine sforcing12(rhs,Qhat,f_diss)
!
! Add a forcing term to rhs.
! Force 3D wave numbers 1 back to the sphere E=1**(-5/3)
! Force 3D wave numbers 2 back to the sphere E=2**(-5/3)
!
use params
use mpi
implicit none
real*8 :: Qhat(g_nz2,nslabx,ny_2dz,n_var) 
real*8 :: rhs(g_nz2,nslabx,ny_2dz,n_var) 
integer km,jm,im,i,j,k,n,wn,ierr
real*8 xw,xfac,f_diss,tauf
real*8 ener(NUMBANDS),ener_target(NUMBANDS),temp(NUMBANDS)


if (0==init_sforcing) then
   call sforcing_init()
endif
if (ntot==0) return

! only CPUS which belong to "comm_sforcing" beyond this point!

f_diss=0
do wn=1,NUMBANDS
   ener_target(wn)=wn**(-5.0/3.0)
   if (wn>2) ener_target(wn)=wn**(-7.0/3.0)

   ener(wn)=0
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
   call MPI_allreduce(temp,ener,NUMBANDS,MPI_REAL8,MPI_SUM,comm_sforcing,ierr)
!   call MPI_allreduce(temp,ener,NUMBANDS,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif


do wn=1,NUMBANDS
   ! Qf = Q*sqrt(ener_target/ener)
   ! forcing = tau (Qf-Q) = tau * (sqrt(ener_target/ener)-1) Q
   tauf=tau*(sqrt(ener_target(wn)/ener(wn))-1)
!   print *,'FORCING:',wn,ener(wn),ener_target(wn)
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

   enddo
enddo
end subroutine 
   
   













   
subroutine sforcing_init
use params
use mpi
implicit none
real*8 :: xw
integer km,jm,im,i,j,k,n,wn,ierr
integer :: color,key
character(len=80) :: message

init_sforcing=1

   if (init_cond_subtype==1) then
      tau=25
   else if (init_cond_subtype==2) then
      tau=1
   else if (init_cond_subtype==3) then
      tau=500
   else
      tau=5
   endif
   write(message,'(a,f8.2)') 'Forcing relaxation parameter tau=',tau
   call print_message(message)

   do n=1,NUMBANDS
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
            do n=1,NUMBANDS
               if (xw>=n-.5 .and. xw<n+.5) then
                  wnforcing(n)%n=wnforcing(n)%n+1
               endif
            enddo

         enddo
      enddo
   enddo
   
   ! allocate storage
   do n=1,NUMBANDS
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
            do n=1,NUMBANDS
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
do wn=1,NUMBANDS
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



end module
