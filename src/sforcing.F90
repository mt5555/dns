subroutine sforcing(rhs,Qhat,f_diss)
!
! Add a forcing term to rhs, in spectral space.
!
use params
implicit none
real*8 :: Qhat
real*8 :: rhs
real*8 :: f_diss
if (forcing_type==1) call sforcing12(rhs,Qhat,f_diss)

end


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
real*8 xw,tau,xfac,f_diss
logical,save :: firstcall=.true.

type wnforcing_d
   integer :: n
   integer, allocatable :: index(:,:)
end type
#define NUMBANDS 2
type(wnforcing_d),save :: wnforcing(NUMBANDS)
real*8 ener(NUMBANDS),ener_target(NUMBANDS),temp(NUMBANDS)


if (firstcall) then
   firstcall=.false.

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
            if (xw>=.5 .and. xw<1.5) then
               wnforcing(1)%n=wnforcing(1)%n+1
            endif
            if (xw>=1.5 .and. xw<2.5) then
               wnforcing(2)%n=wnforcing(2)%n+1
            endif
         enddo
      enddo
   enddo
   
   ! allocate storage
   do n=1,NUMBANDS
      i=wnforcing(n)%n
      allocate(wnforcing(n)%index(i,3))
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
            if (xw>=.5 .and. xw<1.5) then
               wnforcing(1)%n=wnforcing(1)%n+1
               wnforcing(1)%index(wnforcing(1)%n,1)=i
               wnforcing(1)%index(wnforcing(1)%n,2)=j
               wnforcing(1)%index(wnforcing(1)%n,3)=k
            endif
            
            if (xw>=1.5 .and. xw<2.5) then
               wnforcing(2)%n=wnforcing(2)%n+1
               wnforcing(2)%index(wnforcing(2)%n,1)=i
               wnforcing(2)%index(wnforcing(2)%n,2)=j
               wnforcing(2)%index(wnforcing(2)%n,3)=k
            endif
         enddo
      enddo
   enddo
endif



f_diss=0
do wn=1,NUMBANDS
   ener_target(wn)=wn**(-5.0/3.0)
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
#ifdef USE_MPI
   temp=ener
   call MPI_allreduce(temp,ener,NUMBANDS,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
enddo


do wn=1,NUMBANDS
   ! Qf = Q*sqrt(ener_target/ener)
   ! forcing = tau (Qf-Q) = tau * (sqrt(ener_target/ener)-1) Q
   tau=1.0*(sqrt(ener_target(wn)/ener(wn))-1)
!   print *,'FORCING:',wn,ener(wn),ener_target(wn)
   do n=1,wnforcing(wn)%n
      i=wnforcing(wn)%index(n,1)
      j=wnforcing(wn)%index(n,2)
      k=wnforcing(wn)%index(n,3)
      rhs(k,i,j,1) = rhs(k,i,j,1) + tau*Qhat(k,i,j,1)
      rhs(k,i,j,2) = rhs(k,i,j,2) + tau*Qhat(k,i,j,2)
      rhs(k,i,j,3) = rhs(k,i,j,3) + tau*Qhat(k,i,j,3)

      xfac=8
      if (z_kmcord(k)==0) xfac=xfac/2
      if (z_jmcord(j)==0) xfac=xfac/2
      if (z_imcord(i)==0) xfac=xfac/2
      f_diss = f_diss + xfac*tau*(Qhat(k,i,j,1)**2 + &
           Qhat(k,i,j,2)**2 + &
           Qhat(k,i,j,3)**2) 

   enddo
enddo
end 
   
   
   
