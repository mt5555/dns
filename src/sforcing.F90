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
!if (forcing_type==2) call sforcing_random12(rhs,Qhat,f_diss,0)

end subroutine




subroutine gforce(Q,rhs,rhsz,q4,q4z,work,f_diss)
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
real*8 :: f_diss

! local
integer :: n
real*8 :: fdiss

q4=Q
do n=1,3
   call z_fft3d_trashinput(q4(1,1,1,n),rhsz(1,1,1,n),work)
enddo
q4z=0
call sforce(q4z,rhsz,fdiss)
do n=1,3
   call z_ifft3d(q4z(1,1,1,n),rhs(1,1,1,n),work)
enddo
f_diss=fdiss  ! strange bug on SGI, dnsghost, f_diss is non zero here, but
              ! return value in ns_ghost.F90 is 0.  if we introduce 'fdiss',
              ! everthing is ok
return
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
real*8 :: Qhat(g_nz2,nslabx,ny_2dz,3) 
real*8 :: rhs(g_nz2,nslabx,ny_2dz,3) 
integer km,jm,im,i,j,k,n,wn,ierr
real*8 xw,xfac,f_diss,tauf
real*8 ener(NUMBANDS),ener_target(NUMBANDS),temp(NUMBANDS)
character(len=80) :: message

if (0==init_sforcing) then
   call sforcing_init()
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








subroutine sforcing_random12(rhs,Qhat,f_diss,new_f)
!
! Add a forcing term to rhs.
! Random, isotropic, homogenious in first 2 wave nubmers
!
! if new_f=0, add previously computed forcing stored in rmodes() into RHS. 
! if new_f=1, compute a new forcing ONLY, store in rmodes()
!
use params
use mpi
implicit none
integer :: new_f
real*8 :: Qhat(g_nz2,nslabx,ny_2dz,3) 
real*8 :: rhs(g_nz2,nslabx,ny_2dz,3) 
integer km,jm,im,i,j,k,n,wn,ierr
real*8 xw,xfac,f_diss,tauf
real*8 ener(NUMBANDS),ener_target(NUMBANDS),temp(NUMBANDS)

real*8,save :: rmodes(0:2,0:2,0:2,3)       ! value at time tmod
real*8,save :: rmodes_old(0:2,0:2,0:2,3)   ! value at time tmod_old
real*8,save :: tmod,tmod_old
real*8,save :: tscale=.01


if (0==init_sforcing) then
   call sforcing_init()
   rmodes=0
   tmod=0

   ! check that we are not including any wave numbers > 2
   do wn=1,NUMBANDS
   do n=1,wnforcing(wn)%n
      i=wnforcing(wn)%index(n,1)
      j=wnforcing(wn)%index(n,2)
      k=wnforcing(wn)%index(n,3)
      if (abs(z_imcord(i))>2 .or. abs(z_jmcord(j))>2 .or. abs(z_kmcord(k))>2)  then 
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
   call random12(rmodes)
   return
endif



f_diss=0
do wn=1,NUMBANDS

   ener(wn)=0
   do n=1,wnforcing(wn)%n
      i=wnforcing(wn)%index(n,1)
      j=wnforcing(wn)%index(n,2)
      k=wnforcing(wn)%index(n,3)
      rhs(i,j,k,1)=rhs(i,j,k,1) + rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),1)
      rhs(i,j,k,2)=rhs(i,j,k,2) + rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),2)
      rhs(i,j,k,3)=rhs(i,j,k,3) + rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),3)

      xfac=8
      if (z_kmcord(k)==0) xfac=xfac/2
      if (z_jmcord(j)==0) xfac=xfac/2
      if (z_imcord(i)==0) xfac=xfac/2
      f_diss = f_diss + xfac*( &
         Qhat(k,i,j,1)*rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),1) +&
         Qhat(k,i,j,2)*rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),2) +&
         Qhat(k,i,j,3)*rmodes(z_imcord(i),z_jmcord(j),z_kmcord(k),3) )


   enddo
enddo


end subroutine 
   
   





subroutine random12(rmodes)
use params
implicit none
real*8 :: rmodes(0:2,0:2,0:2,3)       

integer km,jm,im,i,j,k,n,wn,ierr
real*8 xw,xfac,f_diss,tauf
real*8 :: R(5*5*5,3,2),Rr,Ri
real*8 :: psix_r(3),psix_i(3)
real*8 :: psiy_r(3),psiy_i(3)
real*8 :: psiz_r(3),psiz_i(3)

rmodes=0
call gaussian(R,5*5*5*3*2)

k=0
do km=-2,2 
do jm=-2,2
do im=-2,2
   k=k+1
   !
   ! Choose R gaussian, theta uniform from [0..1]
   ! vorticty = (R1 + i R2)  exp(im*2pi*x) * exp(jm*2pi*y) * exp(km*2pi*z) 
   !
   ! convert from vorticity to stream function, then take curl:
   ! remove 1 factor of 2*pi from laplacian and derivative, since they candel
   xfac = (im**2 + jm**2 + km**2)
   if (xfac>0) xfac=1/(-2*pi*xfac)
   !if (delt>0) xfac=xfac/sqrt(delt)
   xfac=xfac/100
   
   do n=1,3
      psix_r(n) = -im*R(k,n,2)*xfac
      psix_i(n) =  im*R(k,n,1)*xfac
      psiy_r(n) = -jm*R(k,n,2)*xfac
      psiy_i(n) =  jm*R(k,n,1)*xfac
      psiz_r(n) = -km*R(k,n,2)*xfac
      psiz_i(n) =  km*R(k,n,1)*xfac
   enddo

   !
   !   
   ! convert to sine & cosine modes:
   !
   ! (R1 + i R2) (cosx + i sinx)  (cosy + i siny)  (cosz + i sinz)  
   ! = (real parts only:)
   !   R1  cosx cosy cosz                      R1 (1,1,1)
   ! i R2  cosx cosy sinz  i                  -R2 (1,1,-1) sign(km) 
   ! i R2  cosx siny cosz  i                  -R2 (1,-1,1) sign(jm)
   !   R1  cosx siny sinz  i**2               -R1 (1,-1,-1) sign(km*jm)
   ! i R2  sinx cosy cosz  i                  -R2 (-1,1,1) sign(im) 
   !   R1  sinx cosy sinz  i**2               -R1 (-1,1,-1) sign(im*km)
   !   R1  sinx siny cosz  i**2               -R1 (-1,-1,1) sign(im*jm)
   ! i R2  sinx siny sinz  i**3                R2 (-1,-1,-1) sign(im*jm*km)
   !  
   ! 
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

      rmodes( im, jm, km,n) = rmodes( im, jm, km,n) + Rr 
      rmodes( im, jm,-km,n) = rmodes( im, jm,-km,n) - Ri*sign(1,km) 
      rmodes( im,-jm, km,n) = rmodes( im,-jm, km,n) - Ri*sign(1,jm) 
      rmodes( im,-jm,-km,n) = rmodes( im,-jm,-km,n) - Rr*sign(1,jm*km) 
      rmodes(-im, jm, km,n) = rmodes(-im, jm, km,n) - Ri*sign(1,im)
      rmodes(-im, jm,-km,n) = rmodes(-im, jm,-km,n) - Rr*sign(1,im*km) 
      rmodes(-im,-jm, km,n) = rmodes(-im,-jm, km,n) - Rr*sign(1,im*jm) 
      rmodes(-im,-jm,-km,n) = rmodes(-im,-jm,-km,n) + Ri*sign(1,im*jm*km) 
   enddo

enddo
enddo
enddo
end subroutine



end module
