subroutine sforcing12(rhs,Qhat)
!
! Add a forcing term to rhs.
! Force 3D wave numbers 1 back to the sphere E=1**(-5/3)
! Force 3D wave numbers 2 back to the sphere E=2**(-5/3)
!
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: PSI(nx,ny,nz,n_var)
real*8 :: work(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
integer km,jm,im,i,j,k,n,wn
real*8 xw
logical,save :: firstcall=.true.

type wnforcing_d
   integer :: n
   integer, allocatable :: index(:,:)
end type
type(wnforcing_d) :: wnforcing(2)


if (firstcall) then
   firstcall=.false.

   do n=1,2
      wnforcing(n)%n=0
   enddo
   
   ! count the number of wavenumbers in each band.  
   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)
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
   do n=1,2
      i=wnforcing(n)%n
      allocate(wnforcing(n)%index(i,3))
      wnforcing(n)%n=0  ! reset counter to use again below
   enddo
   
   ! store all the indexes
   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)
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


do wn=1,2
   ener_target=wn**(-5.0/3.0)
   ener=0
   do n=1,wnforcing(wn)%n
      i=wnforcing(wn)%index(n,1)
      j=wnforcing(wn)%index(n,1)
      k=wnforcing(wn)%index(n,1)
      xfac=.5*8
      if (kmcord(k)==0) xfac=xfac/2
      if (jmcord(j)==0) xfac=xfac/2
      if (imcord(i)==0) xfac=xfac/2
      
      ener=ener+xfac*(Q(i,j,k,1)**2+Q(i,j,k,2)**2+Q(i,j,k,3)**2)
   enddo
   ! Qf = Q*sqrt(ener_target/ener)
   ! forcing = tau (Qf-Q) = tau * (sqrt(ener_target/ener)-1) Q
   rhs(i,j,k,1) = rhs(i,j,k,1) + tau*(sqrt(ener_target/ener)-1)*Q(i,j,k,1)
   rhs(i,j,k,2) = rhs(i,j,k,2) + tau*(sqrt(ener_target/ener)-1)*Q(i,j,k,2)
   rhs(i,j,k,3) = rhs(i,j,k,3) + tau*(sqrt(ener_target/ener)-1)*Q(i,j,k,3)



enddo

end 
   
   
   
