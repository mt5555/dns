#include "macros.h"
!
! test cases that require n_var >= 3
!

subroutine init_data_lwisotropic(Q,PSI,work,work2)
!
! low wave number, quasi isotropic initial condition
!
use params
use mpi
use transpose
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: PSI(nx,ny,nz,n_var)
real*8 :: work(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

! local variables
real*8 :: alpha,beta
integer km,jm,im,i,j,k,n,wn,ierr
integer,allocatable :: seed(:)
real*8 xw,ener1,ener2,ener1_target,ener2_target,ener,xfac
character(len=80) message
CPOINTER :: null

!
! 
!   g_xcord(i)=(i-1)*delx	

!random vorticity

! set the seed - otherwise it will be the same for all CPUs,
! producing a bad initial condition
call random_seed(size=k)
allocate(seed(k))
call random_seed(get=seed)
seed=seed+my_pe
call random_seed(put=seed)
deallocate(seed)

do n=1,3
   ! input from random number generator
   ! this gives same I.C independent of cpus
   call input1(PSI(1,1,1,n),work2,work,null,io_pe,.true.)  
enddo



alpha=0
beta=1
do n=1,3
   call helmholtz_periodic_inv(PSI(1,1,1,n),work,alpha,beta)
enddo
call vorticity(Q,PSI,work,work2)


ener1=0
ener2=0
ener1_target=1.0**(-5.0/3.0)
ener2_target=2.0**(-5.0/3.0)

do n=1,3
   call fft3d(Q(1,1,1,n),work) 
   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)
            xw=sqrt(real(km**2+jm**2+im**2))

            xfac = (2*2*2)
            if (km==0) xfac=xfac/2
            if (jm==0) xfac=xfac/2
            if (im==0) xfac=xfac/2
            if (xw>=.5 .and. xw<1.5) then
               ener1=ener1+.5*xfac*Q(i,j,k,n)**2
            else if (xw>=1.5 .and. xw<2.5) then
               ener2=ener2+.5*xfac*Q(i,j,k,n)**2
            else
               Q(i,j,k,n)=0
            endif
         enddo
      enddo
   enddo
enddo

#ifdef USE_MPI
   ener=ener1
   call MPI_allreduce(ener,ener1,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   ener=ener2
   call MPI_allreduce(ener,ener2,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif

if (ener1<=0) then
   call abort("ener1 must be positive")
endif
if (ener2<=0) then
   call abort("ener2 must be positive")
endif


do n=1,3
   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)
            xw=sqrt(real(km**2+jm**2+im**2))

            if (xw>=.5 .and. xw<1.5) then
               Q(i,j,k,n)=Q(i,j,k,n)*sqrt(ener1_target/(ener1))
            else if (xw>=1.5 .and. xw<2.5) then
               Q(i,j,k,n)=Q(i,j,k,n)*sqrt(ener2_target/(ener2))
            endif
         enddo
      enddo
   enddo
enddo

ener=0
ener1=0
ener2=0
do n=1,3
   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)
            xw=sqrt(real(km**2+jm**2+im**2))

            xfac = (2*2*2)
            if (km==0) xfac=xfac/2
            if (jm==0) xfac=xfac/2
            if (im==0) xfac=xfac/2
            if (xw>=.5 .and. xw<1.5) then
               ener1=ener1+.5*xfac*Q(i,j,k,n)**2
            else if (xw>=1.5 .and. xw<2.5) then
               ener2=ener2+.5*xfac*Q(i,j,k,n)**2
            endif
            ener=ener+.5*xfac*Q(i,j,k,n)**2


         enddo
      enddo
   enddo
enddo
do n=1,3
   call ifft3d(Q(1,1,1,n),work) 
enddo
#ifdef USE_MPI
   xfac=ener
   call MPI_allreduce(xfac,ener,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   xfac=ener1
   call MPI_allreduce(xfac,ener1,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   xfac=ener2
   call MPI_allreduce(xfac,ener2,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif

call print_message("Isotropic initial condition in wave numbers:");
write(message,'(a,f8.4,a,f8.4)') "0.5 <= wn < 1.5   E=",ener1,"  E target=",ener1_target
call print_message(message)
write(message,'(a,f8.4,a,f8.4)') "1.5 <= wn < 2.5   E=",ener2,"  E target=",ener2_target
call print_message(message)
write(message,'(a,f8.4,a,f8.4)') "Total E=",ener
call print_message(message)




end subroutine





subroutine init_data_sht(Q,PSI,work1,work2,init)
!
!   init=0   setup paraemeters only
!   init=1   setup paraemeters and compute initial condition
!
!
use params
use transpose
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: PSI(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
integer :: init

! local variables
CPOINTER :: null
integer i,j,k,l
integer :: m,k_0
integer :: n,im,jm,km,ixw,ierr
real*8 :: xw,xfac ,alpha,beta,dummy
real*8 :: E_target,ke,pe,ke2,E_0
real*8 :: E_k(0:max(nx,ny,nz))
real*8 :: E_k2(0:max(nx,ny,nz))
real*8 :: Len,U,R,F
character(len=80) :: message

equations = SHALLOW   ! shallow water equations, not the default NS equations

k_0=14  
m=25
U=1.0        ! velocity scale  U normalized below so that KE=.5 U**2 = .5
Len=1.0/(2*pi*k_0)  ! length scale, determined by mode k_0 above



! set 
! R = Rossby number
! F = Froude number
if (init_cond_subtype==0) then        ! run A
   R=.01
   F=.04
else if (init_cond_subtype==1) then   ! run B
   R=.05
   F=.05
else if (init_cond_subtype==2) then   ! run C
   R=.05
   F=.075
else if (init_cond_subtype==3) then   ! run D
   R=.25
   F=.05
else if (init_cond_subtype==4) then   ! run E
   R=.25
   F=.20
else if (init_cond_subtype==5) then   ! run F
   R=1.00
   F=0.05
else if (init_cond_subtype==6) then   ! run G
   R=1.00
   F=0.30
else if (init_cond_subtype==7) then   ! run H
   R=5.00
   F=0.05
else if (init_cond_subtype==8) then   ! run I
   R=5.00
   F=0.30
else if (init_cond_subtype==9) then   ! run J
   R=20.00
   F=00.30
else if (init_cond_subtype==10) then   ! run K
   R=10.00
   F=00.10
else if (init_cond_subtype==11) then   ! run L
   R=20.00
   F=00.05
else if (init_cond_subtype==12) then   ! run M
   R=2.00
   F=.1
else if (init_cond_subtype==13) then   ! run N
   R=.4
   F=.1
else
   call abort("init_data_swt(): init_cond_subtype set to unsupported value")
endif
   

! set dimensionless constants used in code: f-plane and gravity:
H0 = 1.0
fcor=U/(Len*R)
grav=U**2/(H0*F**2)              

write(message,'(a,2f6.2)') "R, F parameters: ",R,F
call print_message(message)

if (init==0) return


Q=0
! random vorticity
call input1(PSI,work1,work2,null,io_pe,.true.)  

! invert to get PSI
alpha=0
beta=1
call helmholtz_periodic_inv(PSI,work1,alpha,beta)

! initialize amplitude to satisfy:  E(k)/.5k^2
! E(k) = k^(m/2) / (k + k_0)^m
!
!


n=1
call fft3d(PSI,work1) 
E_k=0
do k=nz1,nz2
   km=kmcord(k)
   do j=ny1,ny2
      jm=jmcord(j)
      do i=nx1,nx2
         im=imcord(i)
         xw=sqrt(real(km**2+jm**2+im**2))
         ixw = nint(xw)
         
         xfac = (2*2*2)
         if (km==0) xfac=xfac/2
         if (jm==0) xfac=xfac/2
         if (im==0) xfac=xfac/2

         E_k(ixw)=E_k(ixw)+.5*pi2_squared*xw*xw*xfac*PSI(i,j,k,n)**2
      enddo
   enddo
enddo

#ifdef USE_MPI
   E_k2=E_k
   call MPI_allreduce(E_k2,E_k,max(nx,ny,nz),MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif

E_0=0
do k=nz1,nz2
   km=kmcord(k)
   do j=ny1,ny2
      jm=jmcord(j)
      do i=nx1,nx2
         im=imcord(i)
         xw=sqrt(real(km**2+jm**2+im**2))
         ixw = nint(xw)
         
         xfac = (2*2*2)
         if (km==0) xfac=xfac/2
         if (jm==0) xfac=xfac/2
         if (im==0) xfac=xfac/2

         if (xw>0) then
            E_target=(xw**(m/2.0) / (xw+k_0)**m) 

            PSI(i,j,k,n)=PSI(i,j,k,n)/sqrt(E_k(ixw))
            PSI(i,j,k,n)=PSI(i,j,k,n)*sqrt(E_target)
         else
            PSI(i,j,k,n)=0
         endif

         if (k_0==ixw) then
            E_0 = E_0 + .5*xfac*pi2_squared*xw**2*PSI(i,j,k,n)**2
         endif
      enddo
   enddo
enddo


call ifft3d(PSI,work1) 




! compute velocity = curl PSI
Q=0
call der(PSI,Q,dummy,work1,DX_ONLY,2)
Q(:,:,:,1)=-Q(:,:,:,1)
call der(PSI,Q(1,1,1,2),dummy,work1,DX_ONLY,1)

if (dealias) call dealias_gridspace(Q,work1)


#if 0
! normalize so <u,u>=U**2
ke = .5*sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,1)**2) +  &
     .5*sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,2)**2)  

#ifdef USE_MPI
   ke2=ke
   call MPI_allreduce(ke2,ke,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
ke=ke/g_nx/g_ny
Q=sqrt(.5*U**2) * Q/sqrt(ke)
E_0 = E_0*.5*U**2/ke
print *,'E_0 = ',E_0

#else

! normalize so Energy in band k_0 = 3.23448943961018D-002
Q =   sqrt(3.23448943961018D-002)*Q/sqrt(E_0)

#endif






if (dealias) call dealias_gridspace(Q,work1)

! compute a balenced hight field:
! div( u grad u + f cross u  + g grad h ) = 0
! 
! PSI = u grad u + f cross u 
PSI=0
PSI(:,:,:,1)= fcor*Q(:,:,:,2)
PSI(:,:,:,2)=-fcor*Q(:,:,:,1)

! advection and viscous terms
k=1
do n=1,2

   ! compute u_x
   call der(Q(1,1,1,n),work1,dummy,work2,DX_ONLY,1)
   do j=ny1,ny2
   do i=nx1,nx2
      PSI(i,j,k,n) = PSI(i,j,k,n) - Q(i,j,k,1)*work1(i,j,k) 
   enddo
   enddo


   ! compute u_y
   call der(Q(1,1,1,n),work1,dummy,work2,DX_ONLY,2)
   do j=ny1,ny2
   do i=nx1,nx2
      PSI(i,j,k,n) = PSI(i,j,k,n) - Q(i,j,k,2)*work1(i,j,k) 
   enddo
   enddo

enddo

if (n_var==3) then
   call divfree_gridspace(PSI,Q(1,1,1,n_var),work1,work2)
   Q(:,:,:,n_var)= H0 + Q(:,:,:,n_var)/grav
else
   call abort("init_data_sht() requires n_var==3")
endif



end subroutine







