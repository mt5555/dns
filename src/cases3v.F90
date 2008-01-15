#include "macros.h"
!
! test cases that require n_var >= 3
!

subroutine init_data_lwisotropic(Q,PSI,work,work2,init,rantype)
!
! low wave number, quasi isotropic initial condition
! init=0         set parameters (but there are none to set)
!                but do not compute initial condition
!
! init= 1        compute initial condition
!
! rantype=  0    intialize using a gaussian, in grid space
! rantype=  1    initalize using E=constant and random phase
!
use params
use mpi
use transpose
implicit none
integer :: init,rantype
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: PSI(nx,ny,nz,n_var)
real*8 :: work(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

! local variables
real*8 :: alpha,beta
integer km,jm,im,i,j,k,n,wn,ierr,nb
integer,parameter :: NUMBANDS=15
real*8 xw,enerb(NUMBANDS),enerb_target(NUMBANDS),ener,xfac,theta
real*8 enerb_work(NUMBANDS)
character(len=80) message
CPOINTER :: null


if (init==0) return


! compute U using random vorticity
! rantype=0   reproducable with different parallel decompositions, slow
! rantype=1   fast, not reproducable
call ranvor(Q,PSI,work,work2,rantype)



!
!  rescale data to fit energy profile
!
enerb=0
enerb_target=0
if (init_cond_subtype==0) then
do nb=1,NUMBANDS
   enerb_target(nb)=real(nb)**(-5.0/3.0)
   if (nb>2)  enerb_target(nb)=real(nb)**(-7.0/3.0)
enddo
endif

if (init_cond_subtype==1) then
! Balu 2D initial condition, with low waver number hypo viscosity
do nb=1,12
   enerb_target(nb)=.063
enddo
endif


! rescale E to fit enerb_target():
call rescale_e(Q,work,ener,enerb,enerb_target,NUMBANDS,3)


call print_message("Isotropic initial condition in wave numbers:");
do nb=1,min(NUMBANDS,max(g_nx/2,g_ny/2,g_nz/2))
   write(message,'(a,i4,a,e12.4,a,e12.4)') "wn=",nb,"+/-.5   E=",enerb(nb),&
        "  E target=",enerb_target(nb)
   call print_message(message(1:79))
enddo
write(message,'(a,f8.4,a,f8.4)') "Total E=",ener
call print_message(message)


end subroutine





subroutine init_data_decay(Q,PSI,work,work2,init,rantype,restype)
!
! low wave number, quasi isotropic initial condition
! based on spectrum provided by Menevea
!
! init=0         set parameters 
!                but do not compute initial condition
!
! init= 1        compute initial condition
!
! init=2         rescale energy of initial data in Q
!
!
! rantype=  0    intialize using a gaussian, in grid space
! rantype=  1    initalize using E=constant and random phase
!                  (faster, but not independed of parallel decomp)
!
! restype =      not used
!
! init_cond_subtype: 
!     0               use 2048^3 parameters with kmax eta = 1
!     1               use Livescu initial spectrum, peaked at 10
!     2               use Tzaki initial spectrum
!     3               Livescu spectrum, peaked at 6
!     4               k**-4 spectrum
!     5		      controlled helicity in each mode
!

use params
use mpi
use transpose
implicit none
integer :: init,rantype,restype
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: PSI(nx,ny,nz,n_var)
real*8 :: work(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

! local variables
real*8 :: alpha,beta
integer km,jm,im,i,j,k,n,wn,ierr,nb,NUMBANDS
integer,parameter :: NUMBANDS_MAX=1024
real*8 xw,ener,xfac,theta
real*8 ::  enerb_target(NUMBANDS_MAX),enerb_work(NUMBANDS_MAX),enerb(NUMBANDS_MAX)
real*8 :: xnb,ep23,lenscale,eta,fac1,fac2,fac3,ens,epsilon,mu_m
real*8 :: tmx1,tmx2,p,p1,p2,ku,h_angle

character(len=80) message,fname
CPOINTER :: null



NUMBANDS=NUMBANDS_MAX
if (dealias==1) then
   NUMBANDS = g_nmin/3
else if (dealias==2) then
  NUMBANDS = dealias_sphere_kmax
else if (dealias==3) then
   NUMBANDS = dealias_23sphere_kmax
endif
!write(6,*)'NUMBANDS = ',NUMBANDS

!
! cj2.out run: 128^3
!  R_l = 73   Teddy=.95  Umax=4.25  delx/eta = 1.9
! Teddy/delt = Teddy*Umax
!
if (init_cond_subtype==0) then
   ep23=(2.76029)  **(2d0/3d0) 
   lenscale = .3155
   eta = 9.765e-4
   mu_m=1.359e-4 ! *.7  !   *.89
   do nb=1,NUMBANDS
      xnb = nb
      fac1 = 1.613*ep23*xnb**(-5d0/3d0)
      
      fac2 = ( (xnb*lenscale)**1.2 + .39 )**0.833 
      fac2 = (xnb*lenscale/fac2)**(17d0/3d0)
      fac2 = fac2 * exp(-2.1*xnb*eta)
      
      fac3=.5 + atan( 10*log10(xnb*eta)+12.58 ) / pi
      fac3 = 1 + .522*fac3
      
      enerb_target(nb)=fac1*fac2*fac3
   enddo
else if (init_cond_subtype==1) then
   ! slope: k**2, peak k=10
   call livescu_spectrum(enerb_target,NUMBANDS,0,init_cond_subtype)
else if (init_cond_subtype==2) then
   do nb=1,NUMBANDS
	enerb_target(nb) = nb**(-5./3.)
   enddo
else if (init_cond_subtype==3) then
   ! slope: k**2, peak k=6
   call livescu_spectrum(enerb_target,NUMBANDS,0,init_cond_subtype)
else if (init_cond_subtype==4) then
   ! initial Energy spectrum for Evelyn's 2D forced runs
   do nb=1,NUMBANDS
	enerb_target(nb) = 10.0 * nb**(-4)
   enddo
else if (init_cond_subtype==5) then
! initial energy spectrum, total energy 1,  with controlled helicity 
! in each mode
   ener = 0
   do nb = 1,NUMBANDS
        enerb_target(nb) = nb**2 
	ener = ener + enerb_target(nb)
   enddo
   enerb_target = enerb_target/ener	
else
   call abortdns("init_data_decay: bad init_cond_subtype")
endif

if (init==0) return

if (my_pe==io_pe .and. init_cond_subtype==0) then 
   ! compute some stats:
   ener=0
   ens=0
   do nb=1,NUMBANDS
      ener=ener+enerb_target(nb)
      ens=ens+2 * nb**2 * enerb_target(nb)
   enddo
   epsilon=mu_m*ens
   
   print *,'stats computed from target spectrum and target viscosity:'
   print *,'epsilon: ',epsilon
   print *,'energy:  ',ener
   print *,'u1       ',sqrt(2*ener/3)
   print *,'u1,1     ',epsilon/mu_m/15
   print *,'lambda   ',sqrt(10*ener*mu_m/epsilon)
   print *,'R_l      ',sqrt(10*ener*mu_m/epsilon) * sqrt(2*ener/3)/mu_m
   print *,'eta      ',(mu_m**3 / epsilon ) **.25 
   print *,'965*eta  ',965*(mu_m**3 / epsilon ) **.25 
   print *,'eddy time',2*ener/epsilon

   ! convert to our units, with L=1 instead of 2pi
   enerb_target=enerb_target / (2*pi*2*pi)
endif




! random vorticity initial condition:
if (init==1) then
   !
   !  We are computing a initial condition from scratch:
   !
   call ranvor(Q,PSI,work,work2,rantype)

   ! rescale energy to match enerb_target, preserving the phases
   call rescale_e(Q,work,ener,enerb,enerb_target,NUMBANDS,3)
   
   ! If using controlled helicity initial condition, 
   ! set the helicity angle to h_angle
   if (init_cond_subtype == 5) then
      h_angle = init_cond_param1
      call set_helicity_angle(Q,PSI,work,h_angle,ener)
   endif

else if (init==2) then
   !
   !  This subroutine was called with an initial condition in Q
   !  (from a restart file).  Modify it for this run:  
   !
   if ( 0 <= init_cond_subtype .and. init_cond_subtype <= 4) then
      ! rescale energy to match enerb_target, preserving the phases
      ! used during the intialization procedure for decaying turbulence runs
      call rescale_e(Q,work,ener,enerb,enerb_target,NUMBANDS,3)
   endif

   ! If using controlled helicity initial condition, 
   ! set the helicity angle to h_angle
   if (init_cond_subtype == 5) then
      h_angle = init_cond_param1
      call set_helicity_angle(Q,PSI,work,h_angle,ener)
   endif
endif


call print_message("Isotropic initial condition in wave numbers:");
do nb=1,max(10,NUMBANDS,max(g_nx/2,g_ny/2,g_nz/2))
   write(message,'(a,i4,a,e12.4,a,e12.4)') "wn=",nb,"+/-.5   E=",enerb(nb),&
        "  E target=",enerb_target(nb)
   call print_message(message)
enddo
write(message,'(a,e16.4)') "Total E=",ener
call print_message(message)


if (my_pe==io_pe) then 
! compute some stats:
ener=0
ens=0
do nb=1,NUMBANDS
   ener=ener+enerb(nb)
   ens=ens+2 * (2*pi*nb)**2 * enerb(nb)
enddo
mu_m=mu
epsilon=mu_m*ens

print *,'stats computed from initial condition and our viscosity:'
print *,'epsilon: ',epsilon, epsilon*pi2_squared
print *,'energy:  ',ener, ener*pi2_squared
print *,'u1       ',sqrt(2*ener/3), sqrt(2*ener/3)*pi2 
print *,'u1,1     ',epsilon/mu_m/15
print *,'lambda   ',sqrt(10*ener*mu_m/epsilon),sqrt(10*ener*mu_m/epsilon)*pi2
print *,'R_l      ',sqrt(10*ener*mu_m/epsilon) * sqrt(2*ener/3)/mu_m
print *,'eta      ',    (mu_m**3 / epsilon ) **.25, &
                    pi2*(mu_m**3 / epsilon ) **.25 
print *,'delx/eta ',  (1./g_nmin)/((mu_m**3 / epsilon ) **.25)
print *,'eddy time',2*ener/epsilon
endif



end subroutine







subroutine init_3d_rot(Q,PSI,work1,work2,init)
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
real*8 :: Omega
integer :: init,nb
character(len=240) :: message
integer,parameter :: NUMBANDS=100
real*8  :: enerb_target(NUMBANDS),enerb(NUMBANDS),ener


if (init==0) return


! compute U using random vorticity
call ranvor(Q,PSI,work1,work2,0)
!
!  rescale data to fit energy profile
!
enerb=0
enerb_target=0
if (init_cond_subtype==0) then
do nb=1,NUMBANDS
   enerb_target(nb)=.1*exp(-.5*(nb-forcing_peak_waveno)**2)/sqrt(2*pi)
enddo
endif

! rescale E to fit enerb_target():
call rescale_e(Q,work1,ener,enerb,enerb_target,NUMBANDS,3)


call print_message("Isotropic initial condition in wave numbers:");
do nb=1,min(NUMBANDS,max(g_nx/2,g_ny/2,g_nz/2))
   write(message,'(a,i4,a,e12.4,a,e12.4)') "wn=",nb,"+/-.5   E=",enerb(nb),&
        "  E target=",enerb_target(nb)
   call print_message(message(1:79))
enddo
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
logical :: init_zero=.false.

k_0=14  
if (equations==NS_UVW) then
   k_0=100
endif
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
else if (init_cond_subtype==14) then   ! run O
   R=.01
   F=.01
else if (init_cond_subtype==15) then   ! run P
   R=.001
   F=.001
else if (init_cond_subtype==101) then   ! run N
   R=.05
   F=.05
   init_zero=.true.
else if (init_cond_subtype==105) then   ! run N
   R=1.0
   F=.05
   init_zero=.true.
else
   call abortdns("init_data_swt(): init_cond_subtype set to unsupported value")
endif
   

! set dimensionless constants used in code: f-plane and gravity:
H0 = 1.0
fcor=U/(Len*R)
grav=U**2/(H0*F**2)              

write(message,'(a,2f6.2)') "R, F parameters: ",R,F
call print_message(message)

if (init==0) return


Q=0

! use zero initial condition?
if (init_zero) then
   Q(:,:,:,n_var)= H0 
   return
endif

! random vorticity
call input1(PSI(1,1,1,1),work1,work2,null,io_pe,.true.,-1)  

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
   call mpi_allreduce(E_k2,E_k,max(nx,ny,nz),MPI_REAL8,MPI_SUM,comm_3d,ierr)
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
! u = PSI_y
call der(PSI,Q,dummy,work1,DX_ONLY,2)     
! v = -PSI_x
call der(PSI,Q(1,1,1,2),dummy,work1,DX_ONLY,1)
Q(:,:,:,2)=-Q(:,:,:,2)

if (dealias>0) call dealias_gridspace(Q,work1)



#ifdef USE_MPI
   ke2=E_0
   call mpi_allreduce(ke2,E_0,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
! normalize so Energy in band k_0 = 3.23448943961018D-002
Q =   sqrt(3.23448943961018D-002)*Q/sqrt(E_0)





if (equations==SHALLOW) then
! compute a balenced hight field:
! div( u grad u + f cross u  + g grad h ) = 0
! 
! PSI = u grad u + f cross u 
PSI=0
PSI(:,:,:,1)= fcor*Q(:,:,:,2)
PSI(:,:,:,2)=-fcor*Q(:,:,:,1)

! advection and viscous terms
k=nz1
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

if (n_var>=3) then
   call divfree_gridspace(PSI,Q(1,1,1,n_var),work1,work2)
   do j=ny1,ny2
   do i=nx1,nx2
      Q(i,j,k,n_var)= H0 + Q(i,j,k,n_var)/grav
      if (Q(i,j,k,n_var)<.05*H0)  Q(i,j,k,n_var)=.05*H0
   enddo
   enddo
else
   call abortdns("init_data_sht() with shallow water equations requires n_var==3")
endif
endif



end subroutine







subroutine livescu_spectrum(enerb_target,numb,itype,init_cond_subtype)
!
! compute E(k)
!   itype==0   velocity spectrum
!   itype==1   scalar spectrum
!
!
implicit none
integer :: numb,itype,init_cond_subtype
real*8 :: enerb_target(numb)

real*8 :: p,p1,p2,ku,xnb,fac1,fac2,fac3
integer ::  nb

! default values:
p=2;
ku=10

if (init_cond_subtype==1) then
   p=2;
   ku=10
endif
if (init_cond_subtype==3) then
   p=2;
   ku=6
endif



enerb_target=0
p1=(1+p)/2
p2=p/2


do nb=1,numb
   xnb = nb
   fac1=p2**p1 * xnb**p
   fac2=ku**(p+1)
   fac3=exp(-p2* (xnb/ku)**2)
   enerb_target(nb)=fac1*fac3/fac2
enddo
end subroutine livescu_spectrum


subroutine set_helicity_angle(Q,Qi,work,h_angle,ener)
!
!set the helicity angle in each mode for prescribed spectrum shape.
!
!First convert to R,I formulation of fourier modes (Polifke and Shtilman)
!Set angle between R and I
!Convert back to real formulation of fourier modes
!

use params
use mpi
implicit none

real*8 :: Q(nx,ny,nz,3)
real*8 :: QI(nx,ny,nz,3)
real*8 :: work(nx,ny,nz)
real*8 :: h_angle,ener  ! energy (input argument)
integer :: n,wn,fix,iii,jjj,kkk
integer :: i,j,k
real*8 :: xw2,xw
real*8 RR(3),II(3), IIp(3), RRhat(3), khat(3), yhat(3),RRxk_hat(3), IIxk_hat(3),RRxIIp(3)
real*8 mod_ii, mod_rr, mod_IIp, mod_RRxk, RRdotk, IIdotk, RRdotII, RRdotIIp, mod_IIxk
real*8 costta, tta, tta_sgn, theta, phi 
real*8 :: cos_h_angle,sin_h_angle
character(len=80) :: message

integer :: zerosign
external :: zerosign
integer :: i1,i2,k1,k2,j1,j2

if (ndim /= 3) then
   call abortdns("ERROR: set_helicity_angle() requires ndim=3")
endif

cos_h_angle=cos(h_angle)
sin_h_angle=sin(h_angle)


write(message,*) 'Setting helicity angle = ',h_angle
call print_message(message)


! take FFT, then convert to complex coefficients
! real part is stored in Q array
! complex part is stored in QI array
do n = 1,3
   write(message,*) 'converting gridspace to to complex modes, n=',n   
   call print_message(message)
   call fft3d(Q(1,1,1,n),work)    ! inplace FFT
   work=Q(:,:,:,n)                ! copy fft(q) to work array 
   call sincos_to_complex_field(work,Q(1,1,1,n),QI(1,1,1,n))  ! convert
enddo


!     apply helicity fix:
do kkk=nz1,nz2
do jjj=ny1,ny2
do iii=nx1,nx2
   ! note: wave number is (i,j,k)
   ! index into arrays is (iii,jjj,kkk)
   !
   ! to be more consitent with other loops, it would be better
   ! to have the array index (i,j,k) and wave numbers (im,jm,km) 
   !
   i=imcord_exp(iii)
   j=jmcord_exp(jjj)
   k=kmcord_exp(kkk)
   RR = Q(iii,jjj,kkk,:)
   II = Qi(iii,jjj,kkk,:)

   if (my_pe == io_pe) then 
      if (i == 1 .and. j == 1 .and. k == 1) then
         write(6,*)'Pre i,j,k  ',i,j,k
         write(6,*) RR
	 write(6,*) II
      else if (i == -1 .and. j == -1 .and. k == -1) then    
         write(6,*)'Pre i,j,k  ',i,j,k
         write(6,*) RR
	 write(6,*) II
      endif
   endif

   xw2=i**2 + j**2 + k**2
   xw=sqrt(xw2)

   mod_rr = sqrt(RR(1)*RR(1) + RR(2)*RR(2) + RR(3)*RR(3))
   mod_ii = sqrt(II(1)*II(1) + II(2)*II(2) + II(3)*II(3))

   fix=1
   if (mod_rr == 0 .or. mod_ii == 0 .or. xw==0 .or. h_angle==-1) then
      fix = 0
   endif

   
   ! do the transformation if fix = 1
   if (fix == 1) then         

      RRhat = RR/mod_rr
      khat(1) = i/xw
      khat(2) = j/xw
      khat(3) = k/xw
      
      ! yhat = khat X RRhat (coordinate system with zhat=khat and xhat = RRhat
      
      yhat(1) = khat(2)*RRhat(3) - khat(3)*RRhat(2)
      yhat(2) = - (khat(1)*RRhat(3) - khat(3)*RRhat(1))
      yhat(3) = khat(1)*RRhat(2) - khat(2)*RRhat(1)
      
      ! frst check that RR and II are orthogonal to k (incompressibility).
      RRdotk = (RR(1)*i + RR(2)*j + RR(3)*k)
      IIdotk = (II(1)*i + II(2)*j + II(3)*k)
      if (abs(RRdotk/xw)>(1e-12*sqrt(ener)) .or. abs(IIdotk/xw)>(1e-12*sqrt(ener))) then
         print *,'WARNING:  RR and II not orthogonal to K:',i,j,k,xw
         write(6,*)'RRdotk, norm, KE = ',(RRdotk),mod_rr,sqrt(ener)
         write(6,*)'IIdotk, norm, KE = ',(IIdotk),mod_ii,sqrt(ener)
         write(6,*)'RRdotk_ang = ',acos(RRdotk/mod_rr/xw)
         write(6,*)'IIdotk_ang = ',acos(IIdotk/mod_ii/xw)
      endif
      
      ! Rotate II and RR into plane orthogonal to k
      ! Don't need to do this for Mark's isotropic forcing since incompressibility
      ! already ensures this. Comment out.
      
      !      mod_RRxk = sqrt((RR(2)*k - RR(3)*j)**2 + (RR(1)*k - RR(3)*i)**2 + (RR(1)*j - RR(2)*i)**2)
      !      RRxk_hat(1) = (RR(2)*k - RR(3)*j)/mod_RRxk
      !      RRxk_hat(2) = -(RR(1)*k - RR(3)*i)/mod_RRxk
      !      RRxk_hat(3) = (RR(1)*j - RR(2)*i)/mod_RRxk
      !      IIp = mod_ii * RRxk_hat
      
      !     mod_IIxk = sqrt((II(2)*k - II(3)*j)**2 + (II(1)*k - II(3)*i)**2 + (II(1)*j - II(2)*i)**2)
      !      IIxk_hat(1) = (II(2)*k - II(3)*j)/mod_IIxk
      !      IIxk_hat(2) = -(II(1)*k - II(3)*i)/mod_IIxk
      !      IIxk_hat(3) = (II(1)*j - II(2)*i)/mod_IIxk
      !      RR = mod_rr * IIxk_hat
      
      
      ! write II in terms of new coordinate system (z=khat,x=RRhat,y=yhat)
      IIp = II	
      II(1) = IIp(1)*RRhat(1) + IIp(2)*RRhat(2) + IIp(3)*RRhat(3)
      II(2) = IIp(1)*yhat(1) + IIp(2)*yhat(2) + IIp(3)*yhat(3)
      II(3) = IIp(1)*khat(1) + IIp(2)*khat(2) + IIp(3)*khat(3) !should be zero
      
      ! set II to have angle h_angle wrt RR
      IIp(1) = cos_h_angle*mod_ii
      IIp(2) = sin_h_angle*mod_ii
      IIp(3) = II(3)
      
      
      ! check angle between RR and IIp 
      !         acos(IIp(1)/mod_ii) == h_angle, or IIp(1)/mod_ii == cos(h_angle)
      !         if (abs(IIp(1) - mod_ii*(cos_h_angle)) >1e-10*mod_ii) then
      !            tta = acos(IIp(1)/mod_ii)
      !            print *,'h_angle = ',tta
      !            write(6,*)'postfix angle bet. RR and IIp = ', tta*180/pi
      !         endif
      
      ! now write II in old (i,j,k) coordinate system
      
      II(1) = IIp(1)*RRhat(1) + IIp(2)*yhat(1) + IIp(3)*khat(1)
      II(2) = IIp(1)*RRhat(2) + IIp(2)*yhat(2) + IIp(3)*khat(2)
      II(3) = IIp(1)*RRhat(3) + IIp(2)*yhat(3) + IIp(3)*khat(3)
      
      
      
      ! check angle between final RR and II again
      !         RRdotII = RR(1)*II(1) + RR(2)*II(2) + RR(3)*II(3)
      !         costta = (RRdotII/(mod_rr * mod_ii))
      !	write(6,*),'costta =',costta
      !         if (abs(costta - cos_h_angle) >1e-10) then
      !            tta = acos(costta)
      !           write(6,*),'h_angle = ',h_angle
      !            write(6,*)'postfix angle bet. RR and II = ', tta*180/pi
      !         endif
      
      Q(iii,jjj,kkk,:)=RR
      Qi(iii,jjj,kkk,:)=II
   endif

if (my_pe == io_pe) then
      if (i == 1 .and. j == 1 .and. k == 1) then
         write(6,*)'post i,j,k  ',i,j,k
         write(6,*) RR
	 write(6,*) II
      else if (i == -1 .and. j == -1 .and. k == -1) then
         write(6,*)'post i,j,k  ',i,j,k
         write(6,*) RR
	 write(6,*) II
      endif
   endif


enddo
enddo
enddo

!
!  The helicity adjustment above did not preserve the fact that
!  mode (l,m,n) needs to be the complex conjugate of mode (-l,-m,-n)
!  Reimpose this constraint:
!
#if 1
do n=1,3
   do kkk=nz1,nz2
      do jjj=ny1,ny2
         do iii=nx1,nx2
            i1=iii; i2 = i1 + imsign(iii)
            j1=jjj; j2 = j1 + imsign(jjj)
            k1=kkk; k2 = k1 + imsign(kkk)
            Q(i2,j2,k2,n)  =   Q(i1,j1,k1,n)
            Qi(i2,j2,k2,n)  = -Qi(i1,j1,k1,n)
         enddo
      enddo
   enddo
enddo
#endif

do kkk=nz1,nz2
do jjj=ny1,ny2
do iii=nx1,nx2
   i=imcord_exp(iii)
   j=jmcord_exp(jjj)
   k=kmcord_exp(kkk)

   if (my_pe == io_pe) then
      if (i == 1 .and. j == 1 .and. k == 1) then
         write(6,*) "Post  i,j,k  ",i,j,k
         write(6,*) Q(iii,jjj,kkk,:)
         write(6,*) Qi(iii,jjj,kkk,:)
      else if (i == -1 .and. j == -1 .and. k == -1) then
         write(6,*)'Post i,j,k  ',i,j,k
         write(6,*) Q(iii,jjj,kkk,:) 
         write(6,*) Qi(iii,jjj,kkk,:)
      endif
   endif
enddo
enddo
enddo




!    convert back:
do n = 1,3
   write(message,*) 'converting back to gridspace, n=',n   
   call print_message(message)
   call complex_to_sincos_field(work,Q(1,1,1,n),QI(1,1,1,n))
   Q(:,:,:,n)=work
   call ifft3d(Q(1,1,1,n),work)
enddo


end subroutine set_helicity_angle
