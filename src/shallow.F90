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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine to take one Runge-Kutta 4th order time step
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk4(time,Q_grid,Q,rhs,Q_old,Q_tmp,work1,work2)
use params
implicit none
real*8 :: time
real*8 :: Q(nx,ny,n_var)
real*8 :: Q_grid(nx,ny,n_var)
real*8 :: rhs(nx,ny,n_var)
real*8 :: Q_tmp(nx,ny,n_var)
real*8 :: Q_old(nx,ny,n_var)
real*8 :: work1(nx,ny)
real*8 :: work2(nx,ny)

! local variables
real*8 :: ke,pe
real*8 :: ints_buf(nints),vel,dtf,epfilt=.01
#undef USE_LEAPFROG
#ifdef USE_LEAPFROG
real*8,save :: QS(nx,ny,n_var)
real*8,save :: QM(nx,ny,n_var)
#endif
integer i,j,k,n,ierr
logical,save :: firstcall=.true.
integer,save :: itime=0


if (firstcall) then
   firstcall=.false.
   Q=Q_grid
   do n=1,3
      call fft3d(Q(1,1,n),work1)
      if (dealias>0) call fft_filter_dealias(Q(1,1,n))
   enddo
   if (equations/=SHALLOW) then
      call print_message("Error: shallow water model can only run equations=SHALLOW")
      call abortdns("initial conditions are probably incorrect.")
   endif
   if (ndim/=2) then
      call abortdns("Error: shallow water model cannot run in 3D")
   endif
   if (nz/=1) then
      call abortdns("Error: shallow water model cannot run in 3D")
   endif
endif





#ifndef USE_LEAPFROG

Q_old=Q

! stage 1
call getrhs(1,rhs,Q,Q_grid,time,1,work1,work2)
Q=Q+delt*rhs/6.0

! stage 2
Q_tmp = Q_old + delt*rhs/2.0
Q_grid=Q_tmp
do n=1,3
   call ifft3d(Q_grid(1,1,n),work1)
enddo
call getrhs(2,rhs,Q_tmp,Q_grid,time+delt/2.0,0,work1,work2)
Q=Q+delt*rhs/3.0

! stage 3
Q_tmp = Q_old + delt*rhs/2.0
Q_grid=Q_tmp
do n=1,3
   call ifft3d(Q_grid(1,1,n),work1)
enddo
call getrhs(3,rhs,Q_tmp,Q_grid,time+delt/2.0,0,work1,work2)
Q=Q+delt*rhs/3.0


! stage 4
Q_tmp = Q_old + delt*rhs
Q_grid=Q_tmp
do n=1,3
   call ifft3d(Q_grid(1,1,n),work1)
enddo
call getrhs(4,rhs,Q_tmp,Q_grid,time+delt,0,work1,work2)
Q=Q+delt*rhs/6.0



Q_grid=Q
do n=1,3
   call ifft3d(Q_grid(1,1,n),work1)
enddo

time = time + delt


#else




!  takes one Euler step with DTF=DT/2
!  then one leapfrog step with DTF=DT
!  and then leapfrog & filter with DTF=2*DT
!
!  Q     value at t+delt    (after time step complete)
!  QS    value at t-delt
!  QM    value at t
!

! ignore initial timesteps with delt==0
if (delt>0) itime=itime+1
if (itime<=1) then
   QS=Q
endif

if (itime<=1) dtf=delt/2
if (itime==2) dtf=delt
if (itime>2) dtf=2*delt


QM=Q
call getrhs(1,rhs,Q,Q_grid,time,1,work1,work2)
Q=QS + DTF*rhs

if (itime.eq.1) time=DTF+time
if (itime.gt.1) time=.5*DTF+time


! Robert filter
IF(itime.ge.3) then
   QS = (1-2*epfilt)*QM + epfilt*(Q+QS)
endif

Q_grid=Q
do n=1,3
   call ifft3d(Q_grid(1,1,n),work1)
enddo


#endif






! compute KE, max U  
maxs(1:4)=0
do j=ny1,ny2
do i=nx1,nx2
   do n=1,3
      maxs(n)=max(maxs(n),abs(Q_grid(i,j,n)))   ! max u,v,w
   enddo
   vel = abs(Q_grid(i,j,1))/delx + abs(Q_grid(i,j,2))/dely 
   maxs(4)=max(maxs(4),vel)

   if (Q_grid(i,j,3)<.1*H0) then
      print *,'warning: h is within 10% of 0'
   endif
   if (Q_grid(i,j,3)<0) then
      print *,'error: h is negative'
   endif
enddo
enddo


end subroutine rk4  





subroutine getrhs(rkstage,rhs,Qhat,Q,time,compute_ints,work,work2)
!
! evaluate RHS of N.S. equations:   -u dot grad(u) + mu * laplacian(u)
!
! if compute_ints==1, then we also return the following integrals:
!       ints(3)           ke disspaation from diffusion
!       ints(4)           vorticity z-component
!       ints(5)           helicity
!     
!
! vor(1) = w_y - v_z
! vor(2) = u_z - w_x 
! vor(3) = v_x - u_y
!
! hel = u (w_y-v_z) + v (u_z - w_x)  + w (v_x - u_y)
!
use params
use fft_interface
use spectrum
use sforcing
implicit none

! input
integer :: compute_ints,rkstage
real*8 Q(nx,ny,n_var)
real*8 Qhat(nx,ny,n_var)
real*8 time


! output
real*8 rhs(nx,ny,n_var)

!work
real*8 work(nx,ny)
real*8 work2(nx,ny)

!local
real*8 :: gradu(nx,ny,2)
real*8 :: gradv(nx,ny,2)
real*8 :: gradh(nx,ny,2)
real*8 :: divtau(nx,ny,2)
real*8 dummy,tmx1,tmx2,f_diss,fxx_diss
real*8 :: pe,ke,ke_diss,a_diss,ke_diss2,vor,gradu_diss,normdx,smag_diss
real*8,save :: f_diss_ave
integer n,i,j,k
integer im,jm,numk
real*8 XFAC,hx,hy,normS,hyper_scale,ke1

call wallclock(tmx1)

a_diss=0
ke=0
pe=0
normdx=0
ke_diss=0
vor=0
smag_diss=0

! compute stochastic forcing function which will be used
! for all RK stages in this time step
if (rkstage==1) then
   f_diss_ave=0
   ! compute new forcing function for stochastic,
   ! white in time forcing.  computed at beginning of each RK4 stage
   if (forcing_type==2 .or. forcing_type==4) then
      call sforcing_random12(rhs,Qhat,f_diss,fxx_diss,1)  
   else if (forcing_type==8) then
      ! trashes rhs - used as work array
      call stochastic_highwaveno(rhs,Qhat,f_diss,fxx_diss,1)  
   endif
endif

! compute grad(h)
do n=1,2
   call der(Q(1,1,1),gradu(1,1,n),dummy,work,DX_ONLY,n)
   call der(Q(1,1,2),gradv(1,1,n),dummy,work,DX_ONLY,n)
   call der(Q(1,1,3),gradh(1,1,n),dummy,work,DX_ONLY,n)
enddo


! u dot grad(h)
do j=ny1,ny2
do i=nx1,nx2
   rhs(i,j,3)= -(Q(i,j,1)*gradh(i,j,1) + Q(i,j,2)*gradh(i,j,2))
enddo
enddo


! coriolis
do j=ny1,ny2
do i=nx1,nx2
   rhs(i,j,1)=   fcor*Q(i,j,2)
   rhs(i,j,2)= - fcor*Q(i,j,1)
enddo
enddo



! advection 
do j=ny1,ny2
do i=nx1,nx2
   pe=pe+.5*grav*Q(i,j,3)**2 - .5*grav*H0**2
   do n=1,2
      ke = ke + .5*Q(i,j,3)*Q(i,j,n)**2
      
      normdx=normdx+Q(i,j,3)*(gradu(i,j,n)**2 + gradv(i,j,n)**2)

      rhs(i,j,1) = rhs(i,j,1) - Q(i,j,n)*gradu(i,j,n)
      rhs(i,j,2) = rhs(i,j,2) - Q(i,j,n)*gradv(i,j,n)
   enddo

   rhs(i,j,3) = rhs(i,j,3) - Q(i,j,3)*(gradu(i,j,1)+gradv(i,j,2))

   vor=vor + gradu(i,j,2)-gradv(i,j,1)

enddo
enddo



! apply forcing:
if (forcing_type>0) then
   call sforce(rhs,Qhat,f_diss,fxx_diss)
   ! average over all 4 stages
   f_diss_ave=((rkstage-1)*f_diss_ave+f_diss)/rkstage 
   ! this is not computed correcly, so set to zero
   f_diss_ave=0
   ! we need to compute < u h , f >
endif


if (alpha_value>0) then
   call compute_divtau(Q,divtau,gradu,gradv,gradh,work,work2)

   do j=ny1,ny2
   do i=nx1,nx2
   do n=1,2
      rhs(i,j,n)=rhs(i,j,n)+divtau(i,j,n)
      ! a_diss = uH dot (divtau + grav gradh)
      a_diss=a_diss + Q(i,j,n)*Q(i,j,3)*(divtau(i,j,n) + grav*gradh(i,j,n))
   enddo
   enddo
   enddo

   if (compute_transfer .and. compute_ints==1) then
      spec_model=0
      do n=1,2
         do j=ny1,ny2   
         do i=nx1,nx2
            divtau(i,j,n)=Q(i,j,3)*(divtau(i,j,n) + grav*gradh(i,j,n))
         enddo
         enddo
         call fft3d(divtau(1,1,n),work)
         call compute_spectrum_fft(Qhat(1,1,n),divtau(1,1,n),io_pe,spec_tmp)
         spec_model=spec_model + spec_tmp
      enddo
   endif
   
else
   ! g grad(h) 
   do j=ny1,ny2
   do i=nx1,nx2
      rhs(i,j,1)=rhs(i,j,1)-grav*gradh(i,j,1) 
      rhs(i,j,2)=rhs(i,j,2)-grav*gradh(i,j,2) 
   enddo
   enddo
endif

! smagorinsky terms
!
if (smagorinsky>0) then
   ! divtau used for smag. term
   ! gradh used for work storage
   call compute_smagorinsky(divtau,Q,gradu,gradv,gradh,work,work2,smag_diss)
   rhs(:,:,1:2)=rhs(:,:,1:2)+divtau(:,:,1:2)


   if (compute_transfer .and. compute_ints==1) then
      ! E dissipation from smagorinsky (or alpha term):
      !       h*U dot modeling_term
      ! which we write as:
      !       U dot (h*modeling_term)
      ! and fourier transform each of those to compute the 
      ! spectrum for smagorinsky
      !
      spec_model=0
      do n=1,2
         divtau(:,:,n)=divtau(:,:,n)*Q(:,:,3)
         call fft3d(divtau(1,1,n),work)
         call compute_spectrum_fft(Qhat(1,1,n),divtau(1,1,n),io_pe,spec_tmp)
         spec_model=spec_model + spec_tmp
      enddo
   endif
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  gradh, gradu, gradv have now been used for work storage
!  and no longer have gradient information below this point
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! add in diffusion term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (mu_hyper>0) then
   call ke_shell(Qhat,ke1,numk,dealias_23_kmax2_1,dealias_23_kmax2)
   ! units of E = m**2/s**2
   ! units of E(k) = m**3/s**2
   ! hyper viscosity:  E(kmax)* k**8 / kmax**(8-1.5)  
   ! scaling:  E(kmax)/(kmax**2)(4-.75)
   hyper_scale = sqrt(ke1) * (pi2_squared*2*dealias_23_kmax2)**(-(mu_hyper-.75))

   ! du/dt = (sqrt(E(kmax))  [1/ (kmax**8 kmax**alpha) ] k**8  u  
   ! m/s**2  =  m**1.5/s  kmax**-alpha  m/s
   !  1 = m**1.5 kmax**-alpha     = m**(1.5+alpha)   alpha = -1.5
endif



! use gradu to store -del**4 U to compute viscous KE dissapation
! use gradv to store del U to compute viscous KE dissapation (for alpha)
gradu=0
do j=ny1,ny2
   jm=abs(jmcord(j))
   do i=nx1,nx2
      im=abs(imcord(i))
      
      ! laplacian  del U
      xfac=-((im*im + jm*jm )*pi2_squared)
      gradv(i,j,1)=xfac*Qhat(i,j,1)
      gradv(i,j,2)=xfac*Qhat(i,j,2)
      
      ! -laplacian**4 
      xfac=-(xfac**mu_hyper)
      gradu(i,j,1)=xfac*Qhat(i,j,1)
      gradu(i,j,2)=xfac*Qhat(i,j,2)
      
   enddo
enddo

! back to gridspace:
call ifft3d(gradu(1,1,1),work)
call ifft3d(gradu(1,1,2),work)

do j=ny1,ny2
   do i=nx1,nx2
      rhs(i,j,1)=rhs(i,j,1)+mu_hyper_value*hyper_scale*gradu(i,j,1)/Q(i,j,3)
      rhs(i,j,2)=rhs(i,j,2)+mu_hyper_value*hyper_scale*gradu(i,j,2)/Q(i,j,3)
   enddo
enddo


! compute ke_diss = energy disspation from hyperviscosity:
if (compute_ints==1) then
   call ifft3d(gradv(1,1,1),work)
   call ifft3d(gradv(1,1,2),work)
   do j=ny1,ny2
   do i=nx1,nx2
      ke_diss = ke_diss + (Q(i,j,1)*gradu(i,j,1) &
                        + Q(i,j,2)*gradu(i,j,2)) 
      gradu_diss = gradu_diss + (gradu(i,j,1)*gradv(i,j,1) &
                        + gradu(i,j,2)*gradv(i,j,2)) 
   enddo
   enddo

   ints(1)=gradu_diss/g_nx/g_ny
   ints(2)=normdx/g_nx/g_ny
   !ints(3) = 
   ints(4)=vor/g_nx/g_ny
   ints(5)=ke/g_nx/g_ny
   ints(6)=(ke+pe)/g_nx/g_ny
   !ints(7)
   ints(8)=a_diss/g_nx/g_ny         ! < Hu,div(tau)' >  
   ! ints(9)  = 
   ints(10)=(smag_diss+mu_hyper_value*hyper_scale*ke_diss)/g_nx/g_ny     ! u dot laplacian u



   if (compute_transfer) then
      transfer_comp_time=time
      spec_diff=0
      do n=1,2
         do j=ny1,ny2
            jm=abs(jmcord(j))
            do i=nx1,nx2
               im=abs(imcord(i))
               ! laplacian  del U
               xfac=((im*im + jm*jm )*pi2_squared)
               xfac=(xfac**4)
               work(i,j)=sqrt(xfac)*Qhat(i,j,n)
            enddo
         enddo
         call compute_spectrum_fft(work,work,io_pe,spec_tmp)
         spec_diff=spec_diff - mu_hyper_value*hyper_scale*spec_tmp
      enddo

      

   endif


if (forcing_type==8) then
   ! in this case, f_diss from stage 1 is not very accurate, so use
   ! the average.  Reason: at small scales, forcing has a huge 
   ! effect.  stage1: u & f uncorrelated, <u,f>=small.  But
   ! after a few stages, u & f very correlated, <u,v>=large.  
   ints(3)=f_diss_ave
endif

endif



! back to spectral space:
do n=1,3
   call fft3d(rhs(1,1,n),work)
enddo



! dealias
do j=ny1,ny2
   jm=abs(jmcord(j))
   do i=nx1,nx2
      im=abs(imcord(i))
      
      if ( (jm>g_ny/3)  .or. (im>g_nx/3) )  then
         rhs(i,j,1)=0
         rhs(i,j,2)=0
         rhs(i,j,3)=0
      endif
   enddo
enddo



call wallclock(tmx2)
tims(5)=tims(5)+(tmx2-tmx1)
end








subroutine compute_smagorinsky(rhs,Q,gradu,gradv,s1,work,work2,smag_diss)
!
!  add smagorinsky term into rhs
!
!   S= |       u_x       .5(u_y + v_x)  |
!      |   .5(u_y+v_x)      v_y         |
!
!
! diffusion:    1/rho   div ( mu S)  = 1/rho   div ( rho nu S)
! mu = rho*nu    mu is the quantity that is usually constant?
! 
! incompressible Smag:           div (c1 |S| S )
! compressible Smag:      (1/H)  div (c1 H |S|) S
! 
use params
use sforcing
use transpose
implicit none
real*8 rhs(nx,ny,2)         
real*8 Q(nx,ny,n_var)         
real*8 gradu(nx,ny,2)
real*8 gradv(nx,ny,2)
real*8 s1(nx,ny,2)
real*8 work(nx,ny)
real*8 work2(nx,ny)
real*8 smag_diss

real*8 :: normS,dummy
integer :: i,j


rhs=0
smag_diss=0
normS=0
   do j=ny1,ny2
   do i=nx1,nx2
      normS=normS+gradu(i,j,1)**2 + gradv(i,j,2)**2 + &
         .5*(gradu(i,j,2)+gradv(i,j,1))**2

      s1(i,j,1)=gradu(i,j,1)
      s1(i,j,2)=.5*(gradu(i,j,2)+gradv(i,j,1))
   enddo
   enddo
   normS=sqrt(2*normS)  ! following Marcel Lesieur, Turbulence in Fluids

   s1(:,:,1)=2*(smagorinsky**2)*normS*s1(:,:,1)*Q(:,:,3)
   s1(:,:,2)=2*(smagorinsky**2)*normS*s1(:,:,2)*Q(:,:,3)

   call der(s1(1,1,1),work,dummy,work2,DX_ONLY,1)
   work=work*(delx**2)/Q(:,:,3)
   do j=ny1,ny2
   do i=nx1,nx2
      rhs(i,j,1)=rhs(i,j,1)+work(i,j)      
      smag_diss=smag_diss+Q(i,j,3)*Q(i,j,1)*work(i,j)
   enddo
   enddo

   call der(s1(1,1,2),work,dummy,work2,DX_ONLY,2)
   work=work*(dely**2)/Q(:,:,3)
   do j=ny1,ny2
   do i=nx1,nx2
      rhs(i,j,1)=rhs(i,j,1)+work(i,j)      
      smag_diss=smag_diss+Q(i,j,3)*Q(i,j,1)*work(i,j)
   enddo
   enddo



   do j=ny1,ny2
   do i=nx1,nx2
      ! use s1 to store n=2 row of S:
      s1(i,j,1)=.5*(gradu(i,j,2)+gradv(i,j,1))
      s1(i,j,2)=gradv(i,j,2)
   enddo
   enddo

   s1(:,:,1)=2*(smagorinsky**2)*normS*s1(:,:,1)*Q(:,:,3)
   s1(:,:,2)=2*(smagorinsky**2)*normS*s1(:,:,2)*Q(:,:,3)

   call der(s1(1,1,1),work,dummy,work2,DX_ONLY,1)
   work=work*(delx**2)/Q(:,:,3)
   do j=ny1,ny2
   do i=nx1,nx2
      rhs(i,j,2)=rhs(i,j,2)+work(i,j)      
      smag_diss=smag_diss+Q(i,j,3)*Q(i,j,2)*work(i,j)
   enddo
   enddo

   call der(s1(1,1,2),work,dummy,work2,DX_ONLY,2)
   work=work*(dely**2)/Q(:,:,3)
   do j=ny1,ny2
   do i=nx1,nx2
      rhs(i,j,2)=rhs(i,j,2)+work(i,j)      
      smag_diss=smag_diss+Q(i,j,3)*Q(i,j,2)*work(i,j)
   enddo
   enddo
end subroutine





subroutine compute_divtau(Q,div,gradu,gradv,gradh,work,work2)
! compute one entry of work=DD + DD' - D'D  
! transform work -> p
! compute d(p), apply helmholtz inverse, accumualte into rhs
!
!
use params
use sforcing
use transpose
implicit none
real*8 Q(nx,ny,n_var)         
real*8 div(nx,ny,2)
real*8 gradu(nx,ny,2)
real*8 gradv(nx,ny,2)
real*8 gradh(nx,ny,2)
real*8 work(nx,ny)
real*8 :: work2(nx,ny)
external helmholtz_hform_periodic

!local
integer i,j,k,m1,m2,im,jm,km,l,n
integer nd,ifilt,n1,n1d,n2,n2d,n3,n3d
real*8 :: D(2,2),dummy,xfac


div=0
do m2=1,2
do m1=1,2

   ! compute Tau(m1,m2)   Tau = DD + DD'  - D'D
   do j=ny1,ny2
   do i=nx1,nx2
      do nd=1,2
         D(1,nd) = gradu(i,j,nd)
         D(2,nd) = gradv(i,j,nd)
      enddo
      work(i,j)=0
      do L=1,2
         work(i,j)=work(i,j)+D(m1,L)*D(L,m2)+D(m1,L)*D(m2,L)-D(L,m1)*D(L,m2)
! leray-alpha:
!         work(i,j)=work(i,j)+D(m1,L)*D(L,m2)+D(m1,L)*D(m2,L)
      enddo
   enddo
   enddo

   do j=ny1,ny2
   do i=nx1,nx2
      work(i,j)=work(i,j)*Q(i,j,3)
   enddo
   enddo

   ! compute div(Tau)  
   call der(work,work,dummy,work2,DX_ONLY,m2)
   div(nx1:nx2,ny1:ny2,m1) = div(nx1:nx2,ny1:ny2,m1) -&
                           alpha_value**2 * work(nx1:nx2,ny1:ny2)



enddo
enddo


do m1=1,2
   do j=ny1,ny2
   do i=nx1,nx2
      div(i,j,m1)=div(i,j,m1)/Q(i,j,3)
   enddo
   enddo
enddo



   ! Apply Helmholtz inverse to:   div(tau) - grav grad(h)
   ! 
   ! alpha model:
   !     u_t + ... = helm^-1[div(tau)-grav grad(h)] 
   ! OR:
   !     u_t + ... + grav grad(h) = helm^-1[div(tau)-grav grad(h)] + grav grad(h)
   ! THUS:
   ! 
   ! a_diss should be the KE dissapation from the div(tau) term,
   ! a_diss = <uH, Helmholtz^-1(div(tau)-grav grad(h)>  + <u H , grav grad(h)>
   !        = <Helm(uH),div(tau)-grav grad(h)>  + < u H, grav grad(h)>
   !        = <Helm(uH),div(tau)> -<Helm(uh),grav grad(h)>  + < u H, grav grad(h)>
   !        = <Helm(uH),div(tau)> -<Helm(uh),grav grad(h)>  + < u H, grav grad(h)>
   !        = <Helm(uH),div(tau)> + alpha_value**2 <Lap(uh),grav grad(h)>  
   ! which is rather complicated and not computed yet.  

   do n=1,2
      work=div(:,:,n)-grav*gradh(:,:,n)      
      !call fft3d(work,work2)
      !call fft_filter_dealias(work)
      !call ifft3d(work,work2)
      div(:,:,n)=work  ! use RHS as our initial guess also
      work=work*Q(:,:,3)

      call cgsolver(div(1,1,n),work,1d0,-alpha_value**2,1d-10,Q(1,1,3),&
        helmholtz_hform_periodic,.true.)
!      call cgsolver(divtau(1,1,n),work,1d0,-alpha_value**2,1d-10,Q(1,1,3),&
!        helmholtz_hform_periodic,.false.)
      !call jacobi(divtau(1,1,n),work,1d0,-alpha_value**2,1d-10,Q(1,1,3),&
      !  helmholtz_hform_periodic,.true.)

   enddo


end subroutine

