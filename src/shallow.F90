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
      if (dealias) call fft_filter_dealias(Q(1,1,n))
   enddo
   if (equations/=SHALLOW) then
      call print_message("Error: shallow water model can only run equations=SHALLOW")
      call abort("initial conditions are probably incorrect.")
   endif
   if (ndim/=2) then
      call abort("Error: shallow water model cannot run in 3D")
   endif
   if (nz/=1) then
      call abort("Error: shallow water model cannot run in 3D")
   endif
endif


#ifndef USE_LEAPFROG

Q_old=Q

! stage 1
call getrhs(rhs,Q,Q_grid,time,1,work1,work2)
Q=Q+delt*rhs/6.0

! stage 2
Q_tmp = Q_old + delt*rhs/2.0
Q_grid=Q_tmp
do n=1,3
   call ifft3d(Q_grid(1,1,n),work1)
enddo
call getrhs(rhs,Q_tmp,Q_grid,time+delt/2.0,0,work1,work2)
Q=Q+delt*rhs/3.0

! stage 3
Q_tmp = Q_old + delt*rhs/2.0
Q_grid=Q_tmp
do n=1,3
   call ifft3d(Q_grid(1,1,n),work1)
enddo
call getrhs(rhs,Q_tmp,Q_grid,time+delt/2.0,0,work1,work2)
Q=Q+delt*rhs/3.0


! stage 4
Q_tmp = Q_old + delt*rhs
Q_grid=Q_tmp
do n=1,3
   call ifft3d(Q_grid(1,1,n),work1)
enddo
call getrhs(rhs,Q_tmp,Q_grid,time+delt,0,work1,work2)
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
call getrhs(rhs,Q,Q_grid,time,1,work1,work2)
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
enddo
enddo


end subroutine rk4  





subroutine getrhs(rhs,Qhat,Q,time,compute_ints,work,work2)
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
implicit none

! input
real*8 Q(nx,ny,n_var)
real*8 Qhat(nx,ny,n_var)
real*8 time
integer compute_ints

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
real*8 dummy,tmx1,tmx2
real*8 :: pe,ke,ke_diss,a_diss,ke_diss2,vor,gradu_diss,normdx
integer n,i,j,k
integer im,jm
real*8 XFAC,hx,hy
external :: helmholtz_hform_periodic


call wallclock(tmx1)

a_diss=0
ke=0
pe=0
normdx=0
ke_diss=0
vor=0



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




! advection and viscous terms
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




if (alpha_value>0) then
   call alpha_model_forcing(Q,divtau,gradu,gradv,work,work2)
   ! Apply Helmholtz inverse to:   div(tau) - grav grad(h)

   ! a_diss should be the KE dissapation from the div(tau) term,
   ! a_diss = <uH, Helmholtz^-1(div(tau)-grav grad(h)>  - <u H , grav grad(h)>
   !        = <Helm(uH),div(tau)-grav grad(h)>  - < u H, grav grad(h)>
   !        = ...
   ! which is rather complicated and not computed yet.  

   do n=1,2
      work=divtau(:,:,n)-grav*gradh(:,:,n)      
      !call fft3d(work,work2)
      !call fft_filter_dealias(work)
      !call ifft3d(work,work2)
      divtau(:,:,n)=work  ! use RHS as our initial guess also

      work=work*Q(:,:,3)
      call cgsolver(divtau(1,1,n),work,1d0,-alpha_value**2,1d-8,Q(1,1,3),&
        helmholtz_hform_periodic,.true.)
      !call jacobi(divtau(1,1,n),work,1d0,-alpha_value**2,1d-6,Q(1,1,3),&
      !  helmholtz_hform_periodic,.true.)
   enddo

   do j=ny1,ny2
   do i=nx1,nx2
   do n=1,2
      !a_diss=a_diss+Q(i,j,3)*Q(i,j,n)*divtau(i,j,n)
      rhs(i,j,n)=rhs(i,j,n)+divtau(i,j,n)
   enddo
   enddo
   enddo

else
   ! g grad(h) 
   do j=ny1,ny2
   do i=nx1,nx2
      rhs(i,j,1)=rhs(i,j,1)-grav*gradh(i,j,1) 
      rhs(i,j,2)=rhs(i,j,2)-grav*gradh(i,j,2) 
   enddo
   enddo
endif


! back to spectral space:
do n=1,3
   call fft3d(rhs(1,1,n),work)
enddo

! hyper viscosity and dealias:
!
! use gradu to store -del**4 U to compute viscous KE dissapation
! use gradv to store del U to compute viscous KE dissapation (for alpha)
gradu=0
do j=ny1,ny2
   jm=abs(jmcord(j))
   do i=nx1,nx2
      im=abs(imcord(i))
      
      if ( (jm>g_ny/3)  .or. (im>g_nx/3) )  then
         rhs(i,j,1)=0
         rhs(i,j,2)=0
         rhs(i,j,3)=0
      else
         ! laplacian  del U
         xfac=-((im*im + jm*jm )*pi2_squared)
         gradv(i,j,1)=xfac*Qhat(i,j,1)
         gradv(i,j,2)=xfac*Qhat(i,j,2)

         ! -del**4 
         xfac=-(xfac**4)
         gradu(i,j,1)=xfac*Qhat(i,j,1)
         gradu(i,j,2)=xfac*Qhat(i,j,2)

         ! - mu del**8
         xfac=mu*xfac
         rhs(i,j,1)=rhs(i,j,1) + xfac*Qhat(i,j,1)
         rhs(i,j,2)=rhs(i,j,2) + xfac*Qhat(i,j,2)
      endif
   enddo
enddo

! compute ke_diss = energy disspation from hyperviscosity:
if (compute_ints==1) then
   call ifft3d(gradu(1,1,1),work)
   call ifft3d(gradu(1,1,2),work)
   call ifft3d(gradv(1,1,1),work)
   call ifft3d(gradv(1,1,2),work)
   do j=ny1,ny2
   do i=nx1,nx2
      ke_diss = ke_diss + (Q(i,j,3)*Q(i,j,1)*gradu(i,j,1) &
                        + Q(i,j,3)*Q(i,j,2)*gradu(i,j,2)) 
      gradu_diss = gradu_diss + (Q(i,j,3)*gradu(i,j,1)*gradv(i,j,1) &
                        + Q(i,j,3)*gradu(i,j,2)*gradv(i,j,2)) 
   enddo
   enddo
endif




if (compute_ints==1) then
   ints(1)=gradu_diss/g_nx/g_ny
   ints(2)=normdx/g_nx/g_ny
   !ints(3) = forcing terms

   ints(4)=vor/g_nx/g_ny
   ints(6)=(ke+pe)/g_nx/g_ny
   ints(5)=ke/g_nx/g_ny

   !ints(7)
   ints(8)=a_diss/g_nx/g_ny         ! < Hu,div(tau)' >  
   ! ints(9)  = < u,f >  (alpha model only)
   ints(10)=ke_diss/g_nx/g_ny     ! u dot laplacian u
endif

call wallclock(tmx2)
tims(5)=tims(5)+(tmx2-tmx1)
end











subroutine alpha_model_forcing(Q,div,gradu,gradv,work,work2)
! compute one entry of work=DD + DD' - D'D  
! transform work -> p
! compute d(p), apply helmholtz inverse, accumualte into rhs
!
!
use params
use sforcing
use transpose
implicit none
real*8 Q(nx,ny,nz,n_var)         
real*8 div(nx,ny,nz,2)
real*8 gradu(nx,ny,nz,2)
real*8 gradv(nx,ny,nz,2)
real*8 work(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

!local
integer i,j,k,m1,m2,im,jm,km,l,n
integer nd,ifilt,n1,n1d,n2,n2d,n3,n3d
real*8 :: D(2,2),dummy,xfac


k=1
div=0
do m2=1,2
do m1=1,2

   ! compute Tau(m1,m2)   Tau = DD + DD'  - D'D
   do j=ny1,ny2
   do i=nx1,nx2
      do nd=1,2
         D(1,nd) = gradu(i,j,k,nd)
         D(2,nd) = gradv(i,j,k,nd)
      enddo
      work(i,j,k)=0
      do L=1,2
         work(i,j,k)=work(i,j,k) + D(m1,L)*D(L,m2)+ D(m1,L)*D(m2,L) - D(L,m1)*D(L,m2)
      enddo
   enddo
   enddo
   do j=ny1,ny2
   do i=nx1,nx2
      work(i,j,k)=work(i,j,k)*Q(i,j,k,3)
   enddo
   enddo

   ! compute div(Tau)  
   call der(work,work,dummy,work2,DX_ONLY,m2)
   div(nx1:nx2,ny1:ny2,nz1:nz2,m1) = div(nx1:nx2,ny1:ny2,nz1:nz2,m1) -&
                           alpha_value**2 * work(nx1:nx2,ny1:ny2,nz1:nz2)


enddo
enddo

do m1=1,2
   do j=ny1,ny2
   do i=nx1,nx2
      div(i,j,k,m1)=div(i,j,k,m1)/Q(i,j,k,3)
   enddo
   enddo
enddo



end subroutine

