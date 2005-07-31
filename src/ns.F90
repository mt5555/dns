#include "macros.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine to take one Runge-Kutta 4th order time step

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk4(time,Q_grid,Q,Q_tmp,Q_old,rhs,work1,work2)
use params
implicit none
real*8 :: time
real*8 :: Q_grid(nx,ny,nz,n_var)
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: Q_tmp(nx,ny,nz,n_var)
real*8 :: Q_old(nx,ny,nz,n_var)
real*8 :: rhs(nx,ny,nz,n_var)

call rk4reshape(time,Q_grid,Q,rhs,rhs,Q_tmp,Q_old,work1,work2)
end



subroutine rk4reshape(time,Q_grid,Q,rhs,rhsg,Q_tmp,Q_old,work,work2)
use params
implicit none
real*8 :: time
real*8 :: Q_grid(nx,ny,nz,n_var)
real*8 :: Q(g_nz2,nslabx,ny_2dz,n_var)
real*8 :: Q_tmp(g_nz2,nslabx,ny_2dz,n_var)
real*8 :: Q_old(g_nz2,nslabx,ny_2dz,n_var)
real*8 :: work(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

! overlapped in memory:
real*8 :: rhs(g_nz2,nslabx,ny_2dz,n_var)
real*8 :: rhsg(nx,ny,nz,n_var)


! local variables
real*8 :: ke_old,time_old,vel
integer i,j,k,n,ierr
integer n1,n1d,n2,n2d,n3,n3d,im,jm,km
logical,save :: firstcall=.true.





if (firstcall) then
   firstcall=.false.
   if (smagorinsky>0) then
      call abort("Error: ns3dspectral model does not yet support Smagorinsky")
   endif

   if (dealias==0) then
      call abort("Error: using ns3dspectral model, which must be run dealiased")
   endif
   if (numerical_method/=FOURIER) then
      call abort("Error: ns3dspectral model requires FFT method.")
   endif
   if (equations/=NS_UVW) then
      call print_message("Error: ns3dspectral model can only runs equations==NS_UVW")
      call abort("initial conditions are probably incorrect.")
   endif
#ifdef ALPHA_MODEL
   if (Lz/=1) then
      call abort("Error: NS alpha model with aspect ratio not implemented")
   endif
#else
   if (alpha_value/=0) then
      call abort("Error: alpha>0 but this is not the alpha model!")
   endif
#endif  

   do n=1,n_var
      rhsg(:,:,:,1)=Q_grid(:,:,:,n)
      call z_fft3d_trashinput(rhsg,Q(1,1,1,n),rhsg(1,1,1,2)) ! use rhs as work array
   enddo
endif

! stage 1
call ns3D(rhs,rhs,Q,Q_grid,time,1,work,work2,1)

do n=1,n_var
   do j=1,ny_2dz
   do i=1,nslabx
   do k=1,g_nz
      Q_old(k,i,j,n)=Q(k,i,j,n)
      Q(k,i,j,n)=Q(k,i,j,n)+delt*rhs(k,i,j,n)/6.0
      Q_tmp(k,i,j,n)=Q_old(k,i,j,n) + delt*rhs(k,i,j,n)/2.0
   enddo
   enddo
   enddo
   call z_ifft3d(Q_tmp(1,1,1,n),Q_grid(1,1,1,n),work)
enddo






! stage 2
call ns3D(rhs,rhs,Q_tmp,Q_grid,time+delt/2.0,0,work,work2,2)

do n=1,n_var
   do j=1,ny_2dz
   do i=1,nslabx
   do k=1,g_nz
      Q(k,i,j,n)=Q(k,i,j,n)+delt*rhs(k,i,j,n)/3
      Q_tmp(k,i,j,n)=Q_old(k,i,j,n) + delt*rhs(k,i,j,n)/2
   enddo
   enddo
   enddo


   call z_ifft3d(Q_tmp(1,1,1,n),Q_grid(1,1,1,n),work)
enddo

! stage 3
call ns3D(rhs,rhs,Q_tmp,Q_grid,time+delt/2,0,work,work2,3)

do n=1,n_var
   do j=1,ny_2dz
   do i=1,nslabx
   do k=1,g_nz
      Q(k,i,j,n)=Q(k,i,j,n)+delt*rhs(k,i,j,n)/3
      Q_tmp(k,i,j,n)=Q_old(k,i,j,n)+delt*rhs(k,i,j,n)
   enddo
   enddo
   enddo

   call z_ifft3d(Q_tmp(1,1,1,n),Q_grid(1,1,1,n),work)
enddo

! stage 4
call ns3D(rhs,rhs,Q_tmp,Q_grid,time+delt,0,work,work2,4)


do n=1,n_var
   do j=1,ny_2dz
   do i=1,nslabx
   do k=1,g_nz
      Q(k,i,j,n)=Q(k,i,j,n)+delt*rhs(k,i,j,n)/6
   enddo
   enddo
   enddo
   call z_ifft3d(Q(1,1,1,n),Q_grid(1,1,1,n),work)
enddo

time = time + delt

! compute max U  
maxs(1:4)=0
maxs(10:11)=-9e20
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   do n=1,ndim
      maxs(n)=max(maxs(n),abs(Q_grid(i,j,k,n)))   ! max u,v,w
   enddo
   if (npassive>0) then
      maxs(10)=max(maxs(10),Q_grid(i,j,k,np1))
      maxs(11)=max(maxs(11),-Q_grid(i,j,k,np1))
   endif
   ! used for CFL
   ! physical coodinats   w/delz'   =   w / (Lz/Nz) = w / (Lz*delz)
   vel = abs(Q_grid(i,j,k,1))/delx + abs(Q_grid(i,j,k,2))/dely + abs(Q_grid(i,j,k,3))/delz/Lz
   maxs(4)=max(maxs(4),vel)
enddo
enddo
enddo

end subroutine 







subroutine ns3d(rhs,rhsg,Qhat,Q,time,compute_ints,work,p,rkstage)
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
use sforcing
use spectrum
implicit none

! input
real*8 time
integer compute_ints,rkstage

! input, but data can be trashed if needed
real*8 Qhat(g_nz2,nslabx,ny_2dz,n_var)           ! Fourier data at time t
real*8 Q(nx,ny,nz,n_var)                         ! grid data at time t

! output  (rhsg and rhs are overlapped in memory)
real*8 rhs(g_nz2,nslabx,ny_2dz,n_var)
real*8 rhsg(nx,ny,nz,n_var)    

! work/storage
real*8 work(nx,ny,nz)
! actual dimension: nx,ny,nz, since sometimes used as work array
real*8 p(g_nz2,nslabx,ny_2dz)    
                                 

!local

real*8 xw,xw2,xfac,tmx1,tmx2,xw_viss
real*8 uu,vv,ww,dummy
integer n,i,j,k,im,km,jm,ns
integer n1,n1d,n2,n2d,n3,n3d
real*8 :: ke,uxx2ave,ux2ave,ensave,vorave,helave,maxvor,ke_diss,u2
real*8 :: p_diss(n_var),pke(n_var)
real*8 :: h_diss,ux,uy,uz,vx,vy,vz,wx,wy,wz,hyper_scale(3),ens_diss
real*8 :: f_diss=0,a_diss=0,fxx_diss=0
real*8,save :: f_diss_ave
real*8 :: vor(3)
#ifdef ALPHA_MODEL
real*8,save :: gradu(nx,ny,nz,n_var)
real*8,save :: gradv(nx,ny,nz,n_var)
real*8,save :: gradw(nx,ny,nz,n_var)
#endif

!
! NOTE: for Fourier Coefficients with mode  im=g_nx/2, this is the
! sole "cosine" mode.  Its derivative maps to the g_nx/2 sine mode
! which aliases to 0 on the grid.  So this coefficient should be
! ignored.  This subroutine assumes everything is dealiased, and so
! this mode will be removed no matter what we do with it, and so
! we just dont wory about the case when im=g_nx/2.
!
call wallclock(tmx1)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  PASSIVE SCALARS
!  compute   -u dot grad(s)
!  to compute grad(s):   from s, compute: s_x, s_y, s_z 
!                                needs 2 x-transforms (forward and back),
!                                      2 y-transforms, 
!                                      2 z-transforms
!                        from shat:   shat_x, shat_y, shat_z:
!                                     3 x-transforms (back only)
!                                     3 y-transforms 
!                                     3 z-transforms
!
!  For passive scalars of type "2", we also add (pi-u**2) to the RHS
!  pi will be added later, but u**2 has to be added here
!
!  pi = p +.5u**2, so we are really adding:   p - .5u**2 to the RHS.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do ns=np1,np2
   ! compute u dot grad(s), store (temporally) in rhsg(:,:,:,1)
   call der(Q(1,1,1,ns),work,dummy,p,DX_ONLY,1)  ! s_x
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      rhsg(i,j,k,ns)=-Q(i,j,k,1)*work(i,j,k)
      if (passive_type(ns)==2) rhsg(i,j,k,ns)=rhsg(i,j,k,ns)-Q(i,j,k,1)**2
   enddo
   enddo
   enddo
   do n=2,ndim
   call der(Q(1,1,1,ns),work,dummy,p,DX_ONLY,n)  ! s_y and s_z
   if (n==3) work=work/Lz
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      rhsg(i,j,k,ns)=rhsg(i,j,k,ns)-Q(i,j,k,n)*work(i,j,k)
      if (passive_type(ns)==2) rhsg(i,j,k,ns)=rhsg(i,j,k,ns)-Q(i,j,k,n)**2
   enddo
   enddo
   enddo
   enddo

   ! we cannot dealias here, because below we use rhsg(:,:,:,3),
   ! which coult potentially trash some of rhs(:,:,:,4)
enddo







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute  ds/da(grid) from shat:    1 call to z_ifft3d
! cost: 1 call z_ifft3d:             FFT: 1x,1y,1z  transpose: 1z,2x,2y
! compute  ds/da(grid) from s-grid:  FFT: 2a        transpose: 2a
!
! ns_vorticity1:  compute vorticity-hat, transform into rhsg.
!                 3 calls to z_ifft3d:  FFT: 3x,3y,3z  transpose: 3z,6x,6y
! 
! ns_vorticity2:  Q_grid: vx,wx,uy,wy:  FFT: 4x,4y  tranpose: 4x,4y
!                 Q_hat: uz,vz          FFT: 2x,2y,2z  tranpose 2z,4x,4y
! total:                                FFT: 6x,6y,2z  tranpose: 2z,8x,8y
!               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! for alpha model, we need all 9 terms: ux,uy,uz,vx,vy,vz,wx,wy,wz
! From Qhat: 9 calls to z_ifft3d: 9Z,18x,18y,27FFT
! From Qgrid: x: 6fft, 6x,  y: 6fft, 6y,   z: 6fft, 6z   
!
! z-der from Qhat:      3z,6x,6y,9FFT
! x,yder from Qgrid:       6x,6y,12FFT 
!

#ifdef ALPHA_MODEL
   call ns_alpha_vorticity(gradu,gradv,gradw,Q,work)
#else
   if (ncpu_z==1) then
      call ns_vorticity(rhsg,Qhat,work,p)
      ! call ns_voriticyt2(rhsg,Q,Qhat,work,p)  ! not yet coded!
   else
      call ns_vorticity(rhsg,Qhat,work,p)
   endif
#endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid space computations.  
! 
! form u x (vor+f), overwrite into Q, ifft into rhs 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
vorave=0
ensave=0
helave=0
maxvor=0

! rhs = u x vor
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
#ifdef ALPHA_MODEL
   ! rhsg = gradu
   vor(1)=gradw(i,j,k,2)-gradv(i,j,k,3)
   vor(2)=gradu(i,j,k,3)-gradw(i,j,k,1)
   vor(3)=gradv(i,j,k,1)-gradu(i,j,k,2)
#else
   vor(1)=rhsg(i,j,k,1)
   vor(2)=rhsg(i,j,k,2)
   vor(3)=rhsg(i,j,k,3)
#endif

   vorave = vorave + vor(3)
   ensave = ensave + vor(1)**2 + vor(2)**2 + vor(3)**2
   
   helave = helave + Q(i,j,k,1)*vor(1) + & 
        Q(i,j,k,2)*vor(2) + & 
        Q(i,j,k,3)*vor(3)  
   
   maxvor = max(maxvor,abs(vor(1)))
   maxvor = max(maxvor,abs(vor(2)))
   maxvor = max(maxvor,abs(vor(3)))


   ! add any rotation to vorticity before computing u cross vor
   vor(3)=vor(3) + fcor
   
   !  velocity=(u,v,w)  vorticity=(a,b,c)=(wy-vz,uz-wx,vx-uy)
   !  v*(vx-uy) - w*(uz-wx) = (v vx - v uy + w wx) - w uz
   !  w*(wy-vz) - u*(vx-uy)
   !  u*(uz-wx) - v*(wy-vz)
   uu = ( Q(i,j,k,2)*vor(3) - Q(i,j,k,3)*vor(2) )
   vv = ( Q(i,j,k,3)*vor(1) - Q(i,j,k,1)*vor(3) )
   ww = ( Q(i,j,k,1)*vor(2) - Q(i,j,k,2)*vor(1) )

   
   ! overwrite Q with the result
   Q(i,j,k,1) = uu
   Q(i,j,k,2) = vv
   Q(i,j,k,3) = ww
enddo
enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! back to spectral space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do n=1,3
   call z_fft3d_trashinput(Q(1,1,1,n),rhs(1,1,1,n),work)
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! From this point on, Q() may be used as a work array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



if (compute_ints==1 .and. compute_transfer) then
   spec_diff=0
   do n=1,3
      call compute_spectrum_z_fft(Qhat(1,1,1,n),rhs(1,1,1,n),spec_tmp)
      spec_diff=spec_diff+spec_tmp
   enddo
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! add in diffusion term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (mu_hyper>=2) then
   ! compute hyper viscosity scaling based on energy in last shell:
   call ke_shell_z(Qhat,hyper_scale,ndim)
endif


ke=0
ux2ave=0
ke_diss=0
h_diss=0
ens_diss=0
uxx2ave=0
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nslabx
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

            xw=(im*im + jm*jm + km*km/Lz/Lz)*pi2_squared
            xw_viss=mu*xw
            if (mu_hyper>=2) then
               xw2=hyper_scale(1)*(im*im*pi2_squared)**mu_hyper
               xw2=xw2+hyper_scale(2)*(jm*jm*pi2_squared)**mu_hyper
               xw2=xw2+hyper_scale(3)*(km*km*pi2_squared/(Lz*Lz))**mu_hyper
               xw_viss=xw_viss + mu_hyper_value*xw2
            endif
            if (mu_hypo==1 .and. xw>0) then
               xw_viss=xw_viss + mu_hypo_value/xw
            endif

            rhs(k,i,j,1)=rhs(k,i,j,1) - xw_viss*Qhat(k,i,j,1)
            rhs(k,i,j,2)=rhs(k,i,j,2) - xw_viss*Qhat(k,i,j,2)
            rhs(k,i,j,3)=rhs(k,i,j,3) - xw_viss*Qhat(k,i,j,3)


            if (compute_ints==1) then
! < u (uxx + uyy + uzz) > = < u-hat * (uxx-hat + uyy-hat + uzz-hat) >
!                         = < u-hat*u-hat*( im**2 + jm**2 + km**2)

               xfac = 2*2*2
               if (km==0) xfac=xfac/2
               if (jm==0) xfac=xfac/2
               if (im==0) xfac=xfac/2
               
               u2=Qhat(k,i,j,1)*Qhat(k,i,j,1) + &
                    Qhat(k,i,j,2)*Qhat(k,i,j,2) + &
                    Qhat(k,i,j,3)*Qhat(k,i,j,3)
               
               ke = ke + .5*xfac*u2
               ux2ave = ux2ave + xfac*xw*u2
               ke_diss = ke_diss + xfac*xw_viss*u2
               uxx2ave = uxx2ave + xfac*xw*xw*u2
               
               ! u_x term
               vx = - pi2*im*Qhat(k,i+z_imsign(i),j,2)
               wx = - pi2*im*Qhat(k,i+z_imsign(i),j,3)
               uy = - pi2*jm*Qhat(k,i,j+z_jmsign(j),1)
               wy = - pi2*jm*Qhat(k,i,j+z_jmsign(j),3)
               uz =  - pi2*km*Qhat(k+z_kmsign(k),i,j,1)/Lz
               vz =  - pi2*km*Qhat(k+z_kmsign(k),i,j,2)/Lz
               ! vorcity: ( (wy - vz), (uz - wx), (vx - uy) )
               ! compute 2*k^2 u vor:
               h_diss = h_diss + 2*xfac*mu*xw*&
                    (Qhat(k,i,j,1)*(wy-vz) + &
                     Qhat(k,i,j,2)*(uz-wx) + &
                     Qhat(k,i,j,3)*(vx-uy)) 
               ! compute 2*k^2 vor vor:  
               ens_diss = ens_diss + 2*xfac*mu*xw*&
                    ((wy-vz)**2 + (uz-wx)**2 + (vx-uy)**2) 
               
         endif

      enddo
   enddo
enddo

if (compute_ints==1 .and. compute_transfer) then
   spec_diff=-spec_diff
   spec_f=0
   do n=1,3
      call compute_spectrum_z_fft(Qhat(1,1,1,n),rhs(1,1,1,n),spec_tmp)
      spec_diff=spec_diff+spec_tmp
      spec_f=spec_f+spec_tmp
   enddo
   ! spec_f = (for now) spectrum of advection + diffusion termskmrswith all terms
endif




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! add in forcing to rhs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef ALPHA_MODEL
! For alpha model, we add Helmholtz inverse(div(tau))
! and overwrite Q (used for work arrays)
call alpha_model_forcing(rhs,Qhat,Q,Q,gradu,gradv,gradw,work,p,a_diss)
#endif

! apply forcing:
if (forcing_type>0) then
   
   if (rkstage==1) then
      ! initialize stochastic forcings during rkstage==1:
      f_diss_ave=0
      
      !   may trash Q (used as a work array)
      !   call sforce_rkstage1_init(Q,Qhat,f_diss,fxx_diss,1)
      
      ! compute new forcing function for stochastic,
      ! white in time forcing.  computed at beginning of each RK4 stage
      if (forcing_type==2 .or. forcing_type==4) then
         call sforcing_random12(rhs,Qhat,f_diss,fxx_diss,1)  
      else if (forcing_type==8) then
         ! trashes Q - used as work array
         call stochastic_highwaveno(Q,Qhat,f_diss,fxx_diss,1)  
      endif
   endif
   
   call sforce(rhs,Qhat,f_diss,fxx_diss)
   ! average over all 4 stages
   f_diss_ave=((rkstage-1)*f_diss_ave+f_diss)/rkstage 
endif



#if 0
rhs=0
if (forcing_type>0) call sforce(rhs,Qhat,f_diss,fxx_diss)
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nslabx
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

         ! compute the divergence
         p(k,i,j)= - im*rhs(k,i+z_imsign(i),j,1) &
              - jm*rhs(k,i,j+z_jmsign(j),2) &
              - km*rhs(k+z_kmsign(k),i,j,3)/Lz

      enddo
   enddo
enddo
call z_ifft3d(p,rhsg,work)
print *,'maxval div: ',maxval(abs(rhsg(nx1:nx2,ny1:ny2,nz1:nz2,1)))
stop
#endif




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute information for transfer spectrum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! spec_f   = spectrum of advection + diffusion terms
! spec_tmp = spectrum of advection + diffusion terms + forcing terms
! so compute the difference and store in spec_f to get the spectrum
! of just the forcing terms.  
if (compute_ints==1 .and. compute_transfer) then
   spec_f=-spec_f
   do n=1,3
      call compute_spectrum_z_fft(Qhat(1,1,1,n),rhs(1,1,1,n),spec_tmp)
      spec_f=spec_f+spec_tmp
   enddo
endif

#undef TESTGRIDSPACE
#ifdef TESTGRIDSPACE
print *,'grid space dealiasing'
! fft to grid space:
do n=1,3
   call z_ifft3d(rhs(1,1,1,n),Q(1,1,1,n),work)
enddo
call divfree_gridspace_aspect(Q,p,work,rhs)
! back to spectral space:  
do n=1,3
   call z_fft3d_trashinput(Q(1,1,1,n),rhs(1,1,1,n),work)
enddo
#else


!  make rhs div-free
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nslabx
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

         ! compute the divergence
         p(k,i,j)= - im*rhs(k,i+z_imsign(i),j,1) &
              - jm*rhs(k,i,j+z_jmsign(j),2) &
              - km*rhs(k+z_kmsign(k),i,j,3)/Lz

         ! compute laplacian inverse
         xfac= (im*im +km*km/Lz/Lz + jm*jm)
         if (xfac/=0) xfac = -1/xfac
         p(k,i,j)=xfac*p(k,i,j)
         
      enddo
   enddo
enddo

do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nslabx
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

         ! compute gradient  dp/dx
         uu= - im*p(k,i+z_imsign(i),j) 
         vv= - jm*p(k,i,j+z_jmsign(j))
         ww= - km*p(k+z_kmsign(k),i,j)/Lz
         
         rhs(k,i,j,1)=rhs(k,i,j,1) - uu 
         rhs(k,i,j,2)=rhs(k,i,j,2) - vv 
         rhs(k,i,j,3)=rhs(k,i,j,3) - ww 

         ! dealias           
         if ( dealias_remove(abs(im),abs(jm),abs(km))) then
            rhs(k,i,j,1)=0
            rhs(k,i,j,2)=0
            rhs(k,i,j,3)=0
         endif

#if 0
         if (forcing_type==5) then
         ! balu forcing: keep mode=10 fixed.
         ! in 2D: modes (0,10) and (6,8)
         if (90.25 <= im**2+jm**2 + km**2 .and. im**2+jm**2+km**2<110.25) then
            rhs(k,i,j,1)=0
            rhs(k,i,j,2)=0
            rhs(k,i,j,3)=0
         endif
         endif
#endif

      enddo
   enddo
enddo
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! passive scalars:
! dealias the RHS scalars, and add diffusion:
! for passiv_type=2, also add p computed above p = (pressure + .5*u**2 )
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pke=0
p_diss=0               
do ns=np1,np2
   ! FFT from rhsg -> rhs  
   Q(:,:,:,1)=rhsg(:,:,:,ns)
   call z_fft3d_trashinput(Q,rhs(1,1,1,ns),work)

   ! de-alias, and store in RHS(:,:,:,ns)
   do j=1,ny_2dz
      jm=z_jmcord(j)
      do i=1,nslabx
         im=z_imcord(i)
         do k=1,g_nz
            km=z_kmcord(k)
            if ( dealias_remove(abs(im),abs(jm),abs(km))) then
               rhs(k,i,j,ns)=0
            else
               xw=(im*im + jm*jm + km*km/Lz/Lz)*pi2_squared
               xw_viss=xw*mu/schmidt(ns)
               rhs(k,i,j,ns)=rhs(k,i,j,ns) - xw_viss*Qhat(k,i,j,ns)
               if (passive_type(ns)==2) rhs(k,i,j,ns)=rhs(k,i,j,ns)+p(k,i,j)
            endif
#if 0
            if (compute_ints==1) then
               xfac = 2*2*2
               if (km==0) xfac=xfac/2
               if (jm==0) xfac=xfac/2
               if (im==0) xfac=xfac/2
               u2=Qhat(k,i,j,ns)*Qhat(k,i,j,ns)
               pke(ns) = pke(ns) + .5*xfac*u2
               p_diss(ns) = p_diss(ns) + xfac*xw_viss*u2
         endif
#endif
         enddo
      enddo
   enddo
enddo



if (compute_ints==1 .and. compute_transfer) then
   spec_rhs=0
   do n=1,3
      call compute_spectrum_z_fft(Qhat(1,1,1,n),rhs(1,1,1,n),spec_tmp)
      spec_rhs=spec_rhs+spec_tmp
   enddo
   transfer_comp_time=time
endif



if (compute_ints==1) then
   ! note: dont noramlize quantities computed in spectral space,
   ! but normalize quantities computed in grid space by /g_nx/g_ny/g_nz
   ints(1)=uxx2ave
   ints(2)=ux2ave                     ! <u_x,u_x>
   ints(3)=f_diss    
   ints(4)=vorave/g_nx/g_ny/g_nz
   ints(5)=helave/g_nx/g_ny/g_nz
   ints(6)=ke        
   ints(7)=ensave/g_nx/g_ny/g_nz
   ints(8)=a_diss    
   ints(9)=fxx_diss                     ! < u_xx,f>
   ints(10)=-ke_diss                 ! <u,u_xx>
   ints(11)=h_diss
   ints(12)=ens_diss
   maxs(5)=maxvor

!   ints(6)=pke(4)
!   ints(10)=-p_diss(4)
endif

if (forcing_type==8) then
   ! in this case, f_diss from stage 1 is not very accurate, so use
   ! the average.  Reason: at small scales, forcing has a huge 
   ! effect.  stage1: u & f uncorrelated, <u,f>=small.  But
   ! after a few stages, u & f very correlated, <u,v>=large.  
   ints(3)=f_diss_ave
endif

call wallclock(tmx2)
tims(5)=tims(5)+(tmx2-tmx1)
end











subroutine ns_vorticity(rhsg,Qhat,work,p)
use params
implicit none
real*8 Qhat(g_nz2,nslabx,ny_2dz,n_var)           ! Fourier data at time t
real*8 rhsg(nx,ny,nz,n_var)    
real*8 work(nx,ny,nz)
real*8 p(g_nz2,nslabx,ny_2dz)    

!local
integer n,i,j,k,im,km,jm
real*8 ux,uy,uz,wx,wy,wz,vx,vy,vz,uu,vv,ww



do n=1,3
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nslabx
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

         if (n==1) then
            wy = - jm*Qhat(k,i,j+z_jmsign(j),3)
            vz =  - km*Qhat(k+z_kmsign(k),i,j,2)/Lz
            !rhs(k,i,j,1) = pi2*(wy - vz)
            p(k,i,j) = pi2*(wy - vz)
         endif
         
         if (n==2) then
            wx = - im*Qhat(k,i+z_imsign(i),j,3)
            uz =  - km*Qhat(k+z_kmsign(k),i,j,1)/Lz
            !rhs(k,i,j,2) = pi2*(uz - wx)
            p(k,i,j) = pi2*(uz - wx)
         endif
         
         if (n==3) then
            uy = - jm*Qhat(k,i,j+z_jmsign(j),1)
            vx = - im*Qhat(k,i+z_imsign(i),j,2)
            !rhs(k,i,j,3) = pi2*(vx - uy)
            p(k,i,j) = pi2*(vx - uy)
         endif

      enddo
   enddo
enddo
call z_ifft3d(p,rhsg(1,1,1,n),work)
enddo
end subroutine







subroutine ns_alpha_vorticity(gradu,gradv,gradw,Q,work)
use params
implicit none
real*8 Q(nx,ny,nz,n_var)
real*8 gradu(nx,ny,nz,n_var)
real*8 gradv(nx,ny,nz,n_var)
real*8 gradw(nx,ny,nz,n_var)
real*8 work(nx,ny,nz)

!local
integer nd
real*8 dummy(1)


do nd=1,3
call der(Q(1,1,1,1),gradu(1,1,1,nd),dummy,work,DX_ONLY,nd)
call der(Q(1,1,1,2),gradv(1,1,1,nd),dummy,work,DX_ONLY,nd)
call der(Q(1,1,1,3),gradw(1,1,1,nd),dummy,work,DX_ONLY,nd)
enddo


end subroutine
!
!  
!
!  tau = DD + DD' - D'D
!  NS(u) = -  (1-a**2 Laplacian)^-1 a**2 div(tau) 
!
!
! 3 ways to do this:
! div on grid, inverse in spectral
!          9 der() + 3 z_fft()      transposes: 6+6 x, 6+6 y, 6+3z
! div & inverse in spectral
!         9 z_fft()                 tranposes: 18x,18y,9z
!
! x & y components of div on grid, then transform.  transform z, add then div
! 6 der() + 6 z_fft()               transposes: 6+12 x, 6+12y, 6z
!
subroutine alpha_model_forcing(rhs,Qhat,div,divs,gradu,gradv,gradw,work,p,a_diss)
! compute one entry of work=DD + DD' - D'D  
! transform work -> p
! compute d(p), apply helmholtz inverse, accumualte into rhs
!
! f_diss  KE dissapation from forcing used in RHS  
! a_diss  KE dissapation from alpha term added to RHS 
! normuf  || u dot f ||     
!
use params
use sforcing
implicit none
real*8 rhs(g_nz2,nslabx,ny_2dz,n_var)           ! Fourier data at time t
real*8 Qhat(g_nz2,nslabx,ny_2dz,n_var)          ! Fourier data at time t
real*8 p(g_nz2,nslabx,ny_2dz)    
real*8 work(nx,ny,nz)

! overlapped in memory:
real*8 div(nx,ny,nz,n_var)
real*8 divs(g_nz2,nslabx,ny_2dz,n_var)

real*8 gradu(nx,ny,nz,n_var)
real*8 gradv(nx,ny,nz,n_var)
real*8 gradw(nx,ny,nz,n_var)
real*8 :: a_diss,f_diss,normuf

!local
integer i,j,k,m1,m2,im,jm,km,l,n
integer nd,ifilt,n1,n1d,n2,n2d,n3,n3d
real*8 :: D(3,3),dummy,xfac

! tophat filter code:
integer :: i2,j2,k2,i1,j1,k1,ii,jj,kk
real*8  :: wtot,w,r
real*8  :: gfilt(0:1000)
real*8 p2(g_nz2,nslabx,ny_2dz)    
external helmholtz_periodic



div=0
do m2=1,3
do m1=1,3

   ! compute Tau(m1,m2)   Tau = DD + DD'  - D'D
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      do nd=1,3
         D(1,nd) = gradu(i,j,k,nd)
         D(2,nd) = gradv(i,j,k,nd)
         D(3,nd) = gradw(i,j,k,nd)
      enddo
      work(i,j,k)=0
      do L=1,3
#ifdef LERAY_MODEL
         work(i,j,k)=work(i,j,k) + D(m1,L)*D(L,m2)+ D(m1,L)*D(m2,L)
#else
         work(i,j,k)=work(i,j,k) + D(m1,L)*D(L,m2)+ D(m1,L)*D(m2,L) - D(L,m1)*D(L,m2)
#endif
      enddo
   enddo
   enddo
   enddo

   ! compute div(Tau)  
   call der(work,work,dummy,p,DX_ONLY,m2)
   div(nx1:nx2,ny1:ny2,nz1:nz2,m1) = div(nx1:nx2,ny1:ny2,nz1:nz2,m1) -&
        work(nx1:nx2,ny1:ny2,nz1:nz2)


enddo
enddo
if (infinite_alpha==0) div=alpha_value**2 * div


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Apply Helmholtz inverse to div(tau), accumulate into RHS.
!
! RHS  = RHS + Helmholtz_Inverse(DIV)   (variable DIV contains: -div(tau)  )
! also return a_diss, the KE dissapation from div(tau) term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#undef TOPHAT
#undef ITER

#ifdef TOPHAT

! tophat filter:
ifilt = nint(.5*(sqrt(24.0)*alpha_value)/min(delx,dely,delz))
!           ^^^^ take an extra factor of 2?
!
ifilt=2*ifilt

!
! gaussian filter.  C exp(-.25* (x/alpha)**2 )
!
! x/alpha = 4 is a good place to cutoff filter
! ifilt = 4*alpha
ifilt = nint(4*alpha_value/min(delx,dely,delz))
ifilt=ifilt-1

print *,'gaussian ifilt=',ifilt
do i=0,ifilt
   gfilt(i)=exp(-.25* (i**2)/(alpha_value/min(delx,dely,delz))**2)
enddo


gradu=0
do n=1,3
do k=nz1,nz2
   do j=ny1,ny2
      do i=nx1,nx2


         wtot=0
         do jj=j-ifilt,j+ifilt
         do ii=i-ifilt,i+ifilt
            j1=jj
            if (j1<ny1) j1=j1+nslaby
            if (j1>ny2) j1=j1-nslaby

            i1=ii
            if (i1<nx1) i1=i1+nslabx
            if (i1>nx2) i1=i1-nslabx

            w = gfilt(abs(jj-j))*gfilt(abs(ii-i))
            !w=1

            wtot=wtot+w
            gradu(i,j,k,n)=gradu(i,j,k,n)+w*div(i1,j1,k,n)

         enddo
         enddo
         gradu(i,j,k,n)=gradu(i,j,k,n)/wtot

      enddo
   enddo
enddo
enddo
do n=1,3
   call z_fft3d_trashinput(gradu(1,1,1,n),divs(1,1,1,n),p) 
enddo


a_diss=0
do n=1,3
   do j=1,ny_2dz
      jm=z_jmcord(j)
      do i=1,nslabx
         im=z_imcord(i)
         do k=1,g_nz
            km=z_kmcord(k)
            
            xfac=8
            if (z_kmcord(k)==0) xfac=xfac/2
            if (z_jmcord(j)==0) xfac=xfac/2
            if (z_imcord(i)==0) xfac=xfac/2
            a_diss = a_diss + xfac*Qhat(k,i,j,n)*divs(k,i,j,n)

         enddo
      enddo
   enddo
enddo
rhs=rhs+divs

#elif (defined ITER)

! trash gradu and gradv (using them as work arrays below)
#define JACOBI
gradu=div
do n=1,3
   work=gradu(:,:,:,n)
#ifdef JACOBI
   call jacobi(work,gradu(1,1,1,n),1d0,-alpha_value**2,.15d0,gradv,&
     helmholtz_periodic,.false.)
#else
   call cgsolver(work,gradu(1,1,1,n),1d0,-alpha_value**2,.03d0,gradv,&
     helmholtz_periodic,.false.)
#endif
   call z_fft3d_trashinput(work,divs(1,1,1,n),p) 
enddo
rhs=rhs+divs

a_diss=0
do n=1,3
   do j=1,ny_2dz
      jm=z_jmcord(j)
      do i=1,nslabx
         im=z_imcord(i)
         do k=1,g_nz
            km=z_kmcord(k)
            
            xfac=8
            if (z_kmcord(k)==0) xfac=xfac/2
            if (z_jmcord(j)==0) xfac=xfac/2
            if (z_imcord(i)==0) xfac=xfac/2
            a_diss = a_diss + xfac*Qhat(k,i,j,n)*divs(k,i,j,n)

         enddo
      enddo
   enddo
enddo


#else

gradu=div
do n=1,3
   call z_fft3d_trashinput(gradu(1,1,1,n),divs(1,1,1,n),work) 
enddo
call helminv(rhs,divs,Qhat,a_diss)

#endif



end subroutine











subroutine helm(out,in)
use params
implicit none

real*8 out(g_nz2,nslabx,ny_2dz)
real*8 in(g_nz2,nslabx,ny_2dz)

!local
integer i,j,k,m1,m2,im,jm,km,l,n
real*8 :: xfac


do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nslabx
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)
         

         ! dealias           
         if ( ((abs(km)> g_nz/3) ) .or. &
              ((abs(jm)> g_ny/3) ) .or. &
              ((abs(im)> g_nx/3) ) )  then
            out(k,i,j)=in(k,i,j)
            stop 'fix this with new dealias function'
         else
            ! compute laplacian inverse
            xfac= -(im*im +km*km + jm*jm)*pi2_squared      
            xfac=1 - alpha_value**2 *xfac
            out(k,i,j) = xfac*in(k,i,j)
         endif
      enddo
   enddo
enddo
end subroutine 




subroutine helminv(rhs,divs,Qhat,a_diss)
use params
implicit none

real*8 rhs(g_nz2,nslabx,ny_2dz,n_var)    
real*8 Qhat(g_nz2,nslabx,ny_2dz,n_var)    
real*8 divs(g_nz2,nslabx,ny_2dz,n_var)    
real*8 a_diss

!local
integer i,j,k,m1,m2,im,jm,km,l,n
real*8 :: xfac

a_diss=0
do n=1,3
   do j=1,ny_2dz
      jm=z_jmcord(j)
      do i=1,nslabx
         im=z_imcord(i)
         do k=1,g_nz
            km=z_kmcord(k)
            
            ! compute laplacian inverse
            xfac= (im*im +km*km + jm*jm)*pi2_squared      
            if (infinite_alpha==0) then
               xfac=1 + alpha_value**2 *xfac
            endif
            if (xfac/=0) xfac = 1/xfac
            rhs(k,i,j,n) = rhs(k,i,j,n) + xfac*divs(k,i,j,n)

            xfac=8*xfac
            if (z_kmcord(k)==0) xfac=xfac/2
            if (z_jmcord(j)==0) xfac=xfac/2
            if (z_imcord(i)==0) xfac=xfac/2
            a_diss = a_diss + xfac*Qhat(k,i,j,n)*divs(k,i,j,n)

         enddo
      enddo
   enddo
enddo
end subroutine
