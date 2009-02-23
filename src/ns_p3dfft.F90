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
subroutine rk4(time,Q_grid,Q,Q_tmp,Q_old,rhs,work1,work2)
use params
implicit none
real*8 :: time
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Q_grid
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: Q_tmp(nx,ny,nz,n_var)
real*8 :: Q_old(nx,ny,nz,n_var)
real*8 :: rhs(nx,ny,nz,n_var)
real*8, allocatable, save :: q2(:,:,:,:)
logical,save :: firstcall=.true.
integer :: n

if (firstcall) then
   firstcall=.false.
   call p3_params_init
   if (smagorinsky>0) then
      call abortdns("Error: ns_p3dfft.F90 model does not yet support Smagorinsky")
   endif
#ifndef USE_P3DFFT
   call abortdns("Error: ns_p3dfft.F90 model must be compiled with -DUSE_P3DFFT")
#endif
   if (dealias==0) then
      call abortdns("Error: using ns_p3dfft.F90 model, which must be run dealiased")
   endif
   if (numerical_method/=FOURIER) then
      call abortdns("Error: ns_p3dfft.F90 model requires FFT method.")
   endif
   if (equations/=NS_UVW) then
      call abortdns("Error: ns_p3dfft.F90 model can only runs equations==NS_UVW")
   endif
   if (alpha_value/=0) then
      call abortdns("Error: alpha>0 but this is not the alpha model!")
   endif
   if (use_phaseshift) then
      if (npassive>0) call abortdns("Error: phaseshift not yet coded for passive scalars")
      allocate(q2(nx,ny,nz,ndim))
   endif
   if (n_var<3) call abortdns("Error: ns_p3dfft.F90 requires n_var>=3")
   if (mu_hyper>=2) then
      call abortdns("Error: ns_p3dfft.F90 model not yet coded for hyperviscosity")
   endif
   if (forcing_type>1) then
      call abortdns("Error: ns_p3dfft.F90 model not yet coded for forcing types > 1")
   endif

   ! intialize Q with Fourier Coefficients:
   call p3_fft3d_nvar(Q_grid,Q,work1)
endif
call rk4reshape(time,Q_grid,Q,rhs,rhs,Q_tmp,Q_old,work1,work2,q2)
end




subroutine rk4reshape(time,Q_grid,Q,rhs,rhsg,Q_tmp,Q_old,work,work2,q2)
use params
use p3dfft
implicit none
real*8 :: time
real*8 :: Q_grid(nx,ny,nz,n_var)
complex*16 :: Q(p3_n1,p3_n2,p3_n3,n_var)
complex*16 :: Q_tmp(p3_n1,p3_n2,p3_n3,n_var)
complex*16 :: Q_old(p3_n1,p3_n2,p3_n3,n_var)
real*8 :: work(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: q2(nx,ny,nz,ndim)  ! only allocated if phase shifting turned on

! overlapped in memory:  dont use both in the same n loop
complex*16 :: rhs(p3_n1,p3_n2,p3_n3,n_var)
real*8 :: rhsg(nx,ny,nz,n_var)


! local variables
real*8 :: ke_old,time_old,vel
integer i,j,k,n,ierr
integer n1,n1d,n2,n2d,n3,n3d,im,jm,km

! stage 1
call ns3D(rhs,rhsg,Q,Q_grid,time,1,work,work2,1,q2)
do n=1,n_var
   do k=1,p3_n3
   do j=1,p3_n2
   do i=1,p3_n1
      Q_old(i,j,k,n)=Q(i,j,k,n)
      Q(i,j,k,n)=Q(i,j,k,n)+delt*rhs(i,j,k,n)/6                ! accumulate RHS
      Q_tmp(i,j,k,n)=Q_old(i,j,k,n) + delt*rhs(i,j,k,n)/2      ! Euler timestep with dt=delt/2
   enddo
   enddo
   enddo
enddo
!if (hyper_implicit==1) call hyper_filter(Q_tmp,delt/2)
do n=1,n_var
   call btran_c2r(Q_tmp(1,1,1,n),Q_grid(1,1,1,n))
enddo




! stage 2
call ns3D(rhs,rhsg,Q_tmp,Q_grid,time+delt/2.0,0,work,work2,2,q2)
do n=1,n_var
   do k=1,p3_n3
   do j=1,p3_n2
   do i=1,p3_n1
      Q(i,j,k,n)=Q(i,j,k,n)+delt*rhs(i,j,k,n)/3             ! accumulate RHS
      Q_tmp(i,j,k,n)=Q_old(i,j,k,n) + delt*rhs(i,j,k,n)/2   ! Euler timestep with dt=delt/2
   enddo
   enddo
   enddo
enddo
!if (hyper_implicit==1) call hyper_filter(Q_tmp,delt/2)
do n=1,n_var
   call btran_c2r(Q_tmp(1,1,1,n),Q_grid(1,1,1,n))
enddo


! stage 3
call ns3D(rhs,rhsg,Q_tmp,Q_grid,time+delt/2,0,work,work2,3,q2)
do n=1,n_var
   do k=1,p3_n3
   do j=1,p3_n2
   do i=1,p3_n1
      Q(i,j,k,n)=Q(i,j,k,n)+delt*rhs(i,j,k,n)/3               ! accumulate RHS
      Q_tmp(i,j,k,n)=Q_old(i,j,k,n)+delt*rhs(i,j,k,n)      ! Euler timestep with dt=delt
   enddo
   enddo
   enddo
enddo
!if (hyper_implicit==1) call hyper_filter(Q_tmp,delt)
do n=1,n_var
   call btran_c2r(Q_tmp(1,1,1,n),Q_grid(1,1,1,n))
enddo



! stage 4
call ns3D(rhs,rhsg,Q_tmp,Q_grid,time+delt,0,work,work2,4,q2)
do n=1,n_var
   do k=1,p3_n3
   do j=1,p3_n2
   do i=1,p3_n1
      Q(i,j,k,n)=Q(i,j,k,n)+delt*rhs(i,j,k,n)/6        ! final Euler timestep with delt
   enddo   
   enddo
   enddo
enddo
!if (hyper_implicit==1) call hyper_filter(Q,delt)

do n=1,n_var
   call btran_c2r(Q(1,1,1,n),Q_grid(1,1,1,n))
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







subroutine ns3d(rhs,rhsg,Qhat,Q,time,compute_ints,work,p,rkstage,q2)
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
use p3_forcing
use spectrum
implicit none

! input
real*8 time
integer compute_ints,rkstage

! input, but Q can be overwritten if needed
complex*16 Qhat(p3_n1,p3_n2,p3_n3,n_var)         ! Fourier data at time t
real*8 Q(nx,ny,nz,n_var)                         ! grid data at time t

! output  (rhsg and rhs are overlapped in memory)
complex*16 rhs(p3_n1,p3_n2,p3_n3,n_var)
real*8 rhsg(nx,ny,nz,n_var)    

! work/storage
real*8  :: q2(nx,ny,nz,ndim)  ! only allocated if use_phaseshift
                                  ! phase shifted vorticity 
real*8 work(nx,ny,nz)
! actual dimension: nx,ny,nz, since sometimes used as work array
complex*16 p(p3_n1,p3_n2,p3_n3)    

                                 

!local

real*8 xw,xw2,xfac,tmx1,tmx2,xw_viss,dummy
integer n,i,j,k,im,km,jm,ns
integer n1,n1d,n2,n2d,n3,n3d
real*8 :: ke,uxx2ave,ux2ave,ensave,vorave,helave,maxvor,ke_diss,u2,ens_alpha
real*8 :: p_diss(n_var),pke(n_var)
real*8 :: h_diss,hyper_scale(n_var,n_var),ens_diss2,ens_diss4,ens_diss6
real*8 :: f_diss=0,a_diss=0,fxx_diss=0
real*8,save :: f_diss_ave
real*8 :: vor(3),xi,uu,vv,ww
complex*16 :: ptmp
complex*16  :: ux,uy,uz,vx,vy,vz,wx,wy,wz



call wallclock(tmx1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  PASSIVE SCALARS
!  compute   -u dot grad(s)
!
!  to compute grad(s):   from s, compute: s_x, s_y, s_z 
!                                needs 2 x-transforms (1 forward, 1 back),
!                                      2 y-transforms, 
!                                      2 z-transforms
!                        PLUS: the FFT to compute s from s_hat:  
!                                       1 x-transform
!                                       2 y-transform
!                                       1 z-transform
!                        TOTAL, with y-pencel reference:  
!                            6 transposes, 9 FFTs 
!
!                        "fast" way: compute some of the derivatives
!                        while computing the iFFT:   get s_x 
!                         for 1 extra ifft, no extra transposes
!                                       1 x-transform
!                                       2 y-transform
!                                       1 z-transform
!                           PLUS:  s_y, s_z:  
!                                      2 y-transforms, 
!                                      2 z-transforms
!                        TOTAL, with y-pencel reference:  
!                            4 transposes, 8 FFTs 
!
! 
!                        from s_hat:   s_hat_x, s_hat_y, s_hat_z:
!                                     3 x-transforms (no forward, 3 back)
!                                     6 y-transforms 
!                                     3 z-transforms
!                         PLUS, for BOUSS, we also need s!
!
!
!  For passive scalars of type "2", we also add (pi-u**2) to the RHS
!  pi will be added later, but u**2 has to be added here
!
!  pi = p +.5u**2, so we are really adding:   p - .5u**2 to the RHS.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do ns=np1,np2
   ! compute u dot grad(s), store (temporally) in rhsg(:,:,:,ns)
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
   ! we cannot dealias here, because below we use rhsg(:,:,:,3), below
   ! which coult potentially trash some of rhs(:,:,:,4)
enddo

! check of first scaler is the density and we are running boussinisque
if (passive_type(np1)==4) then
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      rhsg(i,j,k,np1)=rhsg(i,j,k,np1) + bous*Q(i,j,k,3)
   enddo
   enddo
   enddo
endif

call ns_vorticity(rhsg,Qhat,p)


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
   vor(1)=rhsg(i,j,k,1)
   vor(2)=rhsg(i,j,k,2)
   vor(3)=rhsg(i,j,k,3)
   vorave = vorave + vor(3)
   ensave = ensave + vor(1)**2 + vor(2)**2 + vor(3)**2
   
   helave = helave + Q(i,j,k,1)*vor(1) + & 
        Q(i,j,k,2)*vor(2) + & 
        Q(i,j,k,3)*vor(3)  
   
   maxvor = max(maxvor,maxval(abs(vor(:))))

   if (use_phaseshift) vor=vor/2   ! half contribution from this grid
                                   ! half contribution from phase shifted grid
   
   ! add any rotation to vorticity before computing u cross vor
#ifdef BETAPLANE
   vor(3)=vor(3) + ycord(j)*fcor
#else
   vor(3)=vor(3) + fcor
#endif
   
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

   if (npassive>0 .and. passive_type(np1)==4) then
      Q(i,j,k,3)=Q(i,j,k,3)-bous*Q(i,j,k,np1)
   endif

enddo
enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! back to spectral space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do n=1,3
   call ftran_r2c(Q(1,1,1,n),rhs(1,1,1,n))
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! phase shifting:  now compute u x vor with a phase shift,
! and add to RHS:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (use_phaseshift) then
   ! compute phaseshifted (u,v,w), store in Q()
   do n=1,3
      call p3_phaseshift(Qhat(1,1,1,n),1,work)  ! phaseshift Qhat
      call btran_c2r(Qhat(1,1,1,n),Q(1,1,1,n))
   enddo
   ! compute phaseshifted vorticity, store in q2:
   call ns_vorticity(q2,Qhat,p)

   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      uu = ( Q(i,j,k,2)*q2(i,j,k,3) - Q(i,j,k,3)*q2(i,j,k,2) )
      vv = ( Q(i,j,k,3)*q2(i,j,k,1) - Q(i,j,k,1)*q2(i,j,k,3) )
      ww = ( Q(i,j,k,1)*q2(i,j,k,2) - Q(i,j,k,2)*q2(i,j,k,1) )
      ! overwrite Q with the result
      Q(i,j,k,1) = uu
      Q(i,j,k,2) = vv
      Q(i,j,k,3) = ww
   enddo
   enddo
   enddo

   ! back spectral space
   do n=1,3
      call ftran_r2c(Q(1,1,1,n),p)
      call p3_phaseshift(p,-1,work)             ! un-phaseshift p
      rhs(:,:,:,n) = rhs(:,:,:,n) + p/2
      call p3_phaseshift(Qhat(1,1,1,n),-1,work) ! restore Qhat to unphaseshifted version
   enddo
endif

! scaling needed for P3DFFT
rhs=rhs/g_nx/g_ny/g_nz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! From this point on, Q() may be used as a work array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Transfer function diagnostics
!  spec_diff:  spectrum of u dot (u cross omega) 
!              (used as temporary storage for now)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if 0
if (compute_ints==1 .and. compute_transfer) then
   spec_diff=0
   spec_curl_diff=0
   do n=1,ndim
      call compute_spectrum_z_fft(Qhat(1,1,1,n),rhs(1,1,1,n),spec_tmp)
      spec_diff=spec_diff+spec_tmp
   enddo
   call compute_spectrum_curl_z_fft(Qhat,rhs,spec_tmp)
   spec_curl_diff=spec_curl_diff+spec_tmp   
endif
#endif


#if 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! add in hyperviscsoity diffusion term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (mu_hyper_value>0 .and. mu_hyper>=2 .and. hyper_implicit==0) then
   ! compute hyper viscosity scaling based on energy in last shell:
   ! print *,'calling ke_shell  Q=',(qhat(1,1,1,1:3))
   call ke_shell_z(Qhat,hyper_scale)
endif
#endif


ke=0
ux2ave=0
ke_diss=0
h_diss=0
ens_diss2=0
ens_diss4=0
ens_diss6=0
ens_alpha =0
uxx2ave=0



! note: we assume P3DFFT was compiled with STRIDE1 option.
! that means 1st dimension is really Z direction
!            2nd dimension is really X direction
!            3nd dimension is really Y direction
! X wave number im = p3_2mcord(j)
! Y wave number jm = p3_3mcord(k)
! Z wave number km = p3_1mcord(i)
do k=1,p3_n3
   jm=p3_3mcord(k)
   do j=1,p3_n2
      im=p3_2mcord(j)
      do i=1,p3_n1
         km=p3_1mcord(i)

            xw=(im*im + jm*jm + km*km/Lz/Lz)*pi2_squared
            xw_viss=mu*xw
            if (mu_hyper>=2 .and. hyper_implicit==0 ) then
               xw2=hyper_scale(1,1)*(im*im*pi2_squared)
               xw2=xw2+hyper_scale(2,1)*(jm*jm*pi2_squared)
               xw2=xw2+hyper_scale(3,1)*(km*km*pi2_squared/(Lz*Lz))
               xw_viss=xw_viss + mu_hyper_value*xw2**mu_hyper
            endif
            if (mu_hyper==0) then
               xw_viss=xw_viss + mu_hyper_value
            endif
            if (mu_hypo==1 .and. xw>0) then
               xw_viss=xw_viss + mu_hypo_value/xw
            endif

            rhs(i,j,k,1)=rhs(i,j,k,1) - xw_viss*Qhat(i,j,k,1)
            rhs(i,j,k,2)=rhs(i,j,k,2) - xw_viss*Qhat(i,j,k,2)
            rhs(i,j,k,3)=rhs(i,j,k,3) - xw_viss*Qhat(i,j,k,3)


            if (compute_ints==1) then
! < u (uxx + uyy + uzz) > = < u-hat * (uxx-hat + uyy-hat + uzz-hat) >
!                         = < u-hat*u-hat*( im**2 + jm**2 + km**2)

               xfac = 2
               if (im==0) xfac=xfac/2
               
               u2=real( &
                    Qhat(i,j,k,1)*conjg(Qhat(i,j,k,1)) + &
                    Qhat(i,j,k,2)*conjg(Qhat(i,j,k,2)) + &
                    Qhat(i,j,k,3)*conjg(Qhat(i,j,k,3))  )
               
               ke = ke + .5*xfac*u2
               ux2ave = ux2ave + xfac*xw*u2
               ke_diss = ke_diss + xfac*xw_viss*u2
               uxx2ave = uxx2ave + xfac*xw*xw*u2
               

               ! update these for complex*16
               ! u_x term
               vx =  pi2*cmplx(0,im)*Qhat(i,j,k,2)
               wx =  pi2*cmplx(0,im)*Qhat(i,j,k,3)
               uy =  pi2*cmplx(0,jm)*Qhat(i,j,k,1)
               wy =  pi2*cmplx(0,jm)*Qhat(i,j,k,3)
               uz =  pi2*cmplx(0,km)*Qhat(i,j,k,1)/Lz
               vz =  pi2*cmplx(0,km)*Qhat(i,j,k,2)/Lz
               ! vorcity: ( (wy - vz), (uz - wx), (vx - uy) )
               ! compute 2*k^2 u vor:
               h_diss = h_diss + 2*xfac*mu*xw*&
                    imag(Qhat(i,j,k,1)*conjg(wy-vz) + &
                     Qhat(i,j,k,2)*conjg(uz-wx) + &
                     Qhat(i,j,k,3)*conjg(vx-uy)) 
               ! incorrect if using hyperviscosity
               ens_diss2=ens_diss2 + 2*xfac*(mu*xw)*  &
                       (abs(wy-vz)**2 + abs(uz-wx)**2 + abs(vx-uy)**2) 
               ens_diss4=ens_diss4 + 2*xfac*(mu*xw*xw)*  &
                    (abs(wy-vz)**2 + abs(uz-wx)**2 + abs(vx-uy)**2)
               ens_diss6=ens_diss6 + 2*xfac*(mu*xw*(xw)**2)*  &
                    (abs(wy-vz)**2 + abs(uz-wx)**2 + abs(vx-uy)**2)      

             endif           

      enddo
   enddo
enddo


#if 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  spec_tmp:  spectrum of u dot (u cross omega + grad(pi) + diffusion) 
!  since spectrum of u dot grad(pi) is zero, we can now take the
!  differece between spec_tmp and spec_diff to get the diffusion
!  spectrum.  (stored in spec_diff)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (compute_ints==1 .and. compute_transfer) then
   spec_diff=-spec_diff
   spec_f=0
   spec_curl_diff=-spec_curl_diff
   spec_curl_f=0
   do n=1,ndim
      call compute_spectrum_z_fft(Qhat(1,1,1,n),rhs(1,1,1,n),spec_tmp)
      spec_diff=spec_diff+spec_tmp
      spec_f=spec_f+spec_tmp
   enddo
   call compute_spectrum_curl_z_fft(Qhat,rhs,spec_tmp)
   spec_curl_diff=spec_curl_diff+spec_tmp
   spec_curl_f=spec_curl_f+spec_tmp

   ! spec_f = (for now) spectrum of advection + diffusion terms
endif
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! add in forcing to rhs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (forcing_type==1) then
   if (rkstage==1) then
      ! initialize stochastic forcings during rkstage==1:
      f_diss_ave=0
   endif
   call p3_sforce(rhs,Qhat,f_diss,fxx_diss)
   ! average over all 4 stages
   f_diss_ave=((rkstage-1)*f_diss_ave+f_diss)/rkstage 
endif


#if 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute information for transfer spectrum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! spec_f   = spectrum of advection + diffusion terms
! spec_tmp = spectrum of advection + diffusion terms + forcing terms
! so compute the difference and store in spec_f to get the spectrum
! of just the forcing terms.  
if (compute_ints==1 .and. compute_transfer) then
   spec_f=-spec_f
   spec_curl_f=-spec_curl_f
   do n=1,ndim
      call compute_spectrum_z_fft(Qhat(1,1,1,n),rhs(1,1,1,n),spec_tmp)
      spec_f=spec_f+spec_tmp
   enddo
   call compute_spectrum_curl_z_fft(Qhat,rhs,spec_tmp)
   spec_curl_f=spec_curl_f+spec_tmp
endif
#endif

!  make rhs div-free
! note: we assume P3DFFT was compiled with STRIDE1 option.
! that means 1st dimension is really Z direction
!            2nd dimension is really X direction
!            3nd dimension is really Y direction
! X wave number im = p3_2mcord(j)
! Y wave number jm = p3_3mcord(k)
! Z wave number km = p3_1mcord(i)
do k=1,p3_n3
   jm=p3_3mcord(k)
   do j=1,p3_n2
      im=p3_2mcord(j)
      do i=1,p3_n1
         km=p3_1mcord(i)
         
         ! compute laplacian inverse of the divergence
         xfac= (im*im +km*km/Lz/Lz + jm*jm)
         if (xfac/=0) xfac = -1/xfac

         ptmp = xfac*( cmplx(0,im)*rhs(i,j,k,1) &
              + cmplx(0,jm)*rhs(i,j,k,2) &
              + cmplx(0,km)*rhs(i,j,k,3)/Lz  )

         p(i,j,k)=ptmp  ! save p, might be needed below for passive scalars

         ! add -grad(p) term to RHS
         rhs(i,j,k,1)=rhs(i,j,k,1) - cmplx(0,im)*ptmp
         rhs(i,j,k,2)=rhs(i,j,k,2) - cmplx(0,jm)*ptmp
         rhs(i,j,k,3)=rhs(i,j,k,3) - cmplx(0,km)*ptmp/Lz

         ! dealias           
         if ( dealias_remove(abs(im),abs(jm),abs(km))) then
            rhs(i,j,k,1)=0
            rhs(i,j,k,2)=0
            rhs(i,j,k,3)=0
         endif

      enddo
   enddo
enddo


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
   stop 'error!!'
   ! de-alias, and store in RHS(:,:,:,ns)
   do k=1,p3_n3
      jm=p3_3mcord(k)
      do j=1,p3_n2
         im=p3_2mcord(j)
         do i=1,p3_n1
            km=p3_1mcord(i)
            
            if ( dealias_remove(abs(im),abs(jm),abs(km))) then
               rhs(i,j,k,ns)=0
            else
               xw=(im*im + jm*jm + km*km/Lz/Lz)*pi2_squared
               xw_viss=xw*mu/schmidt(ns)
               ! add in hyper viscosity, if used:
               if (mu_hyper>=2 .and. hyper_implicit==0 ) then
	          xw2=hyper_scale(1,ns)*(im*im*pi2_squared)
                  xw2=xw2+hyper_scale(2,ns)*(jm*jm*pi2_squared)
                  xw2=xw2+hyper_scale(3,ns)*(km*km*pi2_squared/(Lz*Lz))
                  xw_viss=xw_viss + mu_hyper_value*(xw2**mu_hyper)

               endif
               if (mu_hyper==0) then
                  xw_viss=xw_viss + mu_hyper_value
               endif
               rhs(i,j,k,ns)=rhs(i,j,k,ns) - xw_viss*Qhat(i,j,k,ns)
               if (passive_type(ns)==2) rhs(i,j,k,ns)=rhs(i,j,k,ns)+p(i,j,k)
            endif
            if (compute_ints==1) then
               xfac = 2
               if (im==0) xfac=xfac/2
               u2=Qhat(i,j,k,ns)*Qhat(i,j,k,ns)
               pke(ns) = pke(ns) + .5*xfac*u2
               p_diss(ns) = p_diss(ns) + xfac*xw_viss*u2
         endif
         enddo
      enddo
   enddo
enddo


! compute spectrum of u dot RHS, store in spec_rhs
#if 0
if (compute_ints==1 .and. compute_transfer) then
   spec_rhs=0
   spec_curl_rhs=0
   do n=1,ndim
      call compute_spectrum_z_fft(Qhat(1,1,1,n),rhs(1,1,1,n),spec_tmp)
      spec_rhs=spec_rhs+spec_tmp
   enddo
   call compute_spectrum_curl_z_fft(Qhat,rhs,spec_tmp)
   spec_curl_rhs=spec_curl_rhs+spec_tmp

   transfer_comp_time=time
endif
#endif


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
   ints(12)=ens_diss2
   ints(13)=ens_diss4
   ints(14)=ens_diss6
   ints(15)=pke(np1)
   ints(16)=-p_diss(np1)

   maxs(5)=maxvor
endif
! on the 4th RK stage, overwrite ints(3) with the average over all
! stages
if (forcing_type==2 .or. forcing_type==4 .or. forcing_type==8) then
   ! in this case, f_diss from stage 1 is not very accurate, so use
   ! the average.  Reason: at small scales, forcing has a huge 
   ! effect.  stage1: u & f uncorrelated, <u,f>=small.  But
   ! after a few stages, u & f very correlated, <u,f>=large.  
   ints(3)=f_diss_ave
endif


call wallclock(tmx2)
tims(5)=tims(5)+(tmx2-tmx1)
end











subroutine ns_vorticity(rhsg,Qhat,p)
use params
implicit none
complex*16 Qhat(p3_n1,p3_n2,p3_n3,n_var)           ! Fourier data at time t
real*8 rhsg(nx,ny,nz,n_var)    
complex*16 p(p3_n1,p3_n2,p3_n3)    

!local
integer n,i,j,k
complex*16 ux,uy,uz,wx,wy,wz,vx,vy,vz
complex*16 im,km,jm

! note: we assume P3DFFT was compiled with STRIDE1 option.
! that means 1st dimension is really Z direction
!            2nd dimension is really X direction
!            3nd dimension is really Y direction
! X wave number im = p3_2mcord(j)
! Y wave number jm = p3_3mcord(k)
! Z wave number km = p3_1mcord(i)
do n=1,3
   do k=1,p3_n3
      jm=cmplx(0,p3_3mcord(k))   ! y derivative
      do j=1,p3_n2
         im=cmplx(0,p3_2mcord(j))   ! x deriviative 
         do i=1,p3_n1  
            km=cmplx(0,p3_1mcord(i))  ! z deriviative


            if (n==1) then
               wy =  jm*Qhat(i,j,k,3)
               vz =   km*Qhat(i,j,k,2)/Lz
               p(i,j,k) = pi2*(wy - vz)
            endif
            
            if (n==2) then
               wx =  im*Qhat(i,j,k,3)
               uz =   km*Qhat(i,j,k,1)/Lz
               p(i,j,k) = pi2*(uz - wx)
            endif
            
            if (n==3) then
               uy =  jm*Qhat(i,j,k,1)
               vx =  im*Qhat(i,j,k,2)
               p(i,j,k) = pi2*(vx - uy)
            endif
         enddo
      enddo
   enddo
   call btran_c2r(p,rhsg(1,1,1,n))
enddo
end subroutine ns_vorticity







subroutine p3_params_init
use params
use p3dfft
implicit none
character(len=80) message
integer :: fail=0
integer :: dims(2)
integer :: isize(3),fsize(3),p3_istart(3),p3_iend(3)
integer :: p3_fstart(3),p3_fend(3)
integer :: i,j,k,iw,jw,kw


! check that our reference decomposition matches the p3dfft gridspace
! decomposition
if (ncpu_x /= 1) then
   fail=1
   call print_message("P3DFFT requires nproc_z = 1")
endif
if (ncpu_y > ncpu_z) then
   fail=1
   call print_message("P3DFFT requires ncpu_y < ncpu_z")
endif
dims(1)=ncpu_y
dims(2)=ncpu_z

! set last arg to .true. to allow overwriting input for the
! inverse FFT.  
call p3dfft_setup(dims,g_nx,g_ny,g_nz,.false.)
call get_dims(p3_istart,p3_iend,isize,1)
call get_dims(p3_fstart,p3_fend,fsize,2)


if ( isize(1)*isize(2)*real(isize(3),r8kind)  > nx*ny*real(nz,r8kind) ) then
   fail=1
   call print_message("P3DFFT error: z padding insufficient")
endif
if ( 2*fsize(1)*fsize(2)*real(fsize(3),r8kind)  > nx*ny*real(nz,r8kind) ) then
   fail=1
   call print_message("P3DFFT error: z padding insufficient")
endif

if (nx1 /= 1 ) then
   fail=1
   call print_message("P3DFFT error: x front padding must be 0")
endif
if (ny1 /= 1 ) then
   fail=1
   call print_message("P3DFFT error: y front padding must be 0")
endif
if (nz1 /= 1 ) then
   fail=1
   call print_message("P3DFFT error: z front padding must be 0")
endif
if (nx2 /= isize(1) ) then
   fail=1
   call print_message("P3DFFT error: x ref decomp dimension does not match p3dfft")
endif
if (ny2 /= isize(2) ) then
   fail=1
   call print_message("P3DFFT error: y ref decomp dimension does not match p3dfft")
endif
if (nz2 /= isize(3) ) then
   fail=1
   call print_message("P3DFFT error: z ref decomp dimension does not match p3dfft")
endif


if (nx /= isize(1) ) then
   fail=1
   call print_message("P3DFFT error: x end padding must be 0")
endif
if (ny /= isize(2) ) then
   fail=1
   call print_message("P3DFFT error: y end padding must be 0")
endif
! we can have padding in z direction
if (nz <  isize(3) ) then
   fail=1
   call print_message("P3DFFT error: z dimmension too small for p3dfft")
endif

#if 0
do i=0,ncpu_x*ncpu_y*ncpu_z-1
if (my_pe==i) then
   write(*,'(i4,a,3i3,a,3i3)') my_pe," decomp=",my_x,my_y,my_z," p3=",&
   (p3_istart(1)-1)/isize(1),&
   (p3_istart(2)-1)/isize(2),&
   (p3_istart(3)-1)/isize(3)
endif
call mpi_barrier(comm_3d,j)
enddo
#endif

! check physical grid
if ( nint(xcord(nx1)/delx) /= p3_istart(1)-1 ) then
   write(*,'(a,4i3,2i5)') 'ERROR: x start index: ',my_pe,my_x,my_y,my_z,&
        nint(xcord(nx1)/delx),p3_istart(1)-1
   fail = 1
endif
if ( nint(xcord(nx2)/delx) /= p3_iend(1)-1 ) then
   print *,'ERROR: x grid end ',nint(xcord(nx2)/delx),p3_iend(1)-1
   fail = 1
endif
if ( nint(ycord(ny1)/dely) /= p3_istart(2)-1 ) then
   write(*,'(a,4i3,2i5)') 'ERROR: y start index: ',my_pe,my_x,my_y,my_z,&
        nint(ycord(ny1)/dely),p3_istart(2)-1
   fail = 1
endif
if ( nint(ycord(ny2)/dely) /= p3_iend(2)-1 ) then
   print *,'ERROR: y grid end ',nint(ycord(ny2)/dely),p3_iend(2)-1
   fail = 1
endif
if ( nint(zcord(nz1)/delz) /= p3_istart(3)-1 ) then
   write(*,'(a,4i3,2i5)') 'ERROR: z start index: ',my_pe,my_x,my_y,my_z,&
        nint(zcord(nz1)/delz),p3_istart(3)-1
   fail = 1
endif
if ( nint(zcord(nz2)/delz) /= p3_iend(3)-1 ) then
   print *,'ERROR: z grid end ',nint(zcord(nz2)/delz),p3_iend(3)-1
   fail = 1
endif


! p3dfft Fourier space loop dimensions
p3_n1 = fsize(1)
p3_n2 = fsize(2)
p3_n3 = fsize(3)

! p3dfft arrays:
allocate(p3_1mcord(p3_n1))
allocate(p3_2mcord(p3_n2))
allocate(p3_3mcord(p3_n3))


! initialize the wave numbers:
! note: we assume P3DFFT was compiled with STRIDE1 option.
! that means 1st dimension is really Z direction
!            2nd dimension is really X direction
!            3nd dimension is really Y direction
do i=1,p3_n1
   iw = p3_fstart(1)+i-1
   p3_1mcord(i) = iw-1
   if (iw > g_nz/2 ) p3_1mcord(i) = -(g_nz-iw+1)
   !print *,i,p3_1mcord(i)
enddo
do j=1,p3_n2
   jw = p3_fstart(2)+j-1
   p3_2mcord(j) = jw-1
   if (jw > g_nx/2 ) p3_2mcord(j) = -(g_nx-jw+1)
   !print *,j,p3_2mcord(j)
enddo
do k=1,p3_n3
   kw = p3_fstart(3)+k-1
   p3_3mcord(k) = kw-1
   if (kw > g_ny/2 ) p3_3mcord(k) = -(g_ny-kw+1)
   !print *,k,p3_3mcord(k)
enddo
if (io_pe==my_pe) then
   write(*,'(a,3i6)') 'p3dfft local grid space dims:    ',isize
   write(*,'(a,3i6)') 'p3dfft local fourier space dims: ',fsize
endif

! 32^3  y direction:  0..15,-16   0,1,2,3,4,5,6,7  8,9,10,11,12,13,14,15,-16
!   print *,my_pe,'imcord; ',p3_1mcord
!   print *,my_pe,'jmcord; ',p3_2mcord
!   print *,my_pe,'kmcord; ',p3_3mcord

call mpi_barrier(comm_3d,i)
if (fail/=0) call abortdns("p3dfft params init dimension settings failure")
end subroutine


subroutine p3_fft3d_nvar(Q_grid,Q,work1)
! convinience function for routines which only know Qn
! with dimensions nx,ny,nz
!
! input:  Q_grid:   state vector in grid space (nx,ny,nz) decomposition
! output: Q         state vector in spectral space, z-pencil decomposition
use params
implicit none
real*8 :: Q_grid(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
complex*16 :: Q(p3_n1,p3_n2,p3_n3,n_var)
integer n,im,jm,km,i,j,k
do n=1,n_var
   call ftran_r2c(Q_grid(1,1,1,n),Q(1,1,1,n))
   Q(:,:,:,n)=Q(:,:,:,n)/g_nx/g_ny/g_nz
!      call btran_c2r(Q(1,1,1,n),work1)
!      print *,'max error = ',maxval(&
!      abs( work1(nx1:nx2,ny1:ny2,nz1:nz2)-Q_grid(nx1:nx2,ny1:ny2,nz1:nz2,n) ) )
enddo

#if 0
! debug code: fft of sine(kx), print modes
n=1
do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
   Q_grid(i,j,k,n) = sin( 5*pi2*zcord(k) ) 
enddo
enddo
enddo
call ftran_r2c(Q_grid(1,1,1,n),Q(1,1,1,n))
Q(:,:,:,n)=Q(:,:,:,n)/g_nx/g_ny/g_nz

do k=1,p3_n3
   jm=p3_3mcord(k)   ! y derivative
   do j=1,p3_n2
      im=p3_2mcord(j)   ! x deriviative 
      do i=1,p3_n1  
         km=p3_1mcord(i)  ! z deriviative
         
         if (abs(Q(i,j,k,n)).gt.1e-8) then
            write(*,'(3i5,2f20.10)') im,jm,km,real(Q(i,j,k,n)),imag(Q(i,j,k,n))
         endif
      enddo
   enddo
enddo
stop
#endif
end subroutine


