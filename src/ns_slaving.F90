#include "macros.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine to take one Runge-Kutta 4th order time step using the
! slaving technique of Frisch and Morf; Flow assumes arbitrary
! hyperviscosity power as defined in the input file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk4(time,Q_grid,Q,Q_tmp,Q_old,rhs,work1,work2)
use params
implicit none
real*8 :: time
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Q_grid(nx,ny,nz,n_var)
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
   if (smagorinsky>0) then
      call abortdns("Error: ns3dspectral model does not yet support Smagorinsky")
   endif

   if (dealias==0) then
      call abortdns("Error: using ns3dspectral model, which must be run dealiased")
   endif
   if (numerical_method/=FOURIER) then
      call abortdns("Error: ns3dspectral model requires FFT method.")
   endif
   if (equations/=NS_UVW) then
      call print_message("Error: ns3dspectral model can only runs equations==NS_UVW")
      call abortdns("initial conditions are probably incorrect.")
   endif
   if (alpha_value/=0) then
      call abortdns("Error: alpha>0 but this is not the alpha model!")
   endif
   if (npassive>0) call abortdns("Error: ns_slaving not yet coded for passive scalars")
   if (n_var<3) call abortdns("Error: ns_slaving requires n_var>=3")
   if (fcor>0) call abortdns("Error: ns_slaving not coded for rotation")
   if (use_phaseshift) then
      allocate(q2(nx,ny,nz,ndim))
   endif

   ! intialize Q with Fourier Coefficients:
   call z_fft3d_nvar(Q_grid,Q,work1,work2) 
endif

call rk4reshape(time,Q_grid,Q,rhs,rhs,Q_tmp,Q_old,work1,work2,q2)
end




subroutine rk4reshape(time,Q_grid,Q,rhs,rhsg,Q_tmp,Q_old,work,work2,q2)
use params
implicit none

!notes on slaving technique to be added here:
!stage 1: 
!
!
!
!
!stage 2:
!
!
!
!
!stage 3:
!
!
!
!
!stage 4:
!

real*8 :: time
real*8 :: Q_grid(nx,ny,nz,n_var)
real*8 :: Q(g_nz2,nx_2dz,ny_2dz,n_var)
real*8 :: Q_tmp(g_nz2,nx_2dz,ny_2dz,n_var),Qb_tmp(g_nz2,nx_2dz,ny_2dz,n_var)
real*8 :: Q_old(g_nz2,nx_2dz,ny_2dz,n_var)
real*8 :: work(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: q2(nx,ny,nz,ndim)  ! only allocated if phase shifting turned on
                             ! only to be used for phase shifting
! overlapped in memory:
real*8 :: rhs(g_nz2,nx_2dz,ny_2dz,n_var),rhsb(g_nz2,nx_2dz,ny_2dz,n_var)
real*8 :: rhsg(nx,ny,nz,n_var)


! local variables
real*8 :: ke_old,time_old,vel
integer i,j,k,n,ierr
integer n1,n1d,n2,n2d,n3,n3d,im,jm,km

real*8 xw,xw2,xw_viss
real*8 exp_f,exph_f,cc,cf,ch,ch_h,small,dth

small=1.0d-4
!delt = 1e-2
dth=delt/2
	

! stage 1
call ns3D(rhs,rhsg,Q,Q_grid,time,1,work,work2,1,q2)

do n=1,n_var
   
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nx_2dz
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

            xw=(im*im + jm*jm + km*km/Lz/Lz)*pi2_squared
            xw_viss=mu*xw
            if (mu_hyper>=2 ) then
               xw2=(im*im*pi2_squared)
               xw2=xw2+(jm*jm*pi2_squared)
               xw2=xw2+(km*km*pi2_squared/(Lz*Lz))
               xw_viss=xw_viss + mu_hyper_value*xw2**mu_hyper
            endif
            if (mu_hyper==0) then
               xw_viss=xw_viss + mu_hyper_value
            endif
            if (mu_hypo==1 .and. xw>0) then
               xw_viss=xw_viss + mu_hypo_value/xw
            endif


            ch = -xw_viss*delt
	    ch_h=ch/2
            exph_f = exp(ch_h)
            exp_f = exp(ch)	   
	    if( -ch .le. small ) then
               cc = (  1.0d0 + ch_h/2.0d0 &
                     &  + ch_h**2/6.0d0 + ch_h**3/24.0d0 + ch_h**4/120.0d0  ) * dth
	       cf = (  840.0d0 + 840.0d0*ch + 378.0d0*ch**2 + 112.0d0*ch**3 + 25.0d0*ch**4  ) &
		     &   /5040.0d0*delt				
            else
                cc = ( exph_f - 1.0d0 )/ch_h*dth
		cf = (  -4.0d0- ch +exp_f*( 4.0d0- 3.0d0*ch + ch**2)  )/ch/ch/ch*delt					  
            endif
      
			      
      Q_old(k,i,j,n)=Q(k,i,j,n)
      Q(k,i,j,n)=exp_f *Q_old(k,i,j,n)+rhs(k,i,j,n)*cf             ! accumulate RHS
	  Q_tmp(k,i,j,n) = exph_f*Q_old(k,i,j,n) + cc*rhs(k,i,j,n)  ! Euler timestep with dt=delt/2
	  Qb_tmp(k,i,j,n)=Q_tmp(k,i,j,n)
	  rhsb(k,i,j,n)=rhs(k,i,j,n)
   enddo
   enddo
   enddo

enddo


do n=1,n_var
   call z_ifft3d(Q_tmp(1,1,1,n),Q_grid(1,1,1,n),work)
enddo



! stage 2
call ns3D(rhs,rhsg,Q_tmp,Q_grid,time+delt/2.0,0,work,work2,2,q2)

do n=1,n_var
   
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nx_2dz
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

            xw=(im*im + jm*jm + km*km/Lz/Lz)*pi2_squared
            xw_viss=mu*xw
            if (mu_hyper>=2) then
               xw2=(im*im*pi2_squared)
               xw2=xw2+(jm*jm*pi2_squared)
               xw2=xw2+(km*km*pi2_squared/(Lz*Lz))
               xw_viss=xw_viss + mu_hyper_value*xw2**mu_hyper
            endif
            if (mu_hyper==0) then
               xw_viss=xw_viss + mu_hyper_value
            endif
            if (mu_hypo==1 .and. xw>0) then
               xw_viss=xw_viss + mu_hypo_value/xw
            endif
			ch = -xw_viss*delt
			ch_h=ch/2
           exph_f = exp(ch_h)
           exp_f = exp(ch)	
		   if( -ch .le. small )then
                      cc = (  1.0d0 + ch_h/2.0d0 &
                &  + ch_h**2/6.0d0 + ch_h**3/24.0d0 + ch_h**4/120.0d0 ) * dth
			cf = (  840.0d0 + 420.0d0*ch + 126.0d0*ch**2 + 28.0d0*ch**3 + 5.0d0*ch**4   ) &
                &   /5040.0d0*delt				
                   else
                      cc = ( exph_f - 1.0d0 )/ch_h*dth
				cf = (  2.0d0+ ch +exp_f*(-2.0d0+ ch)  )/ch/ch/ch*delt					  
                   endif
				   		   
      Q(k,i,j,n)=Q(k,i,j,n)+2*rhs(k,i,j,n)*cf             ! accumulate RHS
      !Q_tmp(k,i,j,n)=Q_old(k,i,j,n) + delt*rhs(k,i,j,n)/2   ! Euler timestep with dt=delt/2
	  Q_tmp(k,i,j,n) = exph_f * Q_old(k,i,j,n) + cc * rhs(k,i,j,n)	  
   enddo
   enddo
   enddo
enddo
do n=1,n_var
   call z_ifft3d(Q_tmp(1,1,1,n),Q_grid(1,1,1,n),work)
enddo


! stage 3
call ns3D(rhs,rhsg,Q_tmp,Q_grid,time+delt/2,0,work,work2,3,q2)

do n=1,n_var
   
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nx_2dz
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

            xw=(im*im + jm*jm + km*km/Lz/Lz)*pi2_squared
            xw_viss=mu*xw
            if (mu_hyper>=2 ) then
               xw2=(im*im*pi2_squared)
               xw2=xw2+(jm*jm*pi2_squared)
               xw2=xw2+(km*km*pi2_squared/(Lz*Lz))
               xw_viss=xw_viss + mu_hyper_value*xw2**mu_hyper
            endif
            if (mu_hyper==0) then
               xw_viss=xw_viss + mu_hyper_value
            endif
            if (mu_hypo==1 .and. xw>0) then
               xw_viss=xw_viss + mu_hypo_value/xw
            endif
			ch = -xw_viss*delt
			ch_h=ch/2
           exph_f = exp(ch_h)
           exp_f = exp(ch)	
		   if( -ch .le. small )then
                      cc = (  1.0d0 + ch_h/2.0d0 &
                &  + ch_h**2/6.0d0 + ch_h**3/24.0d0 + ch_h**4/120.0d0 ) * dth
			cf = (  840.0d0 + 420.0d0*ch + 126.0d0*ch**2 + 28.0d0*ch**3 + 5.0d0*ch**4   ) &
                &   /5040.0d0*delt				
                   else
                      cc = ( exph_f - 1.0d0 )/ch_h*dth
				cf = (  2.0d0+ ch +exp_f*(-2.0d0+ ch)  )/ch/ch/ch*delt					  
                   endif
				   		   
      Q(k,i,j,n)=Q(k,i,j,n)+2*rhs(k,i,j,n)*cf             ! accumulate RHS
      !Q_tmp(k,i,j,n)=Q_old(k,i,j,n) + delt*rhs(k,i,j,n)/2   ! Euler timestep with dt=delt/2
      Q_tmp(k,i,j,n) = exph_f * Qb_tmp(k,i,j,n) + cc*( 2.0d0 * rhs(k,i,j,n) - &
                &		rhsb(k,i,j,n) ) 	  
   enddo
   enddo
   enddo
enddo

do n=1,n_var
   call z_ifft3d(Q_tmp(1,1,1,n),Q_grid(1,1,1,n),work)
enddo

! stage 4
call ns3D(rhs,rhsg,Q_tmp,Q_grid,time+delt,0,work,work2,4,q2)


do n=1,n_var
   
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nx_2dz
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

            xw=(im*im + jm*jm + km*km/Lz/Lz)*pi2_squared
            xw_viss=mu*xw
            if (mu_hyper>=2 ) then
               xw2=(im*im*pi2_squared)
               xw2=xw2+(jm*jm*pi2_squared)
               xw2=xw2+(km*km*pi2_squared/(Lz*Lz))
               xw_viss=xw_viss + mu_hyper_value*xw2**mu_hyper
            endif
            if (mu_hyper==0) then
               xw_viss=xw_viss + mu_hyper_value
            endif
            if (mu_hypo==1 .and. xw>0) then
               xw_viss=xw_viss + mu_hypo_value/xw
            endif
			ch = -xw_viss*delt
			ch_h=ch/2
           exph_f = exp(ch_h)
           exp_f = exp(ch)	
		   if( -ch .le. small )then
			cf = (  840.0d0 - 42.0d0*ch**2 - 14.0d0*ch**3 - 3.0d0*ch**4  ) &
                &   /5040.0d0 * delt			
                   else
				cf = (  -4.0d0-3.0d0*ch - ch**2.0d0+exp_f*(4.0d0-ch)  )/ch/ch/ch* delt				  
                   endif

				   		   
      Q(k,i,j,n)=Q(k,i,j,n)+rhs(k,i,j,n)*cf             ! accumulate RHS
   enddo                                                   ! RHS = weighted avarge of 4 RHS computed above
   enddo
   enddo
enddo
do n=1,n_var
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
use sforcing
use spectrum
implicit none

! input
real*8 time
integer compute_ints,rkstage

! input, but data can be trashed if needed
real*8 Qhat(g_nz2,nx_2dz,ny_2dz,n_var)           ! Fourier data at time t
real*8 Q(nx,ny,nz,n_var)                         ! grid data at time t

! output  (rhsg and rhs are overlapped in memory)
real*8 rhs(g_nz2,nx_2dz,ny_2dz,n_var)
real*8 rhsg(nx,ny,nz,n_var)    

! work/storage
real*8  :: q2(nx,ny,nz,ndim)  ! only allocated if use_phaseshift
real*8 work(nx,ny,nz)
! actual dimension: nx,ny,nz, since sometimes used as work array
real*8 p(g_nz2,nx_2dz,ny_2dz)    



                                 

!local

real*8 xw,xw2,xfac,tmx1,tmx2,xw_viss
real*8 uu,vv,ww,dummy
integer n,i,j,k,im,km,jm,ns
integer n1,n1d,n2,n2d,n3,n3d
real*8 :: ke,uxx2ave,ux2ave,ensave,vorave,helave,maxvor,ke_diss,u2,ens_alpha
real*8 :: p_diss(n_var),pke(n_var)
real*8 :: h_diss,ux,uy,uz,vx,vy,vz,wx,wy,wz,ens_diss2,ens_diss4,ens_diss6
real*8 :: f_diss=0,a_diss=0,fxx_diss=0
real*8,save :: f_diss_ave
real*8 :: vor(3),xi


!
! NOTE: for Fourier Coefficients with mode  im=g_nx/2, this is the
! sole "cosine" mode.  Its derivative maps to the g_nx/2 sine mode
! which aliases to 0 on the grid.  So this coefficient should be
! ignored.  This subroutine assumes everything is dealiased, and so
! this mode will be removed no matter what we do with it, and so
! we just dont wory about the case when im=g_nx/2.
!
call wallclock(tmx1)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid space computations.  
! 
! form u x (vor+f), overwrite into Q, ifft into rhs 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call ns_vorticity(rhsg,Qhat,work,p)

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
   
   maxvor = max(maxvor,abs(vor(1)))
   maxvor = max(maxvor,abs(vor(2)))
   maxvor = max(maxvor,abs(vor(3)))

   if (use_phaseshift) vor=vor/2   ! half contribution from this grid
                                   ! half contribution from phase shifted grid

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
! phase shifting:  now compute u x vor with a phase shift,
! and add to RHS:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (use_phaseshift) then
   ! compute phaseshifted (u,v,w), store in Q()
   do n=1,3
      call z_phaseshift(Qhat(1,1,1,n),1,work)  ! phaseshift Qhat
      call z_ifft3d(Qhat(1,1,1,n),Q(1,1,1,n),work)
   enddo
   ! compute phaseshifted vorticity, store in q2:
   call ns_vorticity(q2,Qhat,work,p)

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
      call z_fft3d_trashinput(Q(1,1,1,n),p,work)
      call z_phaseshift(p,-1,work)             ! un-phaseshift p
      rhs(:,:,:,n) = rhs(:,:,:,n) + p/2
      call z_phaseshift(Qhat(1,1,1,n),-1,work) ! restore Qhat to unphaseshifted version
   enddo
endif









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  compute spectra of RHS (u dot (u cross omega) )
!  used to compute transfer function by diagnostics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (compute_ints==1 .and. compute_transfer) then
   spec_diff=0
   do n=1,ndim
      call compute_spectrum_z_fft(Qhat(1,1,1,n),rhs(1,1,1,n),spec_tmp)
      spec_diff=spec_diff+spec_tmp
   enddo
endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute dissipation related diagnostics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ke=0
ux2ave=0
ke_diss=0
h_diss=0
ens_diss2=0
ens_diss4=0
ens_diss6=0
uxx2ave=0
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nx_2dz
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

            xw=(im*im + jm*jm + km*km/Lz/Lz)*pi2_squared
            xw_viss=mu*xw
            if (mu_hyper>=2  ) then
               xw2=(im*im*pi2_squared)
               xw2=xw2+(jm*jm*pi2_squared)
               xw2=xw2+(km*km*pi2_squared/(Lz*Lz))
               xw_viss=xw_viss + mu_hyper_value*xw2**mu_hyper
            endif
            if (mu_hyper==0) then
               xw_viss=xw_viss + mu_hyper_value
            endif
            if (mu_hypo==1 .and. xw>0) then
               xw_viss=xw_viss + mu_hypo_value/xw
            endif

            if (compute_ints==1) then

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


             endif           

      enddo
   enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  spec_tmp:  spectrum of u dot (u cross omega + grad(pi) + diffusion) 
!  since spectrum of u dot grad(pi) is zero, we can now take the
!  differece between spec_tmp and spec_diff to get the diffusion
!  spectrum.  (stored in spec_diff)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (compute_ints==1 .and. compute_transfer) then
   spec_diff=-spec_diff
   spec_f=0
   do n=1,ndim
      call compute_spectrum_z_fft(Qhat(1,1,1,n),rhs(1,1,1,n),spec_tmp)
      spec_diff=spec_diff+spec_tmp
      spec_f=spec_f+spec_tmp
   enddo
   ! spec_f = (for now) spectrum of advection + diffusion terms
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! add in forcing to rhs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute information for transfer spectrum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! spec_f   = spectrum of advection + diffusion terms
! spec_tmp = spectrum of advection + diffusion terms + forcing terms
! so compute the difference and store in spec_f to get the spectrum
! of just the forcing terms.  
if (compute_ints==1 .and. compute_transfer) then
   spec_f=-spec_f
   do n=1,ndim
      call compute_spectrum_z_fft(Qhat(1,1,1,n),rhs(1,1,1,n),spec_tmp)
      spec_f=spec_f+spec_tmp
   enddo
endif


!  make rhs div-free
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nx_2dz
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
   do i=1,nx_2dz
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

      enddo
   enddo
enddo


! compute spectrum of u dot RHS, store in spec_rhs
if (compute_ints==1 .and. compute_transfer) then
   spec_rhs=0
   do n=1,ndim
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











subroutine ns_vorticity(rhsg,Qhat,work,p)
use params
implicit none
real*8 Qhat(g_nz2,nx_2dz,ny_2dz,n_var)           ! Fourier data at time t
real*8 rhsg(nx,ny,nz,n_var)    
real*8 work(nx,ny,nz)
real*8 p(g_nz2,nx_2dz,ny_2dz)    

!local
integer n,i,j,k,im,km,jm
real*8 ux,uy,uz,wx,wy,wz,vx,vy,vz,uu,vv,ww



do n=1,3
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nx_2dz
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







