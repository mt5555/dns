#include "macros.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine to take one Runge-Kutta 4th order time step
!
!  This is a copy of ns.F90, but for efficiency sake it will
!  work (when in grid space) with the x-pencil decompostion
!  instead of our reference (nx,ny,nz) decompostion.
!
!  This is to avoid the transform_from_x() call in the 3D FFTs
!  which is expensive when using a pencil decopostion.
!
!  This routine should be used when the parallel decomposition is
!
!  N1 x 1 x N2
!
!  In this case, the amount of communcation is reduced by 33%
!
!  When using a slab decomposition, this routine will not save any 
!  communication, but it will save some on-processor memory copies
!
!  ns_xpencil.F90 does not support:
!    passive scalars
!    hyper viscosity
!    rotation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk4(time,Q_grid,Q,Q_tmp,Q_old,rhsg,work,work2)
use params
use transpose
implicit none
real*8 :: time
real*8 :: Q_grid(nx,ny,nz,n_var)
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: work(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: Q_tmp(nx,ny,nz,n_var)
real*8 :: Q_old(nx,ny,nz,n_var)
real*8 :: rhsg(nx,ny,nz,n_var)

logical,save :: firstcall=.true.
integer n


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
   if (alpha_value/=0) then
      call abort("Error: alpha>0 but this is not the alpha model!")
   endif
   if (npassive>0) then
      call abort("Error: dnsp (x-pensil) model with tracers not yet coded")
   endif

   if (data_x_pencils) call abort("ns_xpencil: this is not possible")

   ! intialize Q with Fourier Coefficients:
   call z_fft3d_nvar(Q_grid,Q,work,work2) 

   ! enable more efficient vorticity routine
   if (ncpu_y == 1 ) then
	use_vorticity3=.true. 
	call print_message("Enabling use_vorticity3 (more efficient) option.")
   endif	
endif

if (.not. data_x_pencils) then
   ! initial data, or after output/diagnostics, data in Q_grid may be
   ! stored in nx,ny,nz decompostion instead of x-pencils.  
   ! convert data to x-pencil decompostion:
   Q_tmp=Q_grid;   call transpose_to_x_3d(Q_tmp,Q_grid)
   data_x_pencils=.true.
endif



if (use_vorticity3) then
   ! call fast version that computes Vx and Uy during U iffts:
   ! saves 1 z transforms per RK stage
   call rk4reshape2(time,Q_grid,Q,rhsg,rhsg,Q_tmp,Q_old,work,work2)
else
   call rk4reshape(time,Q_grid,Q,rhsg,rhsg,Q_tmp,Q_old,work,work2)
endif
end






subroutine rk4reshape(time,Q_grid,Q,rhs,rhsg,Q_tmp,Q_old,work,work2)
use params
implicit none
real*8 :: time
!real*8 :: Q_grid(nx,ny,nz,n_var)
real*8 :: Q_grid(g_nx2,nslabz,ny_2dx,n_var)
real*8 :: Q(g_nz2,nx_2dz,ny_2dz,n_var)
real*8 :: Q_tmp(g_nz2,nx_2dz,ny_2dz,n_var)
real*8 :: Q_old(g_nz2,nx_2dz,ny_2dz,n_var)
real*8 :: work(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

! overlapped in memory: (real size: nx,ny,nz)
real*8 :: rhs(g_nz2,nx_2dz,ny_2dz,n_var)
real*8 :: rhsg(g_nx2,nslabz,ny_2dx,n_var)



! local variables
real*8 :: ke_old,time_old,vel,dummy
integer i,j,k,n,ierr
integer n1,n1d,n2,n2d,n3,n3d,im,jm,km


! stage 1
call ns3D(rhs,rhs,Q,Q_grid,time,1,work,work2,1,dummy,dummy)

do n=1,n_var
   do j=1,ny_2dz
   do i=1,nx_2dz
   do k=1,g_nz
      Q_old(k,i,j,n)=Q(k,i,j,n)
      Q(k,i,j,n)=Q(k,i,j,n)+delt*rhs(k,i,j,n)/6.0
      Q_tmp(k,i,j,n)=Q_old(k,i,j,n) + delt*rhs(k,i,j,n)/2.0
   enddo
   enddo
   enddo
   call zx_ifft3d(Q_tmp(1,1,1,n),Q_grid(1,1,1,n),work)
enddo






! stage 2
call ns3D(rhs,rhs,Q_tmp,Q_grid,time+delt/2.0,0,work,work2,2,dummy,dummy)

do n=1,n_var
   do j=1,ny_2dz
   do i=1,nx_2dz
   do k=1,g_nz
      Q(k,i,j,n)=Q(k,i,j,n)+delt*rhs(k,i,j,n)/3
      Q_tmp(k,i,j,n)=Q_old(k,i,j,n) + delt*rhs(k,i,j,n)/2
   enddo
   enddo
   enddo


   call zx_ifft3d(Q_tmp(1,1,1,n),Q_grid(1,1,1,n),work)
enddo

! stage 3
call ns3D(rhs,rhs,Q_tmp,Q_grid,time+delt/2,0,work,work2,3,dummy,dummy)

do n=1,n_var
   do j=1,ny_2dz
   do i=1,nx_2dz
   do k=1,g_nz
      Q(k,i,j,n)=Q(k,i,j,n)+delt*rhs(k,i,j,n)/3
      Q_tmp(k,i,j,n)=Q_old(k,i,j,n)+delt*rhs(k,i,j,n)
   enddo
   enddo
   enddo

   call zx_ifft3d(Q_tmp(1,1,1,n),Q_grid(1,1,1,n),work)
enddo

! stage 4
call ns3D(rhs,rhs,Q_tmp,Q_grid,time+delt,0,work,work2,4,dummy,dummy)


do n=1,n_var
   do j=1,ny_2dz
   do i=1,nx_2dz
   do k=1,g_nz
      Q(k,i,j,n)=Q(k,i,j,n)+delt*rhs(k,i,j,n)/6
   enddo
   enddo
   enddo
   call zx_ifft3d(Q(1,1,1,n),Q_grid(1,1,1,n),work)
enddo

time = time + delt



! compute max U  
maxs(1:4)=0
maxs(10:11)=-9e20
do k=1,ny_2dx
do j=1,nslabz
do i=1,g_nx
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






subroutine rk4reshape2(time,Q_grid,Q,rhs,rhsg,Q_tmp,Q_old,work,work2)
use params
implicit none
real*8 :: time
!real*8 :: Q_grid(nx,ny,nz,n_var)
real*8 :: Q_grid(g_nx2,nslabz,ny_2dx,n_var)
real*8 :: Q(g_nz2,nx_2dz,ny_2dz,n_var)
real*8 :: Q_tmp(g_nz2,nx_2dz,ny_2dz,n_var)
real*8 :: Q_old(g_nz2,nx_2dz,ny_2dz,n_var)
real*8 :: work(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

! overlapped in memory: (real size: nx,ny,nz)
real*8 :: rhs(g_nz2,nx_2dz,ny_2dz,n_var)
real*8 :: rhsg(g_nx2,nslabz,ny_2dx,n_var)


! local variables
real*8 :: ke_old,time_old,vel
integer i,j,k,n,ierr
integer n1,n1d,n2,n2d,n3,n3d,im,jm,km

! storage for first two components of vorticity, needed for
! use_vorticity3 trick:
real*8,save,allocatable  :: uygrid(:,:,:)
real*8,save,allocatable  :: vxgrid(:,:,:)
logical,save :: firstcall=.true.

if (firstcall) then
   firstcall=.false.
   ! initialize uygrid,vxgrid
   call print_message('allocating extra arrays for use_vorticity3')
   allocate(uygrid(nx,ny,nz))
   allocate(vxgrid(nx,ny,nz))
   call zx_ifft3d_and_dy(Q(1,1,1,1),work2,uygrid,work)
   call zx_ifft3d_and_dx(Q(1,1,1,2),work2,vxgrid,work)
endif

! stage 1
call ns3D(rhs,rhsg,Q,Q_grid,time,1,work,work2,1,uygrid,vxgrid)
do n=1,n_var
   do j=1,ny_2dz
   do i=1,nx_2dz
   do k=1,g_nz
      Q_old(k,i,j,n)=Q(k,i,j,n)
      Q(k,i,j,n)=Q(k,i,j,n)+delt*rhs(k,i,j,n)/6.0
      Q_tmp(k,i,j,n)=Q_old(k,i,j,n) + delt*rhs(k,i,j,n)/2.0
   enddo
   enddo
   enddo
enddo

call zx_ifft3d_and_dy(Q_tmp(1,1,1,1),Q_grid(1,1,1,1),uygrid,work)
call zx_ifft3d_and_dx(Q_tmp(1,1,1,2),Q_grid(1,1,1,2),vxgrid,work)
do n=3,n_var
   call zx_ifft3d(Q_tmp(1,1,1,n),Q_grid(1,1,1,n),work)
enddo






! stage 2
call ns3D(rhs,rhsg,Q_tmp,Q_grid,time+delt/2.0,0,work,work2,2,uygrid,vxgrid)

do n=1,n_var
   do j=1,ny_2dz
   do i=1,nx_2dz
   do k=1,g_nz
      Q(k,i,j,n)=Q(k,i,j,n)+delt*rhs(k,i,j,n)/3
      Q_tmp(k,i,j,n)=Q_old(k,i,j,n) + delt*rhs(k,i,j,n)/2
   enddo
   enddo
   enddo
enddo
call zx_ifft3d_and_dy(Q_tmp(1,1,1,1),Q_grid(1,1,1,1),uygrid,work)
call zx_ifft3d_and_dx(Q_tmp(1,1,1,2),Q_grid(1,1,1,2),vxgrid,work)
do n=3,n_var
   call zx_ifft3d(Q_tmp(1,1,1,n),Q_grid(1,1,1,n),work)
enddo



! stage 3
call ns3D(rhs,rhsg,Q_tmp,Q_grid,time+delt/2,0,work,work2,3,uygrid,vxgrid)

do n=1,n_var
   do j=1,ny_2dz
   do i=1,nx_2dz
   do k=1,g_nz
      Q(k,i,j,n)=Q(k,i,j,n)+delt*rhs(k,i,j,n)/3
      Q_tmp(k,i,j,n)=Q_old(k,i,j,n)+delt*rhs(k,i,j,n)
   enddo
   enddo
   enddo
enddo
call zx_ifft3d_and_dy(Q_tmp(1,1,1,1),Q_grid(1,1,1,1),uygrid,work)
call zx_ifft3d_and_dx(Q_tmp(1,1,1,2),Q_grid(1,1,1,2),vxgrid,work)
do n=3,n_var
   call zx_ifft3d(Q_tmp(1,1,1,n),Q_grid(1,1,1,n),work)
enddo


! stage 4
call ns3D(rhs,rhsg,Q_tmp,Q_grid,time+delt,0,work,work2,4,uygrid,vxgrid)


do n=1,n_var
   do j=1,ny_2dz
   do i=1,nx_2dz
   do k=1,g_nz
      Q(k,i,j,n)=Q(k,i,j,n)+delt*rhs(k,i,j,n)/6
   enddo
   enddo
   enddo
enddo
call zx_ifft3d_and_dy(Q(1,1,1,1),Q_grid(1,1,1,1),uygrid,work)
call zx_ifft3d_and_dx(Q(1,1,1,2),Q_grid(1,1,1,2),vxgrid,work)
do n=3,n_var
   call zx_ifft3d(Q(1,1,1,n),Q_grid(1,1,1,n),work)
enddo



time = time + delt

! compute max U  
maxs(1:4)=0
maxs(10:11)=-9e20
do k=1,ny_2dx
do j=1,nslabz
do i=1,g_nx
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










subroutine ns3d(rhs,rhsg,Qhat,Q_grid,time,compute_ints,work,p,rkstage,uygrid,vxgrid)
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
!real*8 Q(nx,ny,nz,n_var)                         ! grid data at time t
real*8 Q_grid(g_nx2,nslabz,ny_2dx,n_var)
real*8 uygrid(nx,ny,nz)
real*8 vxgrid(nx,ny,nz)


! output  (rhsg and rhs are overlapped in memory)
! true size must be nx,ny,nz,n_var
real*8 rhs(g_nz2,nx_2dz,ny_2dz,n_var)
real*8 rhsg(g_nx2,nslabz,ny_2dx,n_var)

! work/storage
real*8 work(nx,ny,nz)
! actual dimension: nx,ny,nz, since sometimes used as work array
real*8 p(g_nz2,nx_2dz,ny_2dz)    

!local

real*8 xw,xw2,xfac,tmx1,tmx2,xw_viss
real*8 uu,vv,ww,dummy
integer n,i,j,k,im,km,jm,ns
integer n1,n1d,n2,n2d,n3,n3d
real*8 :: ke,uxx2ave,ux2ave,ensave,vorave,helave,maxvor,ke_diss,u2
real*8 :: p_diss(n_var),pke(n_var),enstrophy
real*8 :: h_diss,ux,uy,uz,vx,vy,vz,wx,wy,wz,ens_diss2,ens_diss4,ens_diss6
real*8 :: f_diss=0,a_diss=0,fxx_diss=0
real*8,save :: f_diss_ave
real*8 :: vor(3)

call wallclock(tmx1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! see notes in ns.F90 about pros/cons on how to compute vorticity
!
if (use_vorticity3) then
   k=ny_2dx
   j=nslabz
   i=g_nx
   call ns_vorticity3(rhsg,Qhat,work,p,uygrid,vxgrid)  
else
   call ns_vorticity(rhsg,Qhat,work,p)
endif



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
do k=1,ny_2dx
do j=1,nslabz
do i=1,g_nx
   vor(1)=rhsg(i,j,k,1)
   vor(2)=rhsg(i,j,k,2)
   vor(3)=rhsg(i,j,k,3)

   if (compute_ints==1) then
      vorave = vorave + vor(3)
      helave = helave + Q_grid(i,j,k,1)*vor(1) + & 
           Q_grid(i,j,k,2)*vor(2) + & 
           Q_grid(i,j,k,3)*vor(3)  
      maxvor = max(maxvor,abs(vor(1)))
      maxvor = max(maxvor,abs(vor(2)))
      maxvor = max(maxvor,abs(vor(3)))
   endif


   !  velocity=(u,v,w)  vorticity=(a,b,c)=(wy-vz,uz-wx,vx-uy)
   !  v*(vx-uy) - w*(uz-wx) = (v vx - v uy + w wx) - w uz
   !  w*(wy-vz) - u*(vx-uy)
   !  u*(uz-wx) - v*(wy-vz)
   uu = ( Q_grid(i,j,k,2)*vor(3) - Q_grid(i,j,k,3)*vor(2) )
   vv = ( Q_grid(i,j,k,3)*vor(1) - Q_grid(i,j,k,1)*vor(3) )
   ww = ( Q_grid(i,j,k,1)*vor(2) - Q_grid(i,j,k,2)*vor(1) )

   
   ! overwrite Q with the result
   Q_grid(i,j,k,1) = uu
   Q_grid(i,j,k,2) = vv
   Q_grid(i,j,k,3) = ww

enddo
enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! back to spectral space
! we are done with rhsg - from now on only use rhs (spectral space version)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do n=1,3
   call zx_fft3d_trashinput(Q_grid(1,1,1,n),rhs(1,1,1,n),work)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! From this point on, Q_grid() may be used as a work array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  spec_diff:  spectrum of u dot (u cross omega) 
!              (used as temporary storage for now)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (rkstage==1 .and. compute_transfer) then
   spec_diff=0
   do n=1,ndim
      call compute_spectrum_z_fft(Qhat(1,1,1,n),rhs(1,1,1,n),spec_tmp)
      spec_diff=spec_diff+spec_tmp
   enddo
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! add in diffusion term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ke=0
ux2ave=0
ke_diss=0
h_diss=0
enstrophy=0
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
               
                 enstrophy = enstrophy+ xfac* &           
                       ((wy-vz)**2 + (uz-wx)**2 + (vx-uy)**2)   


               ! incorrect if using hyperviscosity
                  ens_diss2=ens_diss2 + 2*xfac*(mu*xw)*  &
                       ((wy-vz)**2 + (uz-wx)**2 + (vx-uy)**2) 
                  ens_diss4=ens_diss4 + 2*xfac*(mu*xw*xw)*  &
                       ((wy-vz)**2 + (uz-wx)**2 + (vx-uy)**2)
                  ens_diss6=ens_diss6 + 2*xfac*(mu*xw*(xw)**2)*  &
                       ((wy-vz)**2 + (uz-wx)**2 + (vx-uy)**2)      
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
if (rkstage==1 .and. compute_transfer) then
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
         ! trashes Q_grid - used as work array
         call stochastic_highwaveno(Q_grid,Qhat,f_diss,fxx_diss,1)  
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
if (rkstage==1 .and. compute_transfer) then
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
if (rkstage==1 .and. compute_transfer) then
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
   ints(7) = enstrophy
   ints(8)=a_diss    
   ints(9)=fxx_diss                     ! < u_xx,f>
   ints(10)=-ke_diss                 ! <u,u_xx>
   ints(11)=h_diss
   ints(12)=ens_diss2
   ints(13)=ens_diss4
   ints(14)=ens_diss6
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
real*8 Qhat(g_nz2,nx_2dz,ny_2dz,n_var)           ! Fourier data at time t
!real*8 rhsg(nx,ny,nz,n_var)    
real*8 rhsg(g_nx2,nslabz,ny_2dx,n_var)
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
call zx_ifft3d(p,rhsg(1,1,1,n),work)
enddo
end subroutine






subroutine ns_vorticity3(rhsg,Qhat,work,p,uygrid,vxgrid)
!
!  fast version, where vx and uy have already been computed
!
use params
implicit none
real*8 Qhat(g_nz2,nx_2dz,ny_2dz,n_var)           ! Fourier data at time t
real*8 rhsg(g_nx2,nslabz,ny_2dx,n_var)
real*8 work(nx,ny,nz)
real*8 uygrid(g_nz2,nslabz,ny_2dz)    
real*8 vxgrid(g_nz2,nslabz,ny_2dz)    
real*8 p(g_nz2,nx_2dz,ny_2dz)    

!local
integer n,i,j,k,im,km,jm
real*8 ux,uy,uz,wx,wy,wz,vx,vy,vz,uu,vv,ww


! third component of vorticity   vx-uy
do k=1,ny_2dx
do j=1,nslabz
do i=1,g_nx
   rhsg(i,j,k,3)=vxgrid(i,j,k)-uygrid(i,j,k)
enddo
enddo
enddo


! compute 1st and 2nd components in the usual method:
n=1
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nx_2dz
!      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)
         wy = - jm*Qhat(k,i,j+z_jmsign(j),3)
         vz =  - km*Qhat(k+z_kmsign(k),i,j,2)/Lz
         !rhs(k,i,j,1) = pi2*(wy - vz)
         p(k,i,j) = pi2*(wy - vz)
      enddo
   enddo
enddo
call zx_ifft3d(p,rhsg(1,1,1,n),work)

n=2
do j=1,ny_2dz
!   jm=z_jmcord(j)
   do i=1,nx_2dz
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)
         wx = - im*Qhat(k,i+z_imsign(i),j,3)
         uz =  - km*Qhat(k+z_kmsign(k),i,j,1)/Lz
         !rhs(k,i,j,2) = pi2*(uz - wx)
         p(k,i,j) = pi2*(uz - wx)
      enddo
   enddo
enddo
call zx_ifft3d(p,rhsg(1,1,1,n),work)

end subroutine







