#include "macros.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine to take one Runge-Kutta 4th order time step
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk4(time,Q_grid)
use params
use transpose
implicit none
real*8 :: time
real*8 :: Q_grid(nx,ny,nz,n_var)

! local variables
real*8 :: Q_old(nx,ny,nz,n_var)
real*8 :: Q_tmp(nx,ny,nz,n_var)
real*8 :: rhs(nx,ny,nz,n_var)

real*8 :: z_rhs(g_nz2,nslabx,ny_2dz,n_var)
real*8 :: z_Q(g_nz2,nslabx,ny_2dz,n_var)


real*8 :: ke_old,time_old,vel
real*8,save :: Q(nx,ny,nz,n_var)

real*8,save,allocatable :: Q2(:,:,:,:)



integer i,j,k,n,ierr
integer n1,n1d,n2,n2d,n3,n3d
logical,save :: firstcall=.true.



if (firstcall) then
   firstcall=.false.
   allocate(Q2(g_nz2,nslabx,ny_2dz,n_var))

   Q_tmp=Q_grid
   do n=1,3
      call z_fft3d(Q_tmp(1,1,1,n),Q2(1,1,1,n),rhs)  ! use rhs as a work array
      call transpose_from_z(Q2(1,1,1,n),Q(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
   enddo

   if (.not. dealias) then
      call abort("Error: using ns3dspectral model, which must be run dealiased")
   endif
endif







! stage 1

call ns3D(z_rhs,Q2,Q_grid,time,1)
do n=1,3
   call transpose_from_z(z_rhs(1,1,1,n),rhs(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo


do n=1,3
   do j=1,ny_2dz
   do i=1,nslabx
   do k=1,g_nz
      Q2(k,i,j,n)=Q2(k,i,j,n)+delt*z_rhs(k,i,j,n)/6.0
   enddo
   enddo
   enddo



   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      Q_old(i,j,k,n)=Q(i,j,k,n)
      Q(i,j,k,n)=Q(i,j,k,n)+delt*rhs(i,j,k,n)/6.0
      Q_tmp(i,j,k,n) = Q_old(i,j,k,n) + delt*rhs(i,j,k,n)/2.0
      Q_grid(i,j,k,n)=Q_tmp(i,j,k,n)
   enddo
   enddo
   enddo
   call ifft3d(Q_grid(1,1,1,n),rhs)   ! use rhs as a work array
enddo


! stage 2
do n=1,3
   call transpose_to_z(Q_tmp(1,1,1,n),z_Q(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo
call ns3D(z_rhs,z_Q,Q_grid,time+delt/2.0,0)
do n=1,3
   call transpose_from_z(z_rhs(1,1,1,n),rhs(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo



do n=1,3
   do j=1,ny_2dz
   do i=1,nslabx
   do k=1,g_nz
      Q2(k,i,j,n)=Q2(k,i,j,n)+delt*z_rhs(k,i,j,n)/3.0
   enddo
   enddo
   enddo

   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      Q(i,j,k,n)=Q(i,j,k,n)+delt*rhs(i,j,k,n)/3.0
      Q_tmp(i,j,k,n) = Q_old(i,j,k,n) + delt*rhs(i,j,k,n)/2.0
      Q_grid(i,j,k,n)=Q_tmp(i,j,k,n)
   enddo
   enddo
   enddo
   call ifft3d(Q_grid(1,1,1,n),rhs)   ! use rhs as a work array
enddo

! stage 3
do n=1,3
   call transpose_to_z(Q_tmp(1,1,1,n),z_Q(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo
call ns3D(z_rhs,z_Q,Q_grid,time+delt/2.0,0)
do n=1,3
   call transpose_from_z(z_rhs(1,1,1,n),rhs(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo


do n=1,3
   do j=1,ny_2dz
   do i=1,nslabx
   do k=1,g_nz
      Q2(k,i,j,n)=Q2(k,i,j,n)+delt*z_rhs(k,i,j,n)/3.0
   enddo
   enddo
   enddo

   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      Q(i,j,k,n)=Q(i,j,k,n)+delt*rhs(i,j,k,n)/3.0
      Q_tmp(i,j,k,n) = Q_old(i,j,k,n) + delt*rhs(i,j,k,n)
      Q_grid(i,j,k,n)=Q_tmp(i,j,k,n)
   enddo
   enddo
   enddo
   call ifft3d(Q_grid(1,1,1,n),rhs)   ! use rhs as a work array
enddo


! stage 4
do n=1,3
   call transpose_to_z(Q_tmp(1,1,1,n),z_Q(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo
call ns3D(z_rhs,z_Q,Q_grid,time+delt,0)
do n=1,3
   call transpose_from_z(z_rhs(1,1,1,n),rhs(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo


do n=1,3
   do j=1,ny_2dz
   do i=1,nslabx
   do k=1,g_nz
      Q2(k,i,j,n)=Q2(k,i,j,n)+delt*z_rhs(k,i,j,n)/6.0
   enddo
   enddo
   enddo

   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      Q(i,j,k,n)=Q(i,j,k,n)+delt*rhs(i,j,k,n)/6.0
      Q_grid(i,j,k,n)=Q(i,j,k,n)
   enddo
   enddo
   enddo
   call ifft3d(Q_grid(1,1,1,n),rhs)   ! use rhs as a work array
enddo




time = time + delt


! compute KE, max U  
ints_timeU=time
ints(1)=0
maxs(1:4)=0
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   do n=1,3
      ints(1)=ints(1)+.5*Q_grid(i,j,k,n)**2  ! KE
      maxs(n)=max(maxs(n),abs(Q_grid(i,j,k,n)))   ! max u,v,w
   enddo
   vel = abs(Q_grid(i,j,k,1))/delx + abs(Q_grid(i,j,k,2))/dely + abs(Q_grid(i,j,k,3))/delz
   maxs(4)=max(maxs(4),vel)
enddo
enddo
enddo




end subroutine rk4  







subroutine ns3d(rhs,z_Qhat,Q,time,compute_ints)
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
use transpose
implicit none

! input
real*8 z_Qhat(g_nz2,nslabx,ny_2dz,n_var)
real*8 Q(nx,ny,nz,n_var)
real*8 time
integer compute_ints

! output
real*8 rhs(g_nz2,nslabx,ny_2dz,n_var)


!local
real*8 grid(nx,ny,nz,n_var)


real*8 xfac,tmx1,tmx2
real*8 ux,uy,uz,wx,wy,wz,vx,vy,vz,uu,vv,ww
integer n,i,j,k,im,km,jm
integer n1,n1d,n2,n2d,n3,n3d
real*8 :: ke_diss,vor,hel,maxvor



real*8 p(g_nz2,nslabx,ny_2dz,n_var)






call wallclock(tmx1)


!
! NOTE: for Fourier Coefficients with mode  im=g_nx/2, this is the
! sole "cosine" mode.  Its derivative maps to the g_nx/2 sine mode
! which aliases to 0 on the grid.  So this coefficient should be
! ignored.  This subroutine assumes everything is dealiased, and so
! this mode will be removed no matter what we do with it, and so
! we just dont wory about the case when im=g_nx/2.
!


! compute vorticity-hat, store in RHS



do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nslabx
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)
         
         ! u_x term
         vx = - im*z_Qhat(k,i+z_imsign(i),j,2)
         wx = - im*z_Qhat(k,i+z_imsign(i),j,3)
         
         uy = - jm*z_Qhat(k,i,j+z_jmsign(j),1)
         wy = - jm*z_Qhat(k,i,j+z_jmsign(j),3)
         
         uz =  - km*z_Qhat(k+z_kmsign(k),i,j,1)
         vz =  - km*z_Qhat(k+z_kmsign(k),i,j,2)
         
         rhs(k,i,j,1) = pi2*(wy - vz)
         rhs(k,i,j,2) = pi2*(uz - wx)
         rhs(k,i,j,3) = pi2*(vx - uy)

         
      enddo
   enddo
enddo
! 3 forward&back z-transforms, 9 ffts.  
do n=1,3
   call z_ifft3d(rhs(1,1,1,n),grid(1,1,1,n),p)
enddo


! alternative:
! compute all terms above: vx,wx  uy,wy, uz,vz
! 12 ffts, but only 2 f&b z transforms.  








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Grid space computations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

vor=0
hel=0
maxvor=0

   ! rhs = u x vor
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      vor = vor + grid(i,j,k,3)

      hel = hel + Q(i,j,k,1)*grid(i,j,k,1) + & 
                  Q(i,j,k,2)*grid(i,j,k,2) + & 
                  Q(i,j,k,3)*grid(i,j,k,3)  

      maxvor = max(maxvor,abs(grid(i,j,k,1)))
      maxvor = max(maxvor,abs(grid(i,j,k,2)))
      maxvor = max(maxvor,abs(grid(i,j,k,3)))

      uu = ( Q(i,j,k,2)*grid(i,j,k,3) - Q(i,j,k,3)*grid(i,j,k,2) )
      vv = ( Q(i,j,k,3)*grid(i,j,k,1) - Q(i,j,k,1)*grid(i,j,k,3) )
      ww = ( Q(i,j,k,1)*grid(i,j,k,2) - Q(i,j,k,2)*grid(i,j,k,1) )

      grid(i,j,k,1) = uu
      grid(i,j,k,2) = vv
      grid(i,j,k,3) = ww
   enddo
   enddo
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! back to spectral space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do n=1,3
   call z_fft3d(grid(1,1,1,n),rhs(1,1,1,n),p)
enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! add in diffusion term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ke_diss=0
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nslabx
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

            xfac=-mu*(im*im + jm*jm + km*km)*pi2_squared
            rhs(k,i,j,1)=rhs(k,i,j,1) + xfac*z_Qhat(k,i,j,1)
            rhs(k,i,j,2)=rhs(k,i,j,2) + xfac*z_Qhat(k,i,j,2)
            rhs(k,i,j,3)=rhs(k,i,j,3) + xfac*z_Qhat(k,i,j,3)

! < u (uxx + uyy + uzz) > = < u-hat * (uxx-hat + uyy-hat + uzz-hat) >
!                         = < u-hat*u-hat*( im**2 + jm**2 + km**2)

            xfac = 2*2*2*xfac*(g_nx*g_ny*g_nz)
            if (km==0) xfac=xfac/2
            if (jm==0) xfac=xfac/2
            if (im==0) xfac=xfac/2

            ke_diss = ke_diss + xfac*z_Qhat(k,i,j,1)**2 + &
                                xfac*z_Qhat(k,i,j,2)**2 + &
                                xfac*z_Qhat(k,i,j,3)**2 

         

      enddo
   enddo
enddo


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
              - km*rhs(k+z_kmsign(k),i,j,3)
         
         
         ! compute laplacian inverse
         xfac= (im*im +km*km + jm*jm)
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
         uu= - im*p(k,i+imsign(i),j) 
         vv= - jm*p(k,i,j+jmsign(j))
         ww= - km*p(k+kmsign(k),i,j)
         
         rhs(k,i,j,1)=rhs(k,i,j,1) - uu 
         rhs(k,i,j,2)=rhs(k,i,j,2) - vv 
         rhs(k,i,j,3)=rhs(k,i,j,3) - ww 
         
         
         ! dealias           
         if ( ((abs(km)> g_nz/3) .and. (km/=0)) .or. &
              ((abs(jm)> g_ny/3) .and. (jm/=0)) .or. &
              ((abs(im)> g_nx/3) .and. (im/=0)) )  then
            rhs(k,i,j,1)=0
            rhs(k,i,j,2)=0
            rhs(k,i,j,3)=0
         endif

      enddo
   enddo
enddo



if (compute_ints==1) then
   ints_timeDU=time
   ints(3)=ke_diss
   ints(4)=vor
   ints(5)=hel
   maxs(5)=maxvor
endif

call wallclock(tmx2)
tims(5)=tims(5)+(tmx2-tmx1)
end











