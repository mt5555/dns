#include "macros.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine to take one Runge-Kutta 4th order time step
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk4(time,Q_grid,Q,rhs)
use params
implicit none
real*8 :: time
real*8 :: Q_grid(nx,ny,nz,n_var)
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: rhs(nx,ny,nz,n_var)

call rk4reshape(time,Q_grid,Q,rhs,rhs)
end



subroutine rk4reshape(time,Q_grid,Q,rhs,rhsg)
use params
use structf
implicit none
real*8 :: time
real*8 :: Q_grid(nx,ny,nz,n_var)
real*8 :: Q(g_nz2,nslabx,ny_2dz,n_var)
! overlapped in memory:
real*8 :: rhs(g_nz2,nslabx,ny_2dz,n_var)
real*8 :: rhsg(nx,ny,nz,n_var)


! local variables
real*8 :: Q_tmp(g_nz2,nslabx,ny_2dz,n_var)
real*8 :: Q_old(g_nz2,nslabx,ny_2dz,n_var)
real*8 :: work(nx,ny,nz)

real*8 :: ke_old,time_old,vel
integer i,j,k,n,ierr
integer n1,n1d,n2,n2d,n3,n3d
logical,save :: firstcall=.true.
integer,save :: countx=-1,county=-1,countz=-1


if (struct_nx>0) countx=mod(countx+1,struct_nx)  
if (struct_ny>0) county=mod(county+1,struct_ny)  
if (struct_nz>0) countz=mod(countz+1,struct_nz)  

if (firstcall) then
   firstcall=.false.
   do n=1,3
      rhsg(:,:,:,1)=Q_grid(:,:,:,n)
      call z_fft3d_trashinput(rhsg,Q(1,1,1,n),rhsg(1,1,1,2)) ! use rhs as work array
   enddo

   if (.not. dealias) then
      call abort("Error: using ns3dspectral model, which must be run dealiased")
   endif
endif



! stage 1
call ns3D(rhs,rhs,Q,Q_grid,time,1,work)

do n=1,3
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
call ns3D(rhs,rhs,Q_tmp,Q_grid,time+delt/2.0,0,work)

do n=1,3
   do j=1,ny_2dz
   do i=1,nslabx
   do k=1,g_nz
      Q(k,i,j,n)=Q(k,i,j,n)+delt*rhs(k,i,j,n)/3.0
      Q_tmp(k,i,j,n)=Q_old(k,i,j,n) + delt*rhs(k,i,j,n)/2.0
   enddo
   enddo
   enddo


   call z_ifft3d(Q_tmp(1,1,1,n),Q_grid(1,1,1,n),work)
enddo

! stage 3
call ns3D(rhs,rhs,Q_tmp,Q_grid,time+delt/2.0,0,work)


do n=1,3
   do j=1,ny_2dz
   do i=1,nslabx
   do k=1,g_nz
      Q(k,i,j,n)=Q(k,i,j,n)+delt*rhs(k,i,j,n)/3.0
      Q_tmp(k,i,j,n)=Q_old(k,i,j,n)+delt*rhs(k,i,j,n)
   enddo
   enddo
   enddo

   call z_ifft3d(Q_tmp(1,1,1,n),Q_grid(1,1,1,n),work)
enddo

! stage 4
call ns3D(rhs,rhs,Q_tmp,Q_grid,time+delt,0,work)


do n=1,3
   do j=1,ny_2dz
   do i=1,nslabx
   do k=1,g_nz
      Q(k,i,j,n)=Q(k,i,j,n)+delt*rhs(k,i,j,n)/6.0
   enddo
   enddo
   enddo
   !call z_ifft3d(Q(1,1,1,n),Q_grid(1,1,1,n),work)
enddo
call z_ifft3d_str(Q,Q_grid,rhs,rhs,(countx==0),(county==1),(countz==0))


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
ints(1)=ints(1)/g_nx/g_ny/g_nz

end subroutine 







subroutine ns3d(rhs,rhsg,Qhat,Q,time,compute_ints,work)
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
real*8 time
integer compute_ints

! input, but data can be trashed if needed
real*8 Qhat(g_nz2,nslabx,ny_2dz,n_var)           ! Fourier data at time t
real*8 Q(nx,ny,nz,n_var)                         ! grid data at time t

! output  (rhsg and rhs are overlapped in memory)
real*8 rhs(g_nz2,nslabx,ny_2dz,n_var)
real*8 rhsg(nx,ny,nz,n_var)    

! work/storage
real*8 work(nx,ny,nz)
                                 

!local
real*8 p(g_nz2,nslabx,ny_2dz)    
real*8 xfac,tmx1,tmx2
real*8 ux,uy,uz,wx,wy,wz,vx,vy,vz,uu,vv,ww
integer n,i,j,k,im,km,jm
integer n1,n1d,n2,n2d,n3,n3d
real*8 :: ke_diss,vor,hel,maxvor,f_diss







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
! 3 z-transforms, 9 ffts.  
do n=1,3
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nslabx
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

         if (n==1) then
            wy = - jm*Qhat(k,i,j+z_jmsign(j),3)
            vz =  - km*Qhat(k+z_kmsign(k),i,j,2)
            !rhs(k,i,j,1) = pi2*(wy - vz)
            p(k,i,j) = pi2*(wy - vz)
         endif
         
         if (n==2) then
            wx = - im*Qhat(k,i+z_imsign(i),j,3)
            uz =  - km*Qhat(k+z_kmsign(k),i,j,1)
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
! alternative:
! compute from Q_grid: vx,wx  uy,wy      8 ffts, 0 z-transforms 
!              Qhat:   vz,uz             6 ffts, 2 z-transforms
!




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
   vor = vor + rhsg(i,j,k,3)
   
   hel = hel + Q(i,j,k,1)*rhsg(i,j,k,1) + & 
        Q(i,j,k,2)*rhsg(i,j,k,2) + & 
        Q(i,j,k,3)*rhsg(i,j,k,3)  
   
   maxvor = max(maxvor,abs(rhsg(i,j,k,1)))
   maxvor = max(maxvor,abs(rhsg(i,j,k,2)))
   maxvor = max(maxvor,abs(rhsg(i,j,k,3)))
   
   !  velocity=(u,v,w)  vorticity=(a,b,c)=(wy-vz,uz-wx,vx-uy)
   !  v*(vx-uy) - w*(uz-wx) = (v vx - v uy + w wx) - w uz
   !  w*(wy-vz) - u*(vx-uy)
   !  u*(uz-wx) - v*(wy-vz)
   uu = ( Q(i,j,k,2)*rhsg(i,j,k,3) - Q(i,j,k,3)*rhsg(i,j,k,2) )
   vv = ( Q(i,j,k,3)*rhsg(i,j,k,1) - Q(i,j,k,1)*rhsg(i,j,k,3) )
   ww = ( Q(i,j,k,1)*rhsg(i,j,k,2) - Q(i,j,k,2)*rhsg(i,j,k,1) )
   
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
            rhs(k,i,j,1)=rhs(k,i,j,1) + xfac*Qhat(k,i,j,1)
            rhs(k,i,j,2)=rhs(k,i,j,2) + xfac*Qhat(k,i,j,2)
            rhs(k,i,j,3)=rhs(k,i,j,3) + xfac*Qhat(k,i,j,3)

! < u (uxx + uyy + uzz) > = < u-hat * (uxx-hat + uyy-hat + uzz-hat) >
!                         = < u-hat*u-hat*( im**2 + jm**2 + km**2)

            xfac = 2*2*2*xfac
            if (km==0) xfac=xfac/2
            if (jm==0) xfac=xfac/2
            if (im==0) xfac=xfac/2

            ke_diss = ke_diss + xfac*(Qhat(k,i,j,1)**2 + &
                                Qhat(k,i,j,2)**2 + &
                                Qhat(k,i,j,3)**2) 

         

      enddo
   enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! add in forcing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (forcing_type>0) call sforcing(rhs,Qhat,f_diss)

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
         uu= - im*p(k,i+z_imsign(i),j) 
         vv= - jm*p(k,i,j+z_jmsign(j))
         ww= - km*p(k+z_kmsign(k),i,j)
         
         rhs(k,i,j,1)=rhs(k,i,j,1) - uu 
         rhs(k,i,j,2)=rhs(k,i,j,2) - vv 
         rhs(k,i,j,3)=rhs(k,i,j,3) - ww 

         ! dealias           
         if ( ((abs(km)> g_nz/3) ) .or. &
              ((abs(jm)> g_ny/3) ) .or. &
              ((abs(im)> g_nx/3) ) )  then
            rhs(k,i,j,1)=0
            rhs(k,i,j,2)=0
            rhs(k,i,j,3)=0
         endif

      enddo
   enddo
enddo



if (compute_ints==1) then
   ints_timeDU=time
   ints(2)=f_diss    ! computed in spectral space - dont have to normalize
   ints(3)=ke_diss   ! computed in spectral space - dont have to normalize
   ints(4)=vor/g_nx/g_ny/g_nz
   ints(5)=hel/g_nx/g_ny/g_nz
   maxs(5)=maxvor
endif

call wallclock(tmx2)
tims(5)=tims(5)+(tmx2-tmx1)
end











