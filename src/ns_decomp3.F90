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

! local variables
real*8 :: Q_tmp(nx,ny,nz,n_var)
real*8 :: ke_old,time_old,vel
real*8 :: Q_old(nx,ny,nz,n_var)
integer i,j,k,n,ierr
logical,save :: firstcall=.true.



if (firstcall) then
   firstcall=.false.
   Q=Q_grid
   do n=1,3
      call fft3d(Q(1,1,1,n),rhs)  ! use rhs as a work array
   enddo
   if (.not. dealias) then
      call abort("Error: using ns3dspectral model, which must be run dealiased")
   endif
   if (forcing_type==1) then
      call abort("Error: SPEC_SLOW model cannot use forcing_type=1")
   endif
endif







! stage 1
call ns3D(rhs,Q,Q_grid,time,1)


do n=1,3
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
call ns3D(rhs,Q_tmp,Q_grid,time+delt/2.0,0)


do n=1,3
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
call ns3D(rhs,Q_tmp,Q_grid,time+delt/2.0,0)

do n=1,3
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
call ns3D(rhs,Q_tmp,Q_grid,time+delt,0)


do n=1,3
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
ints(1)=ints(1)/g_nx/g_ny/g_nz



end subroutine rk4  







subroutine ns3d(rhs,Qhat,Q,time,compute_ints)
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
real*8 Qhat(nx,ny,nz,n_var)
real*8 Q(nx,ny,nz,n_var)
real*8 time
integer compute_ints

! output
real*8 rhs(nx,ny,nz,n_var)


!local
real*8 p(nx,ny,nz)
real*8 xfac,tmx1,tmx2
real*8 ux,uy,uz,wx,wy,wz,vx,vy,vz,uu,vv,ww
integer n,i,j,k,im,km,jm
integer n1,n1d,n2,n2d,n3,n3d
real*8 :: ke_diss,vor,hel,maxvor,f_diss


real*8 tmp(g_nz2,nslabx,ny_2dz,n_var)
real*8 tmp2(g_nz2,nslabx,ny_2dz,n_var)





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


do k=nz1,nz2
   km=kmcord(k)
   do j=ny1,ny2
      jm=jmcord(j)
      do i=nx1,nx2
         im=imcord(i)
         
         ! u_x term
         !ux = - im*Qhat(i+imsign(i),j,k,1)
         vx = - im*Qhat(i+imsign(i),j,k,2)
         wx = - im*Qhat(i+imsign(i),j,k,3)
         
         
         uy = - jm*Qhat(i,j+jmsign(j),k,1)
         !vy = - jm*Qhat(i,j+jmsign(j),k,2)
         wy = - jm*Qhat(i,j+jmsign(j),k,3)
         
         
         uz =  - km*Qhat(i,j,k+kmsign(k),1)
         vz =  - km*Qhat(i,j,k+kmsign(k),2)
         !wz =  - km*Qhat(i,j,k+kmsign(k),3)
         
         rhs(i,j,k,1) = pi2*(wy - vz)
         rhs(i,j,k,2) = pi2*(uz - wx)
         rhs(i,j,k,3) = pi2*(vx - uy)

      enddo
   enddo
enddo
! 3 forward&back z-transforms, 9 ffts.  
do n=1,3
   call ifft3d(rhs(1,1,1,n),p,n1,n1d,n2,n2d,n3,n3d)
enddo






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
      vor = vor + rhs(i,j,k,3)

      hel = hel + Q(i,j,k,1)*rhs(i,j,k,1) + & 
                  Q(i,j,k,2)*rhs(i,j,k,2) + & 
                  Q(i,j,k,3)*rhs(i,j,k,3)  

      maxvor = max(maxvor,abs(rhs(i,j,k,1)))
      maxvor = max(maxvor,abs(rhs(i,j,k,2)))
      maxvor = max(maxvor,abs(rhs(i,j,k,3)))

      uu = ( Q(i,j,k,2)*rhs(i,j,k,3) - Q(i,j,k,3)*rhs(i,j,k,2) )
      vv = ( Q(i,j,k,3)*rhs(i,j,k,1) - Q(i,j,k,1)*rhs(i,j,k,3) )
      ww = ( Q(i,j,k,1)*rhs(i,j,k,2) - Q(i,j,k,2)*rhs(i,j,k,1) )

      rhs(i,j,k,1) = uu
      rhs(i,j,k,2) = vv
      rhs(i,j,k,3) = ww
   enddo
   enddo
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! back to spectral space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do n=1,3
   call fft3d(rhs(1,1,1,n),p)
enddo

ke_diss = 0

!  diffusion
   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)

            xfac=-mu*(im*im + jm*jm + km*km)*pi2_squared
            rhs(i,j,k,1)=rhs(i,j,k,1) + xfac*Qhat(i,j,k,1)
            rhs(i,j,k,2)=rhs(i,j,k,2) + xfac*Qhat(i,j,k,2)
            rhs(i,j,k,3)=rhs(i,j,k,3) + xfac*Qhat(i,j,k,3)

! < u (uxx + uyy + uzz) > = < u-hat * (uxx-hat + uyy-hat + uzz-hat) >
!                         = < u-hat*u-hat*( im**2 + jm**2 + km**2)

            xfac = 2*2*2*xfac
            if (kmcord(k)==0) xfac=xfac/2
            if (jmcord(j)==0) xfac=xfac/2
            if (imcord(i)==0) xfac=xfac/2

            ke_diss = ke_diss + xfac*Qhat(i,j,k,1)**2 + &
                                xfac*Qhat(i,j,k,2)**2 + &
                                xfac*Qhat(i,j,k,3)**2 

         enddo
      enddo
   enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! add in forcing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (forcing_type>0) call sforcing(rhs,Qhat,f_diss)



!  make rhs div-free
   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)

            ! compute the divergence
            p(i,j,k)= - im*rhs(i+imsign(i),j,k,1) &
                      - jm*rhs(i,j+jmsign(j),k,2) &
                      - km*rhs(i,j,k+kmsign(k),3)


            ! compute laplacian inverse
            xfac= (im*im +km*km + jm*jm)
            if (xfac/=0) xfac = -1/xfac
            p(i,j,k)=xfac*p(i,j,k)


         enddo
      enddo
   enddo

   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)

            ! compute gradient  dp/dx
            uu= - im*p(i+imsign(i),j,k) 
            vv= - jm*p(i,j+jmsign(j),k)
            ww= - km*p(i,j,k+kmsign(k))

            rhs(i,j,k,1)=rhs(i,j,k,1) - uu 
            rhs(i,j,k,2)=rhs(i,j,k,2) - vv 
            rhs(i,j,k,3)=rhs(i,j,k,3) - ww 


            ! dealias           
            if ( ((abs(km)> g_nz/3) .and. (km/=0)) .or. &
                 ((abs(jm)> g_ny/3) .and. (jm/=0)) .or. &
                 ((abs(im)> g_nx/3) .and. (im/=0)) )  then
               rhs(i,j,k,1)=0
               rhs(i,j,k,2)=0
               rhs(i,j,k,3)=0
            endif


         enddo
      enddo
   enddo



if (compute_ints==1) then
   ints_timeDU=time
   ints(2)=ke_diss
   ints(3)=f_diss
   ints(4)=vor
   ints(5)=hel
   maxs(5)=maxvor
endif

call wallclock(tmx2)
tims(5)=tims(5)+(tmx2-tmx1)
end











