#include "macros.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  compute derivative along index index=1,2 or 3
!  numder = 1  compute p_x, return in px.   (pxx is not accessed)
!  numder = 2  compute p_xx, return in pxx.  
!
!  Note: it is safe, but maybe bad Fortran, to have p=px or pxx.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine der(p,px,pxx,pt,numder,index)
use params
use fft_interface
use transpose
implicit none

!input:
integer numder,index
real*8 p(nx,ny,nz)    ! original data
real*8 pt(nx,ny,nz)   ! work array

!output:
real*8 pxx(nx,ny,nz)
real*8 px(nx,ny,nz)

integer n1,n1d,n2,n2d,n3,n3d


if (index==1) then

   call transpose_to_x(p,pt,n1,n1d,n2,n2d,n3,n3d)
   call fft_derivatives(pt,px,numder,n1,n1d,n2,n2d,n3,n3d)
   if (numder==2) then
      call transpose_from_x(px,pxx,n1,n1d,n2,n2d,n3,n3d)
   endif
   call transpose_from_x(pt,px,n1,n1d,n2,n2d,n3,n3d)


else if (index==2) then

   call transpose_to_y(p,pt,n1,n1d,n2,n2d,n3,n3d)
   ! 1st derivative returned in pt, 2nd derivative returned in px
   call fft_derivatives(pt,px,numder,n1,n1d,n2,n2d,n3,n3d)
   if (numder==2) then
      call transpose_from_y(px,pxx,n1,n1d,n2,n2d,n3,n3d)
   endif
   call transpose_from_y(pt,px,n1,n1d,n2,n2d,n3,n3d)

else if (index==3) then
   call transpose_to_z(p,pt,n1,n1d,n2,n2d,n3,n3d)
   call fft_derivatives(pt,px,numder,n1,n1d,n2,n2d,n3,n3d)
   if (numder==2) then
      call transpose_from_z(px,pxx,n1,n1d,n2,n2d,n3,n3d)
   endif
   call transpose_from_z(pt,px,n1,n1d,n2,n2d,n3,n3d)

endif


end subroutine



subroutine divfree(u,p,work,work2)
!
! make u divergence free
!    solve:  div(u) = laplacian(p)
!    then:   unew = u - grad(p)
!    
! 
!
use params
use fft_interface
implicit none
real*8 u(nx,ny,nz,3)
real*8 work(nx,ny,nz)
real*8 work2(nx,ny,nz)
real*8 p(nx,ny,nz)

!local variables
real*8 :: dummy(1)
real*8 :: alpha=0
real*8 :: beta=1
integer i,j,k

!
! two methods.  both require 18 total FFT's but second method has
! the advantage that we can also dealias at no additional cost
!
#if 0
call divergence(p,u,work,work2)

! solve laplacian(p)=div(u)
call poisson(p,work,alpha,beta)

! compute u=u-grad(p)
do i=1,3
   call der(p,work,dummy,work2,1,i)
   u(:,:,:,i) = u(:,:,:,i) - work
enddo

if (dealias) then
   do i=1,3
      call fft3d(u(1,1,1,i),work)
      call fft_filter_dealias(u(1,1,1,i))
      call ifft3d(u(1,1,1,i),work)
   enddo
endif


#else
integer im,jm,km,i2,j2,k2
real*8 :: uu,vv,ww,xfac


ASSERT("divfree(): nslabx must be even ",mod(nslabx,2)==0)
ASSERT("divfree(): nslaby must be even ",mod(nslaby,2)==0)
ASSERT("divfree(): nslabz must be even ",(mod(nslabz,2)==0 .or. nslabz==1))

do i=1,3
   call fft3d(u(1,1,1,i),work)
enddo

   do k=nz1,nz2
      km=kmcord(k)
      if (km==g_nz/2) km=0
      do j=ny1,ny2
         jm=jmcord(j)
         if (jm==g_ny/2) jm=0
         do i=nx1,nx2
            im=imcord(i)
            if (im==g_nx/2) im=0

            ! compute the divergence
            p(i,j,k)=0
            if (mod(i-nx1+1,2)==1) then
               p(i,j,k)=p(i,j,k) - im*u(i+1,j,k,1)
            else
               p(i,j,k)=p(i,j,k) + im*u(i-1,j,k,1)
            endif

            if (mod(j-ny1+1,2)==1) then
               p(i,j,k)=p(i,j,k) - jm*u(i,j+1,k,2)
            else
               p(i,j,k)=p(i,j,k) + jm*u(i,j-1,k,2)
            endif

            if (g_nz>1) then
            if (mod(k-nz1+1,2)==1) then
               p(i,j,k)=p(i,j,k) - km*u(i,j,k+1,3)
            else
               p(i,j,k)=p(i,j,k) + km*u(i,j,k-1,3)
            endif
            endif

            ! compute laplacian inverse
            xfac= (im*im +km*km + jm*jm)
            if (xfac/=0) xfac = -1/xfac
            p(i,j,k)=xfac*p(i,j,k)


         enddo
      enddo
   enddo

   do k=nz1,nz2
      km=kmcord(k)
      if (km==g_nz/2) km=0
      do j=ny1,ny2
         jm=jmcord(j)
         if (jm==g_ny/2) jm=0
         do i=nx1,nx2
            im=imcord(i)
            if (im==g_nx/2) im=0

            ! compute gradient  dp/dx
            if (mod(i-nx1+1,2)==1) then
               uu= - im*p(i+1,j,k)

            else
               uu= + im*p(i-1,j,k)
            endif
            if (mod(j-ny1+1,2)==1) then
               vv= - jm*p(i,j+1,k)
            else
               vv= + jm*p(i,j-1,k)
            endif
            if (mod(k-nz1+1,2)==1) then
               ww= - km*p(i,j,k+1)
            else
               ww= + km*p(i,j,k-1)
            endif

            u(i,j,k,1)=u(i,j,k,1) - uu
            u(i,j,k,2)=u(i,j,k,2) - vv
            u(i,j,k,3)=u(i,j,k,3) - ww

         enddo
      enddo
   enddo


do i=1,3
   if (dealias) call fft_filter_dealias(u(1,1,1,i))
   call ifft3d(u(1,1,1,i),work)
enddo
#endif





end subroutine


subroutine vorticity(vor,u,d1,work)
use params
use fft_interface
use transpose
implicit none
real*8 u(nx,ny,nz,3)    ! input
real*8 vor(nx,ny,nz,3)    ! output
real*8 d1(nx,ny,nz) 
real*8 work(nx,ny,nz) 

! local variables
integer i
real*8 dummy(1)

vor=0
do i=1,3

   ! compute u_x, u_xx
   call der(u(1,1,1,i),d1,dummy,work,DX_ONLY,1)
   if (i==3) vor(:,:,:,2) = vor(:,:,:,2) - d1
   if (i==2) vor(:,:,:,3) = vor(:,:,:,3) + d1

   ! compute u_y, u_yy
   call der(u(1,1,1,i),d1,dummy,work,DX_ONLY,2)
   if (i==3) vor(:,:,:,1) = vor(:,:,:,1) + d1
   if (i==1) vor(:,:,:,3) = vor(:,:,:,3) -d1

   ! compute u_z, u_zz
   call der(u(1,1,1,i),d1,dummy,work,DX_ONLY,3)
   if (i==2) vor(:,:,:,1) = vor(:,:,:,1) -d1
   if (i==1) vor(:,:,:,2) = vor(:,:,:,2) +d1

enddo


end subroutine




subroutine divergence(div,u,work1,work2)
use params
use fft_interface
use transpose
implicit none
real*8 u(nx,ny,nz,3)    ! input
real*8 div(nx,ny,nz)    ! output
real*8 work1(nx,ny,nz) ! wk array
real*8 work2(nx,ny,nz) ! wk array

! local variables
integer i
real*8 dummy(1)

i=1
call der(u(1,1,1,i),div,dummy,work2,DX_ONLY,i)

i=2
call der(u(1,1,1,i),work1,dummy,work2,DX_ONLY,i)
div = div+work1

i=3
call der(u(1,1,1,i),work1,dummy,work2,DX_ONLY,i)
div = div+work1


end subroutine









subroutine poisson(f,work,alpha,beta)
!
!  solve laplacian(p) = f
!  input:  f 
!  ouput:  f   will be overwritten with the solution p
!
use params
use fft_interface
use transpose
implicit none
real*8 f(nx,ny,nz)    ! input/output
real*8 work(nx,ny,nz) ! work array
real*8 :: alpha
real*8 :: beta


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k


call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d) 
call fft1(work,n1,n1d,n2,n2d,n3,n3d)     
call transpose_from_x(work,f,n1,n1d,n2,n2d,n3,n3d) 


call transpose_to_y(f,work,n1,n1d,n2,n2d,n3,n3d)  ! x,y,z -> y,x,z
call fft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(work,f,n1,n1d,n2,n2d,n3,n3d) 

call transpose_to_z(f,work,n1,n1d,n2,n2d,n3,n3d)  
call fft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(work,f,n1,n1d,n2,n2d,n3,n3d)  

! solve [alpha + beta*Laplacian] p = f.  f overwritten with output  p
call fft_laplace_inverse(f,alpha,beta)

call transpose_to_z(f,work,n1,n1d,n2,n2d,n3,n3d)       
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(work,f,n1,n1d,n2,n2d,n3,n3d)       

call transpose_to_y(f,work,n1,n1d,n2,n2d,n3,n3d)       
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(work,f,n1,n1d,n2,n2d,n3,n3d)         ! y,x,z -> x,y,z

call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d) 
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_x(work,f,n1,n1d,n2,n2d,n3,n3d )



end




subroutine fft3d(f,work)
!
!  compute the spectrum, ouput in f
!
use params
use fft_interface
use transpose
implicit none
real*8 f(nx,ny,nz)    ! input/output
real*8 work(nx,ny,nz) ! work array
integer n1,n1d,n2,n2d,n3,n3d


call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d)
call fft1(work,n1,n1d,n2,n2d,n3,n3d)     
call transpose_from_x(work,f,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_y(f,work,n1,n1d,n2,n2d,n3,n3d)
call fft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(work,f,n1,n1d,n2,n2d,n3,n3d)


call transpose_to_z(f,work,n1,n1d,n2,n2d,n3,n3d)
call fft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(work,f,n1,n1d,n2,n2d,n3,n3d)


end




subroutine ifft3d(f,work)
!
!  compute inverse fft 3d of f, return in f
!
use params
use fft_interface
use transpose
implicit none
real*8 f(nx,ny,nz)    ! input/output
real*8 work(nx,ny,nz) ! work array


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k

call transpose_to_z(f,work,n1,n1d,n2,n2d,n3,n3d)
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(work,f,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_y(f,work,n1,n1d,n2,n2d,n3,n3d)
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(work,f,n1,n1d,n2,n2d,n3,n3d)


call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d)
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_x(work,f,n1,n1d,n2,n2d,n3,n3d)


end







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Solve  [alpha + beta*Laplacian] p = rhs
!
! on input,  p = fourier coefficients of rhs
! on output, p = fourier coefficients of solution
!
! highest mode tweaked so that laplacian = div grad
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_laplace_inverse(p,alpha,beta)
use params
implicit none
real*8 p(nx,ny,nz)
real*8 alpha,beta

integer i,j,k,im,jm,km
real*8 xfac

   do k=nz1,nz2
      km=kmcord(k)
      if (km==g_nz/2) km=0
      do j=ny1,ny2
         jm=jmcord(j)
         if (jm==g_ny/2) jm=0
         do i=nx1,nx2
            im=imcord(i)
            if (im==g_nx/2) im=0
            xfac= alpha + beta*(-im*im -km*km - jm*jm)*pi2_squared      
            if (xfac/=0) xfac = 1/xfac
            p(i,j,k)=p(i,j,k)*xfac
         enddo
      enddo
   enddo

end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Filter out the highest cosine mode
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_filter_last(p)
use params
implicit none
real*8 p(nx,ny,nz)

integer i,j,k,im,jm,km
real*8 xfac

   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)

            if ( ((km==g_nz/2) .and. (km>0)) .or. &
                 ((jm==g_ny/2) .and. (jm>0)) .or. &
                 ((im==g_nx/2) .and. (im>0)) )  then
               p(i,j,k)=0
            endif
         enddo
      enddo
   enddo

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Filter out last 1/3 of spectrum (de-alias quadratic terms)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_filter_dealias(p)
use params
implicit none
real*8 p(nx,ny,nz)

integer i,j,k,im,jm,km
real*8 xfac

   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)

            if ( ((km>=g_nz/3) .and. (km>0)) .or. &
                 ((jm>=g_ny/3) .and. (jm>0)) .or. &
                 ((im>=g_nx/3) .and. (im>0)) )  then
               p(i,j,k)=0
            endif
         enddo
      enddo
   enddo

end subroutine

