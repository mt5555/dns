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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! three methods for divergence free projection.  All require 18 total FFTs
! (1) grid space method using derivative operators.  Cant dealias without extra FFTs
! (2) spectral space method, can dealias for free
! (3) same as (2) but with loops are fused.
! does either (2) or (3) have better performance? 
!
! Note: method (2) and (3) depend on FFT data structure, and are thus 
!       in the fft_*_interface.F90 files.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine divfree_gridspace(u,p)
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
real*8          :: p(nx,ny,nz)

real*8 :: work(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: dummy(1)
real*8 :: alpha=0
real*8 :: beta=1


!local variables
integer i,j,k

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







subroutine z_fft3d(f,fout)
!
!  compute fft of f, return in fout.
!  f,fout can overlap in memory
!
use params
use fft_interface
use transpose
implicit none
real*8 f(nx,ny,nz)    ! input
real*8 fout(*)        ! output
real*8 work(nx,ny,nz) ! work array1
real*8 work2(nx,ny,nz) ! work array2
integer n1,n1d,n2,n2d,n3,n3d


call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d)
call fft1(work,n1,n1d,n2,n2d,n3,n3d)     
call transpose_from_x(work,work2,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_y(work2,work,n1,n1d,n2,n2d,n3,n3d)
call fft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(work,work2,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_z(work2,fout,n1,n1d,n2,n2d,n3,n3d)
call fft1(fout,n1,n1d,n2,n2d,n3,n3d)

end




subroutine z_ifft3d(fin,f)
!
!  compute inverse fft 3d of fin, return in f
!  fin and f can overlap in memory
!
use params
use fft_interface
use transpose
implicit none
real*8 fin(g_nz2,nslabx,ny_2dz)  ! input
real*8 f(nx,ny,nz)    ! output
real*8 work(nx,ny,nz) ! work array
real*8 work2(g_nz2,nslabx,ny_2dz) ! work array


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k

n1=g_nz
n1d=g_nz2   	
n2=nslabx
n2d=nslabx
n3=ny_2dz
n3d=ny_2dz

work2=fin
call ifft1(work2,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(work2,f,n1,n1d,n2,n2d,n3,n3d)

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
use fft_interface ! for pi2_squared
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
!  Compute FFT derivates.
!
!  input: px 
!  output:
!     if numder=1   return d/dx along first direction in px
!                   (and pxx is not accessed) 
!     if numder=2   return d2/dx2 along first direction in pxx
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_derivatives(px,pxx,numder,n1,n1d,n2,n2d,n3,n3d)
use params ! just needed for pi
use fft_interface
implicit none

real*8 px(n1d,n2d,n3d)
real*8 pxx(n1d,n2d,n3d)
integer numder,n1,n1d,n2,n2d,n3,n3d

integer i,j,k,m
real*8 temp
ASSERT("fft99_interface.F90: numder<=2",numder<=2)
ASSERT("fft99_interface.F90: numder>=1",numder>=1)

call fft1(px,n1,n1d,n2,n2d,n3,n3d)

if (numder>=2) then
   do k=1,n3
   do j=1,n2
      do i=1,n1
         m=(i-1)/2
         pxx(i,j,k) = -m*m * pi2_squared * px(i,j,k)
      enddo

      ! note: for i=n1, we are working with the cos(n1/2) mode
      ! but d/dx of this mode goes to sin(n1/2) = 0 on our grid, so we
      ! just take m=0 to make sure that dxx = dx dx

   enddo
   enddo
   call ifft1(pxx,n1,n1d,n2,n2d,n3,n3d)
endif

   do k=1,n3
   do j=1,n2
      ! note: for i=2, m=0, we are actually working with the cos(n1/2) mode
      ! but d/dx of this mode goes to sin(n1/2) = 0, so just take m=0 
      do i = 1, n1,2
         m=(i-1)/2
         temp =  pi2* m * px(i,j,k)
         px(i,j,k) = -pi2 *m * px(i+1,j,k)
         px(i+1,j,k) = temp
      enddo

   enddo
   enddo
   call ifft1(px,n1,n1d,n2,n2d,n3,n3d)

end subroutine





subroutine divfree(u,p)
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
real*8 p(nx,ny,nz)

!local variables
integer i,j,k
integer im,jm,km,i2,j2,k2
real*8 :: uu,vv,ww,xfac

do i=1,3
   call fft3d(u(1,1,1,i),p)
enddo

   do k=nz1,nz2
      km=(kmcord(k))
      if (km==g_nz/2) km=0
      do j=ny1,ny2
         jm=(jmcord(j))
         if (jm==g_ny/2) jm=0
         do i=nx1,nx2
            im=(imcord(i))
            if (im==g_nx/2) im=0

            ! compute the divergence
            p(i,j,k)= - im*u(i+imsign(i),j,k,1) &
                      - jm*u(i,j+jmsign(j),k,2) &
                      - km*u(i,j,k+kmsign(k),3)


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
            uu= - im*p(i+imsign(i),j,k)  ! cosine mode
            vv= - jm*p(i,j+jmsign(j),k)
            ww= - km*p(i,j,k+kmsign(k))

            u(i,j,k,1)=u(i,j,k,1) - uu
            u(i,j,k,2)=u(i,j,k,2) - vv
            u(i,j,k,3)=u(i,j,k,3) - ww

         enddo
      enddo
   enddo

do i=1,3
   if (dealias) call fft_filter_dealias(u(1,1,1,i))
   call ifft3d(u(1,1,1,i),p)
enddo
end subroutine







subroutine divfree_loopfused(u,p)
!
! make u divergence free
!    solve:  div(u) = laplacian(p)
!    then:   unew = u - grad(p)
!
!  u  input/output 
!  p is used as a work array    
! 
!
use params
implicit none
real*8 u(nx,ny,nz,3)
real*8 p(nx,ny,nz)

!local variables
integer i,j,k
integer im,jm,km,i0,i1,j0,j1,k0,k1,n
real*8 :: xfac
integer ii,jj,kk

real*8 pi0j0k0
real*8 pi1j0k0
real*8 pi0j1k0
real*8 pi1j1k0
real*8 pi0j0k1
real*8 pi1j0k1
real*8 pi0j1k1
real*8 pi1j1k1 

do i=1,3
   call fft3d(u(1,1,1,i),p)
enddo

p=0
   do k=nz1,nz2,2
      km=kmcord(k)
      if (km==g_nz/2) km=0
      do j=ny1,ny2,2
         jm=jmcord(j)
         if (jm==g_ny/2) jm=0
         do i=nx1,nx2,2
            im=imcord(i)
            if (im==g_nx/2) im=0

            ! compute the divergence
            ! compute laplacian inverse
            xfac= (im*im +km*km + jm*jm)
            if (xfac/=0) xfac = -1/xfac

            i0=i
            i1=i+1
            j0=j
            j1=j+1
            k0=k
            k1=k+1
            if (g_nz==1) k1=k ! ignore k, we are doing a 2d problem.  km=0

            !                  u_x                  v_y               w_z
            pi0j0k0 = (-im*u(i1,j0,k0,1)  -  jm*u(i0,j1,k0,2) - km*u(i0,j0,k1,3))*xfac
            pi1j0k0 = ( im*u(i0,j0,k0,1)  -  jm*u(i1,j1,k0,2) - km*u(i1,j0,k1,3))*xfac
            pi0j1k0 = (-im*u(i1,j1,k0,1)  +  jm*u(i0,j0,k0,2) - km*u(i0,j1,k1,3))*xfac
            pi1j1k0 = ( im*u(i0,j1,k0,1)  +  jm*u(i1,j0,k0,2) - km*u(i1,j1,k1,3))*xfac

            if (g_nz>1) then
            pi0j0k1 = (-im*u(i1,j0,k1,1)  -  jm*u(i0,j1,k1,2) + km*u(i0,j0,k0,3))*xfac
            pi1j0k1 = ( im*u(i0,j0,k1,1)  -  jm*u(i1,j1,k1,2) + km*u(i1,j0,k0,3))*xfac
            pi0j1k1 = (-im*u(i1,j1,k1,1)  +  jm*u(i0,j0,k1,2) + km*u(i0,j1,k0,3))*xfac
            pi1j1k1 = ( im*u(i0,j1,k1,1)  +  jm*u(i1,j0,k1,2) + km*u(i1,j1,k0,3))*xfac
            endif

            ! u = u - px
            u(i0,j0,k0,1) = u(i0,j0,k0,1) + im*pi1j0k0
            u(i1,j0,k0,1) = u(i1,j0,k0,1) - im*pi0j0k0
            u(i0,j1,k0,1) = u(i0,j1,k0,1) + im*pi1j1k0
            u(i1,j1,k0,1) = u(i1,j1,k0,1) - im*pi0j1k0
            if (g_nz>1) then
            u(i0,j0,k1,1) = u(i0,j0,k1,1) + im*pi1j0k1
            u(i1,j0,k1,1) = u(i1,j0,k1,1) - im*pi0j0k1
            u(i0,j1,k1,1) = u(i0,j1,k1,1) + im*pi1j1k1
            u(i1,j1,k1,1) = u(i1,j1,k1,1) - im*pi0j1k1
            endif

            ! v = v - py
            u(i0,j0,k0,2) = u(i0,j0,k0,2) + jm*pi0j1k0
            u(i1,j0,k0,2) = u(i1,j0,k0,2) + jm*pi1j1k0
            u(i0,j1,k0,2) = u(i0,j1,k0,2) - jm*pi0j0k0
            u(i1,j1,k0,2) = u(i1,j1,k0,2) - jm*pi1j0k0
            if (g_nz>1) then
            u(i0,j0,k1,2) = u(i0,j0,k1,2) + jm*pi0j1k1
            u(i1,j0,k1,2) = u(i1,j0,k1,2) + jm*pi1j1k1
            u(i0,j1,k1,2) = u(i0,j1,k1,2) - jm*pi0j0k1
            u(i1,j1,k1,2) = u(i1,j1,k1,2) - jm*pi1j0k1
            endif


            ! v = v - pz
            if (g_nz>1) then
            u(i0,j0,k0,3) = u(i0,j0,k0,3) + km*pi0j0k1
            u(i1,j0,k0,3) = u(i1,j0,k0,3) + km*pi1j0k1
            u(i0,j1,k0,3) = u(i0,j1,k0,3) + km*pi0j1k1
            u(i1,j1,k0,3) = u(i1,j1,k0,3) + km*pi1j1k1
            u(i0,j0,k1,3) = u(i0,j0,k1,3) - km*pi0j0k0
            u(i1,j0,k1,3) = u(i1,j0,k1,3) - km*pi1j0k0
            u(i0,j1,k1,3) = u(i0,j1,k1,3) - km*pi0j1k0
            u(i1,j1,k1,3) = u(i1,j1,k1,3) - km*pi1j1k0
            endif

            if (dealias) then
            if ( ((km>=g_nz/3) .and. (km>0)) .or. &
                 ((jm>=g_ny/3) .and. (jm>0)) .or. &
                 ((im>=g_nx/3) .and. (im>0)) )  then
               do n=1,3
                  u(i0,j0,k0,n) = 0
                  u(i1,j0,k0,n) = 0
                  u(i0,j1,k0,n) = 0
                  u(i1,j1,k0,n) = 0
                  u(i0,j0,k1,n) = 0
                  u(i1,j0,k1,n) = 0
                  u(i0,j1,k1,n) = 0
                  u(i1,j1,k1,n) = 0     
               enddo
            endif            
            ! dont forget the last cosine mode, stored in 2nd array position.
            ! it has its mode number set to zero for the above computations
            if (im==0) then
               do n=1,3
                  u(i1,j0,k0,n) = 0
                  u(i1,j1,k0,n) = 0
                  u(i1,j0,k1,n) = 0
                  u(i1,j1,k1,n) = 0
               enddo
            endif
            if (jm==0) then
               do n=1,3
                  u(i0,j1,k0,n) = 0
                  u(i1,j1,k0,n) = 0
                  u(i0,j1,k1,n) = 0
                  u(i1,j1,k1,n) = 0
               enddo	
            endif
            if (km==0 .and. g_nz>1) then
               do n=1,3
                  u(i0,j0,k1,n) = 0
                  u(i1,j0,k1,n) = 0
                  u(i0,j1,k1,n) = 0
                  u(i1,j1,k1,n) = 0
               enddo
            endif
            endif

         enddo
      enddo
   enddo


do n=1,3
   call ifft3d(u(1,1,1,n),p)
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
      km=abs(kmcord(k))
      do j=ny1,ny2
         jm=abs(jmcord(j))
         do i=nx1,nx2
            im=abs(imcord(i))

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
      km=abs(kmcord(k))
      do j=ny1,ny2
         jm=abs(jmcord(j))
         do i=nx1,nx2
            im=abs(imcord(i))

            if ( (km>g_nz/3)  .or. &
                 (jm>g_ny/3)  .or. &
                 (im>g_nx/3) )  then
               p(i,j,k)=0
            endif
         enddo
      enddo
   enddo

end subroutine

