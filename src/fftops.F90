#include "macros.h"

! used to try out 2nd order scheme:
#undef CENTER2H
! used to try out 2nd order, but 2 delx scheme:
#undef CENTER4H        

! use DX applied twice for DXX in der()  
! this effects diffusion operator and iterative solves
! (because they call Helmholtz_periodic() which calls der())
#undef DX4DX4
!
! use DX applied twice for DXX in HELMHOLTZ_INVERSE
#define HINV_DX4DX4
!
! Notes on DX4DX4 settings.  These effect the routines used in
! periodic cases only:
!
!      helmholtz_periodic_inv()
!      der()
!      helmholtz_periodic()       (it uses d/dxx from der()
!
!
! Fourier methods, NS_UVW or NS_PSIVOR:  
!      direct solve:  DX4DX4 ignored, 
!                     NS_UVW requiires HINV_DX4DX4 so div(grad)=Helmholtz_periodic_inv
!                     NS_PSIVOR: both values of HINV_DX4DX4 will work
!      iterative solve, no preconditioner:
!              (only used for testing!)  HINV_DX4DX4 ignored
!              NS_UVW requires DX4DX4 so that div(grad)=Helmholtz_periodic
!              NS_PSIVOR: both values of DX4DX4 should work.   

! Shallow Water:
!      uses iterative solve with helmholtz_hform_periodic()
!      which computes using calls to der() for first derivatives only
!      so DX4DX4 doesn't matter.  
!      Predonditioner uses helmholtz_periodic_inv(), but should work
!      with either seting of HINV_DX4DX4.  
!
! 4th order methods:   
!   UVW form:  LAPLACE must be the same as div(grad)
!              (which is 2 successive calls to der() )
!
!              direct solve:  
!                    HINV_DX4DX4 set
!                    DX4DX4 doesn't matter
!              iterative solve, no preconditioner:
!                    DX4DX4       set
!                    HINV_DX4DX4  doesn't matter (not used)
!              iterative solve, with preconditioner:
!                    both set
!
!
!   PSI-VOR    direct solve: 
!                    HINV_DX4DX4 both work, different answers
!                    DX4DX4 doesn't matter 
!              iterative solve, no preconditioner:
!                    DX4DX4       both work
!                    HINV_DX4DX4  doesn't matter
!              iterative solve, with preconditioner:
!                    both set, or both unset
!   
!   
! For non-periodic cases, we are using:
! helmholtz_dirichlet_inv()
! helmholtz_dirichlet()
! which are hard coded to use 3 point stencils for d/dxx
!                



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

if (bdy_y1/=PERIODIC) then
   call abort('der() can only handle periodic boundaries')
endif
if (bdy_y1/=PERIODIC) then
   call abort('der() can only handle periodic boundaries')
endif
if (bdy_z1/=PERIODIC) then
   call abort('der() can only handle periodic boundaries')
endif

if (numerical_method==FORTH_ORDER) then

if (index==1) then

   call transpose_to_x(p,pt,n1,n1d,n2,n2d,n3,n3d)
   call fd_derivatives(pt,px,numder,n1,n1d,n2,n2d,n3,n3d,delx)
   if (numder==2) then
      call transpose_from_x(px,pxx,n1,n1d,n2,n2d,n3,n3d)
   endif
   call transpose_from_x(pt,px,n1,n1d,n2,n2d,n3,n3d)


else if (index==2) then

   call transpose_to_y(p,pt,n1,n1d,n2,n2d,n3,n3d)
   ! 1st derivative returned in pt, 2nd derivative returned in px
   call fd_derivatives(pt,px,numder,n1,n1d,n2,n2d,n3,n3d,dely)
   if (numder==2) then
      call transpose_from_y(px,pxx,n1,n1d,n2,n2d,n3,n3d)
   endif
   call transpose_from_y(pt,px,n1,n1d,n2,n2d,n3,n3d)

else if (index==3) then
   call transpose_to_z(p,pt,n1,n1d,n2,n2d,n3,n3d)
   call fd_derivatives(pt,px,numder,n1,n1d,n2,n2d,n3,n3d,delz)
   if (numder==2) then
      call transpose_from_z(px,pxx,n1,n1d,n2,n2d,n3,n3d)
   endif
   call transpose_from_z(pt,px,n1,n1d,n2,n2d,n3,n3d)

endif


else ! FFT method

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
subroutine divfree_gridspace(u,p,work,work2)
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
real*8 :: u(nx,ny,nz,3)
real*8 :: p(nx,ny,nz)
real*8 :: work(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

!local
real*8 :: dummy(1),tol
real*8 :: alpha=0
real*8 :: beta=1
integer i,j,k,n
external helmholtz_periodic,helmholtz_dirichlet


! solve laplacian(p)=div(u)

if (bdy_x1==PERIODIC .and. bdy_y1==PERIODIC .and. bdy_z1==PERIODIC) then
   call divergence(p,u,work,work2)
   call helmholtz_periodic_inv(p,work,alpha,beta)

   !work=p  ! RHS
   !p=0  ! initial guess
   !tol=1e-10
   !call cgsolver(p,work,alpha,beta,tol,work2,helmholtz_periodic,.false.)

else
   stop 'divfree_gridspace: only supports periodic case'
   call divergence(work,u,p,work2)
   p=0  ! initial guess
   tol=1e-10
   call cgsolver(p,work,alpha,beta,tol,work2,helmholtz_dirichlet,.false.)
endif




! compute u=u-grad(p)
do n=1,3
   call der(p,work,dummy,work2,1,n)
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      u(i,j,k,n)=u(i,j,k,n)-work(i,j,k)
   enddo
   enddo
   enddo
enddo

if (dealias) then
   call dealias_gridspace(u,work)
endif

end subroutine





subroutine dealias_gridspace(u,work)
use params
use fft_interface
implicit none
real*8 :: u(nx,ny,nz,3)
real*8 :: work(nx,ny,nz)
integer :: i

do i=1,3
   call fft3d(u(1,1,1,i),work)
   call fft_filter_dealias(u(1,1,1,i))
   call ifft3d(u(1,1,1,i),work)
enddo

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
integer i,j,k,n
real*8 dummy(1)

vor=0
do n=1,ndim

   ! compute u_x, u_xx
   call der(u(1,1,1,n),d1,dummy,work,DX_ONLY,1)
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      if (n==3) vor(i,j,k,2) = vor(i,j,k,2) - d1(i,j,k)
      if (n==2) vor(i,j,k,3) = vor(i,j,k,3) + d1(i,j,k)
   enddo
   enddo
   enddo

   ! compute u_y, u_yy
   call der(u(1,1,1,n),d1,dummy,work,DX_ONLY,2)
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      if (n==3) vor(i,j,k,1) = vor(i,j,k,1) + d1(i,j,k)
      if (n==1) vor(i,j,k,3) = vor(i,j,k,3) -d1(i,j,k)
   enddo
   enddo
   enddo

   if (ndim==3) then
   ! compute u_z, u_zz
   call der(u(1,1,1,n),d1,dummy,work,DX_ONLY,3)
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      if (n==2) vor(i,j,k,1) = vor(i,j,k,1) -d1(i,j,k)
      if (n==1) vor(i,j,k,2) = vor(i,j,k,2) +d1(i,j,k)
   enddo
   enddo
   enddo
   endif

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
integer i,j,k,n
real*8 dummy(1)

n=1
call der(u(1,1,1,n),div,dummy,work2,DX_ONLY,n)

n=2
call der(u(1,1,1,n),work1,dummy,work2,DX_ONLY,n)
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   div(i,j,k) = div(i,j,k)+work1(i,j,k)
enddo
enddo
enddo

if (ndim==3) then
n=3
call der(u(1,1,1,n),work1,dummy,work2,DX_ONLY,n)
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   div(i,j,k) = div(i,j,k)+work1(i,j,k)
enddo
enddo
enddo
endif


end subroutine








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







subroutine z_fft3d_trashinput(f,fout,work)
!
!  compute fft of f, return in fout.
!  f,fout can overlap in memory
!  data in f is ovewritten
!
use params
use fft_interface
use transpose

use ghost

implicit none
real*8 f(nx,ny,nz)    ! input
real*8 fout(g_nz2,nslabx,ny_2dz)  ! output
real*8 work(nx,ny,nz) ! work array1
integer n1,n1d,n2,n2d,n3,n3d



call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d)
call fft1(work,n1,n1d,n2,n2d,n3,n3d)     
call transpose_from_x(work,f,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_y(f,work,n1,n1d,n2,n2d,n3,n3d)
call fft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(work,f,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_z(f,fout,n1,n1d,n2,n2d,n3,n3d)
call fft1(fout,n1,n1d,n2,n2d,n3,n3d)


end




subroutine z_ifft3d(fin,f,work)
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
! true size must be nx,ny,nz:
real*8 work(g_nz2,nslabx,ny_2dz) ! work array


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k

n1=g_nz
n1d=g_nz2   	
n2=nslabx
n2d=nslabx
n3=ny_2dz
n3d=ny_2dz

work=fin
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
use fft_interface ! for pi2_squared
implicit none
real*8 p(nx,ny,nz)
real*8 alpha,beta

!local
integer i,j,k,im,jm,km
real*8 xfac,xm,ym,zm

if (numerical_method==FORTH_ORDER) then
do k=nz1,nz2
   do j=ny1,ny2
      do i=nx1,nx2

         ! u(x+h)-u(x-h)    ->   2i sin(k*pi2*h)
         ! u(x+2h)-u(x-2h)  ->   2i sin(k*pi2*2*h)
         ! applied twice: ->   -4 sin(k*pi2*h)^2

#ifdef CENTER2H
         xm=2*sin(imcord(i)*pi2*delx)/(2*delx)
         ym=2*sin(jmcord(j)*pi2*dely)/(2*dely)
         zm=2*sin(kmcord(k)*pi2*delz)/(2*delz)
         xfac=-(xm*xm + ym*ym  +zm*zm)
#elif (defined CENTER4H) 
         xm=2*sin(imcord(i)*pi2*2*delx)/(4*delx)
         ym=2*sin(jmcord(j)*pi2*2*dely)/(4*dely)
         zm=2*sin(kmcord(k)*pi2*2*delz)/(4*delz)
         xfac=-(xm*xm + ym*ym  +zm*zm)
#elif (defined HINV_DX4DX4)
         xm=4*sin(imcord(i)*pi2*delx)/3
         ym=4*sin(jmcord(j)*pi2*dely)/3
         zm=4*sin(kmcord(k)*pi2*delz)/3
         xm=xm - sin(imcord(i)*pi2*2*delx)/6
         ym=ym - sin(jmcord(j)*pi2*2*dely)/6
         zm=zm - sin(kmcord(k)*pi2*2*delz)/6
         xm=xm/delx
         ym=ym/dely
         zm=zm/delz
         xfac=-(xm*xm + ym*ym  +zm*zm)
#else
         ! -u(x+2h)+16(x+h)-30u(x)+16u(x-h)-u(x-2h)/(12*h*h)
         ! u(x+h)+u(x-h)   ->  2 cos(k*pi2*h)
         ! u(x+2h)+u(x-2h) ->  2 cos(k*pi2*2h)
         xm=-30
         ym=-30
         zm=-30
         xm=xm + 2*16*cos(imcord(i)*pi2*delx)
         ym=ym + 2*16*cos(jmcord(j)*pi2*dely)
         zm=zm + 2*16*cos(kmcord(k)*pi2*delz)
         xm=xm - 2*cos(imcord(i)*pi2*2*delx)
         ym=ym - 2*cos(jmcord(j)*pi2*2*dely)
         zm=zm - 2*cos(kmcord(k)*pi2*2*delz)
         xm=xm/(12*delx*delx)
         ym=ym/(12*dely*dely)
         zm=zm/(12*delz*delz)
         xfac=xm+ym+zm
#endif

         if (abs(xfac)<1e-12) xfac=0
         xfac= alpha + beta*xfac
         if (xfac/=0) xfac = 1/xfac
         p(i,j,k)=p(i,j,k)*xfac
      enddo
   enddo
enddo

else

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
endif


end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Compute 4th order derivates.
!  This routine will be obsolete when we get the ghost cells working,
!  and should be deleted
!
!  input: px 
!  output:
!     if numder=1   return d/dx along first direction in px
!                   (and pxx is not accessed) 
!     if numder=2   return d2/dx2 along first direction in pxx
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fd_derivatives(px,pxx,numder,n1,n1d,n2,n2d,n3,n3d,h)

implicit none

integer numder,n1,n1d,n2,n2d,n3,n3d
real*8 :: px(n1d,n2d,n3d)
real*8 :: pxx(n1d,n2d,n3d)
real*8 :: work(n1d)
real*8 :: h

!local
integer i,j,k,i0,i1,i2,i3

if (n1==1) then
   px=0
   pxx=0
   return
endif



do k=1,n3
do j=1,n2

   do i=1,n1
      work(i)=px(i,j,k)
   enddo
   i0=n1-1
   i1=n1
   i2=2
   i3=3
   do i=1,n1
#ifdef CENTER2H
      px(i,j,k)= (work(i2)-work(i1))/(2*h)
#elif (defined CENTER4H)
      px(i,j,k)= (work(i3)-work(i0))/(4*h)
#else
      px(i,j,k)= (2*(work(i2)-work(i1))/3 - (work(i3)-work(i0))/12 )/h
#endif
      i0=i1
      i1=i
      i2=i3
      i3=i3+1
      if (i3>n1) i3=i3-n1
   enddo
   if (numder>=2) then
      i0=n1-1
      i1=n1
      i2=2
      i3=3
      do i=1,n1
#ifdef CENTER2H
         pxx(i,j,k)= (px(i2,j,k)-px(i1,j,k))/(2*h)
#elif (defined CENTER4H)
         pxx(i,j,k)= (px(i3,j,k)-px(i0,j,k))/(4*h)
#elif (defined DX4DX4)
         pxx(i,j,k)= (2*(px(i2,j,k)-px(i1,j,k))/3 - (px(i3,j,k)-px(i0,j,k))/12 )/h
#else
         pxx(i,j,k)= (-work(i0)+16*work(i1)-30*work(i)+16*work(i2)-work(i3) )&
               /(12*h*h)
#endif
         i0=i1
         i1=i
         i2=i3
         i3=i3+1
         if (i3>n1) i3=i3-n1
      enddo
   endif
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

integer numder,n1,n1d,n2,n2d,n3,n3d
real*8 px(n1d,n2d,n3d)
real*8 pxx(n1d,n2d,n3d)


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




#if 0
example code, cant handle non fft methods 

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
#endif





#if 0
example code, cant handle non fft methods 
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

            if (ndim==3) then
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
            if (ndim==3) then
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
            if (ndim==3) then
            u(i0,j0,k1,2) = u(i0,j0,k1,2) + jm*pi0j1k1
            u(i1,j0,k1,2) = u(i1,j0,k1,2) + jm*pi1j1k1
            u(i0,j1,k1,2) = u(i0,j1,k1,2) - jm*pi0j0k1
            u(i1,j1,k1,2) = u(i1,j1,k1,2) - jm*pi1j0k1
            endif


            ! v = v - pz
            if (ndim==3) then
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
            if (km==0 .and. ndim==3) then
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
#endif







subroutine helmholtz_periodic(f,lf,alpha,beta,work)
!
!  input: f
!  output: lf
!
!     lf = [alpha + beta*laplacian](f)
!
use params
use fft_interface
use transpose
implicit none
real*8 f(nx,ny,nz)    ! input
real*8 lf(nx,ny,nz)    ! output
real*8 work(nx,ny,nz)    ! work array
real*8 :: alpha
real*8 :: beta

!local
real*8 gradf(nx,ny,nz,2) ! work array
real*8 fxx(nx,ny,nz)     ! work array
real*8 dummy
integer n

! linear helmholtz operator
lf=alpha*f
if (beta==0) return
do n=1,3
   call der(f,gradf,fxx,work,DX_AND_DXX,n)
   lf=lf+beta*fxx
enddo

end subroutine





#define COMPACT
subroutine helmholtz_dirichlet_setup(f,p,work,setbdy)
!
! for compact: replace f with:   f + h**2/12 (fxx + fyy)
! in the interior!
!
! THEN, if setbdy==1 set boundary values in f to those given in p.
!
! be sure to have called ghost_update(f) before calling this routine!
!
use params
use ghost
implicit none
real*8 f(nx,ny,nz)    ! input
real*8 p(nx,ny,nz)    ! output
real*8 work(nx,ny,nz)   
integer i,j,k,setbdy


if (ndim==2) then
   k=1

#ifdef COMPACT
   do j=ny1,ny2
   do i=nx1,nx2

      work(i,j,k)= &
            (                                       &
               (f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k)) + &
               (f(i,j+1,k)-2*f(i,j,k)+f(i,j-1,k))   &
             ) / 12

   enddo
   enddo
   do j=ny1,ny2
   do i=nx1,nx2
      f(i,j,k)=f(i,j,k)+work(i,j,k)
   enddo
   enddo
#endif

      
   if (setbdy==1) then
   ! now overwrite boundary with same data as in p
   if (my_x==0) then
      do j=ny1,ny2
         f(nx1,j,k)=p(nx1,j,k)
      enddo
   endif
   if (my_x==ncpu_x-1) then
      do j=ny1,ny2
         f(nx2,j,k)=p(nx2,j,k)
      enddo
   endif
   if (my_y==0) then
      do i=nx1,nx2
         f(i,ny1,k)=p(i,ny1,k)
      enddo
   endif
   if (my_y==ncpu_y-1) then
      do i=nx1,nx2
         f(i,ny2,k)=p(i,ny2,k)
      enddo
   endif
   endif

else
   stop 'helm_rhs_correction: 3D not yet coded'
endif


end subroutine




subroutine zero_boundary(f)
!
!  input/output: f
!  f set to zero on the boundary
!
use params
implicit none
real*8 f(nx,ny,nz)    ! input
!local
integer n,i,j,k

   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      if ( (my_x==0 .and. i==nx1)  .or. &
           (my_x==ncpu_x-1 .and. i==nx2)  .or. &
           (my_y==0 .and. j==ny1) .or. &
           (my_y==ncpu_y-1 .and. j==ny2) .or.&
           (my_z==0 .and. k==nz1 .and. ndim==3) .or. &
           (my_z==ncpu_z-1 .and. k==nz2 .and. ndim==3) ) then
         f(i,j,k)=0
      endif
   enddo
   enddo
   enddo

end subroutine



subroutine helmholtz_dirichlet(f,lf,alpha,beta,work)
!
!  input: f
!  output: lf
!
!     lf = [alpha + beta*laplacian](f)
!
!     lf on the boundary copied from f.
!
!  2nd order uses the regular stencil:     
!  x=1/delx**2
!  y=1/dely**2
!                 y
!           x  -2x-2y   x
!                 y
!
!  4th order compact uses above + a  * D2Y D2
!  where a= (1/x+1/y)/12
!
!
!                 y                   xy  -2xy    xy
!           x  -2x-2y   x   +   a   -2xy   4xy  -2xy 
!                 y                   xy  -2xy    xy
!
! which is:
!
!              axy       y  -2axy         axy
!           x-2axy    -2x-2y+4axy      x-2axy
!              axy       y  -2axy         axy
!
use params
use ghost
implicit none
real*8 f(nx,ny,nz)    ! input
real*8 lf(nx,ny,nz)    ! output
real*8 work(nx,ny,nz)
real*8 :: alpha
real*8 :: beta

!local
integer n,i,j,k


if (ndim==2) then
   call ghost_update_x(f,1)

   k=1

#ifdef COMPACT
   ! work = Dx(f)
   do j=ny1,ny2
   do i=nx1,nx2
      work(i,j,k)=(f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k))/(delx*delx)
   enddo
   enddo
   call ghost_update_y(work,1)
#endif

   call ghost_update_y(f,1)

   do j=ny1,ny2
   do i=nx1,nx2
      if ( (my_x==0 .and. i==nx1)  .or. &
           (my_x==ncpu_x-1 .and. i==nx2)  .or. &
           (my_y==0 .and. j==ny1) .or. &
           (my_y==ncpu_y-1 .and. j==ny2) ) then
         lf(i,j,k)=f(i,j,k)
         cycle
      endif

      lf(i,j,k)=(f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k))/(delx*delx) + &
                (f(i,j+1,k)-2*f(i,j,k)+f(i,j-1,k))/(dely*dely)


      ! add the DyDx(f) term:
#ifdef COMPACT
      lf(i,j,k)=lf(i,j,k) + ((delx*delx+dely*dely)/12)* &
             (work(i,j+1,k)-2*work(i,j,k)+work(i,j-1,k))/(dely*dely)
#endif

      if (alpha==0) then
         lf(i,j,k)=beta*lf(i,j,k)
      else
         lf(i,j,k)=alpha*f(i,j,k)+beta*lf(i,j,k) 
#ifdef COMPACT
         lf(i,j,k)=lf(i,j,k) + &
            alpha*(                                       &
                     (f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k)) + &
                     (f(i,j+1,k)-2*f(i,j,k)+f(i,j-1,k))   &
                   ) / 12
#endif

      endif
      
   enddo
   enddo
else if (ndim==3) then
   call ghost_update_x(f,1)
   call ghost_update_y(f,1)
   call ghost_update_z(f,1)
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      if ( (my_x==0 .and. i==nx1)  .or. &
           (my_x==ncpu_x-1 .and. i==nx2)  .or. &
           (my_y==0 .and. j==ny1) .or. &
           (my_y==ncpu_y-1 .and. j==ny2) .or.&
           (my_z==0 .and. k==nz1) .or. &
           (my_z==ncpu_z-1 .and. k==nz2) ) then
         lf(i,j,k)=f(i,j,k)
         cycle
      endif

      
      lf(i,j,k)=alpha*f(i,j,k)+ &
            beta*(f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k))/(delx*delx) +&
            beta*(f(i,j+1,k)-2*f(i,j,k)+f(i,j-1,k))/(dely*dely) +&
            beta*(f(i,j,k+1)-2*f(i,j,k)+f(i,j,k-1))/(delz*delz) 
   enddo
   enddo
   enddo

endif


end subroutine


subroutine helmholtz_periodic_ghost(f,lf,alpha,beta,work)
!
!  no boundary conditions:   useing ghost cell data
!               works for PERIODIC and RELFECT, REFLECT_ODD
!  input: f
!  output: lf
!
!     lf = [alpha + beta*laplacian](f)
!
!  2nd order uses the regular stencil:     
!  x=1/delx**2
!  y=1/dely**2
!                 y
!           x  -2x-2y   x
!                 y
!
!  4th order compact uses above + a  * D2Y D2
!  where a= (1/x+1/y)/12
!
!
!                 y                   xy  -2xy    xy
!           x  -2x-2y   x   +   a   -2xy   4xy  -2xy 
!                 y                   xy  -2xy    xy
!
! which is:
!
!              axy       y  -2axy         axy
!           x-2axy    -2x-2y+4axy      x-2axy
!              axy       y  -2axy         axy
!
use params
use ghost
implicit none
real*8 f(nx,ny,nz)    ! input
real*8 lf(nx,ny,nz)    ! output
real*8 work(nx,ny,nz)
real*8 :: alpha
real*8 :: beta

!local
integer n,i,j,k


if (ndim==2) then
   call ghost_update_x(f,1)

   k=1

#ifdef COMPACT
   ! work = Dx(f)
   do j=ny1,ny2
   do i=nx1,nx2
      work(i,j,k)=(f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k))/(delx*delx)
   enddo
   enddo
   call ghost_update_y(work,1)
#endif

   call ghost_update_y(f,1)

   do j=ny1,ny2
   do i=nx1,nx2
      lf(i,j,k)=(f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k))/(delx*delx) + &
                (f(i,j+1,k)-2*f(i,j,k)+f(i,j-1,k))/(dely*dely)

#if 0
      if (i==3 .and. j==130) then
         print *,i,j,lf(i,j,k)
         print *,'x: ',f(i-1,j,k),f(i,j,k),f(i+1,j,k)
         print *,'y: ',f(i,j-1,k),f(i,j,k),f(i,j+1,k)
         stop
      endif
#endif

      ! add the DyDx(f) term:
#ifdef COMPACT
      lf(i,j,k)=lf(i,j,k) + ((delx*delx+dely*dely)/12)* &
             (work(i,j+1,k)-2*work(i,j,k)+work(i,j-1,k))/(dely*dely)
#endif

      if (alpha==0) then
         lf(i,j,k)=beta*lf(i,j,k)
      else
         lf(i,j,k)=alpha*f(i,j,k)+beta*lf(i,j,k) 
#ifdef COMPACT
         lf(i,j,k)=lf(i,j,k) + &
            alpha*(                                       &
                     (f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k)) + &
                     (f(i,j+1,k)-2*f(i,j,k)+f(i,j-1,k))   &
                   ) / 12
#endif

      endif
      
   enddo
   enddo
else if (ndim==3) then
   call ghost_update_x(f,1)
   call ghost_update_y(f,1)
   call ghost_update_z(f,1)
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      lf(i,j,k)=alpha*f(i,j,k)+ &
            beta*(f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k))/(delx*delx) +&
            beta*(f(i,j+1,k)-2*f(i,j,k)+f(i,j-1,k))/(dely*dely) +&
            beta*(f(i,j,k+1)-2*f(i,j,k)+f(i,j,k-1))/(delz*delz) 
   enddo
   enddo
   enddo

endif


end subroutine






subroutine helmholtz_hform_periodic(f,lf,alpha,beta,h)
!
!  input: f
!  output: lf
!
!      lf = [h*alpha + beta*div h grad](f)
!
use params
use fft_interface
use transpose
implicit none
real*8 f(nx,ny,nz)    ! input
real*8 lf(nx,ny,nz)    ! output
real*8 h(nx,ny,nz)    ! height field (used only for shallow water equations)
real*8 :: alpha
real*8 :: beta

!local
real*8 work(nx,ny,nz)    ! work array
real*8 gradf(nx,ny,nz,2) ! work array
real*8 fxx(nx,ny,nz)     ! work array
real*8 dummy
integer n

do n=1,2
   call der(f,gradf(1,1,1,n),dummy,work,DX_ONLY,n)
   gradf(:,:,:,n)=h*gradf(:,:,:,n)
enddo
call divergence(lf,gradf,fxx,work)
lf=alpha*f*h + beta*lf

end subroutine








subroutine helmholtz_dirichlet_inv(f,work,alpha,beta)
!
!  solve [alpha + beta*laplacian](p) = f
!  input:  f 
!  ouput:  f   will be overwritten with the solution p
!  b.c. for f specified in f. 
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
real*8 phi(nx,ny,nz)    
real*8 lphi(nx,ny,nz)    
real*8 xm,ym,zm,xfac
real*8 :: axy,x_axy,y_axy
real*8 :: xb(ny,2),yb(nx,2)

#ifdef COMPACT
if (ndim/=2) then
   call abort("helmholtz_dirichlet_inv() not coded for 3D compact")
endif
#endif

if (beta==0) then
   f=f/alpha
   return
endif
!
!  let phi = f on the boundary, 0 inside
!      lphi = Helmholtz(phi)
!
!  solve:  Helmholtz(psi) = b - lphi  with psi=0 on boundary
!          then set psi=psi+phi
!
!
!NOTE: we dont need phi, lphi.  when this code is debugged, replace by:
! save boundary data in single 1D array
! set f = 0 on boundary
! using boundary data alone: compute f = f - lphi 
!      along points just inside boundary
!
! solve as before
! restor boundary conditions in f from 1D array.

! copy boundary data into temp array:
! replace (interior only) f with f - lphi
if (my_x==0) then
   i=nx1
   k=1
   do j=ny1,ny2
      xb(j,1)=f(i,j,k)
   enddo
endif
if (my_x==ncpu_x-1) then
   i=nx2
   k=1
   do j=ny1,ny2
      xb(j,2)=f(i,j,k)
   enddo
endif
if (my_y==0) then
   j=ny1
   k=1
   do i=nx1,nx2
      yb(i,1)=f(i,j,k)
   enddo
endif
if (my_y==ncpu_y-1) then
   j=ny2
   k=1
   do i=nx1,nx2
      yb(i,2)=f(i,j,k)
   enddo
endif

! correct corner points:
if (my_x==0 .and. my_y==0) then
endif
if (my_x==ncpu_x-1 .and. my_y==0) then
endif
if (my_x==0 .and. my_y==ncpu_y-1) then
endif
if (my_x==ncpu_x-1 .and. my_y==ncpu_y-1) then
endif




! phi = f on boundary, zero inside
phi=f
call zero_boundary(phi)
phi=f-phi
call helmholtz_dirichlet(phi,lphi,alpha,beta,work)


! convert problem to zero b.c. problem
f=f-lphi
call zero_boundary(f)
! solve Helm(f) = b with 0 b.c.


call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d) 
call sinfft1(work,n1,n1d,n2,n2d,n3,n3d)     
call transpose_from_x(work,f,n1,n1d,n2,n2d,n3,n3d) 


call transpose_to_y(f,work,n1,n1d,n2,n2d,n3,n3d)  ! x,y,z -> y,x,z
call sinfft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(work,f,n1,n1d,n2,n2d,n3,n3d) 

call transpose_to_z(f,work,n1,n1d,n2,n2d,n3,n3d)  
call sinfft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(work,f,n1,n1d,n2,n2d,n3,n3d)  


!  2nd order uses the regular stencil:     
!  x=1/delx**2
!  y=1/dely**2
!                 y
!           x  -2x-2y   x
!                 y
!                                                      
!  4th order compact uses above + a  * D2Y D2X         
!  where a= (1/x+1/y)/12                               
!
!
!                 y                   xy  -2xy    xy
!           x  -2x-2y   x   +   a   -2xy   4xy  -2xy 
!                 y                   xy  -2xy    xy
!
! which is:
!
!              axy       y  -2axy         axy
!           x-2axy    -2x-2y+4axy      x-2axy
!              axy       y  -2axy         axy
!
axy=(delx**2+dely**2)/(12*delx*delx*dely*dely)
x_axy=1/(delx*delx) - 2*axy
y_axy=1/(dely*dely) - 2*axy

do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
#ifdef COMPACT
   ! mode imsine(i),jmsine(j),kmsine(k)
   ! sin(x+h)sin(y+h) + sin(x-h)sin(y+h) + 
   ! sin(x+h)sin(y-h) + sin(x-h)sin(y-h) =
   ! 
   !   [  2*cos(h)  ]  sin(x) [sin(y+h)+sin(y-h) ]  = 
   !
   !   [  4*cos(hx) cos(hy) ]  sin(x) sin(y)
   
   xfac = x_axy*2*cos(pi*delx*imsine(i))         ! x term
   xfac = xfac + y_axy*2*cos(pi*dely*jmsine(j))  ! y term
   ! diagonal terms:
   xfac = xfac + axy*4*cos(pi*delx*imsine(i))*cos(pi*dely*jmsine(j))
   ! center term:
   xfac = xfac -2/(delx*delx)-2/(dely*dely)+4*axy   
#else
   ! [sin(x+h) - 2 sin(x) + sin(x-h)]   = 
   !
   ! [ -2 + 2*cos(h)  ]  sin(x)
   ! 
   xm = -2 + 2*cos(pi*delx*imsine(i))
   xm=xm/(delx*delx)
   ym = -2 + 2*cos(pi*dely*jmsine(j))
   ym=ym/(dely*dely)
   zm = -2 + 2*cos(pi*delz*kmsine(k))
   zm=zm/(delz*delz)
   xfac=xm+ym+zm
#endif



!   if (abs(f(i,j,k))>1e-6) then
!      print *,imsine(i),jmsine(j),kmsine(k),f(i,j,k)
!   endif
   xfac = alpha + beta*xfac
   if (imsine(i)+jmsine(j)+kmsine(k)==0) then
      f(i,j,k) = 0
   else
      f(i,j,k) = f(i,j,k)/xfac
   endif
     
 
enddo
enddo
enddo

call transpose_to_z(f,work,n1,n1d,n2,n2d,n3,n3d)       
call isinfft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(work,f,n1,n1d,n2,n2d,n3,n3d)       

call transpose_to_y(f,work,n1,n1d,n2,n2d,n3,n3d)       
call isinfft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(work,f,n1,n1d,n2,n2d,n3,n3d)         ! y,x,z -> x,y,z

call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d) 
call isinfft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_x(work,f,n1,n1d,n2,n2d,n3,n3d )


!print *,maxval(abs(f(nx1:nx2,ny1:ny2,1)-lphi(nx1:nx2,ny1:ny2,1)))


!restor boundary values
if (my_x==0) then
   i=nx1
   k=1
   do j=ny1,ny2
      f(i,j,k)=xb(j,1)
   enddo
endif
if (my_x==ncpu_x-1) then
   i=nx2
   k=1
   do j=ny1,ny2
      f(i,j,k)=xb(j,2)
   enddo
endif
if (my_y==0) then
   j=ny1
   k=1
   do i=nx1,nx2
      f(i,j,k)=yb(i,1)
   enddo
endif
if (my_y==ncpu_y-1) then
   j=ny2
   k=1
   do i=nx1,nx2
      f(i,j,k)=yb(i,2)
   enddo
endif

!f=f+phi
end subroutine












subroutine helmholtz_periodic_inv(f,work,alpha,beta)
!
!  solve [alpha + beta*laplacian](p) = f
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

if (beta==0) then
   f=f/alpha
   return
endif


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
































!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute 3D sin transform along first dimension
!
! input data:   (for example, if n1=5)
!   1 2 3 4 5      x(1)=x(5)=0
! odd extension:
!   1 2 3 4 5 -4 -3 -2  0 0      (fft does not include last periodic point)
!                                and we add 2 for padding
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sinfft1(p,n1,n1d,n2,n2d,n3,n3d)
use fft_interface
implicit none
integer n1,n1d,n2,n2d,n3,n3d
real*8 :: p(n1d,n2d,n3d)
real*8 :: w(2*n1,n2)

integer i,j,k
if (n1==1) return

do k=1,n3
   do j=1,n2
      ! make an odd extension
      do i=1,n1-1
         w(i,j)=p(i,j,k)               !  1 2 3 4 
         w(i+n1-1,j)=-p(n1-i+1,j,k)    !          -5 -4 -2
      enddo
   enddo
   call fft1(w,2*n1-2,2*n1,n2,n2,1,1)
   do j=1,n2
      ! save only the sine modes:
      do i=1,n1-1
         p(i,j,k)=w(2*i,j)
      enddo
   enddo
enddo
end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute 3D sin transform along first dimension
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine isinfft1(p,n1,n1d,n2,n2d,n3,n3d)
use fft_interface
implicit none
integer n1,n1d,n2,n2d,n3,n3d
real*8 :: p(n1d,n2d,n3d)
real*8 :: w(2*n1,n2)

integer i,j,k
if (n1==1) return


do k=1,n3
   do j=1,n2
      ! unpack the sine modes into a sine/cosine array:
      ! save only the sine modes:
      do i=1,n1-1
         w(2*i-1,j)=0
         w(2*i,j)=p(i,j,k)
      enddo
   enddo
   call ifft1(w,2*n1-2,2*n1,n2,n2,1,1)
   do j=1,n2
      do i=1,n1
         p(i,j,k)=w(i,j)
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

