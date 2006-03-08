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

if (g_bdy_x1/=PERIODIC ) then
   call abort('der() can only handle periodic boundaries')
endif
if (g_bdy_y1/=PERIODIC) then
   call abort('der() can only handle periodic boundaries')
endif
if (g_bdy_z1/=PERIODIC ) then
   call abort('der() can only handle periodic boundaries')
endif

if (numerical_method==FOURTH_ORDER) then

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
real*8 :: one=1
integer i,j,k,n
external helmholtz_periodic,helmholtz_dirichlet


! solve laplacian(p)=div(u)

if ( g_bdy_x1==PERIODIC .and. &
     g_bdy_y1==PERIODIC .and. &
     g_bdy_z1==PERIODIC) then
   call divergence(p,u,work,work2)
   call helmholtz_periodic_inv(p,work,alpha,beta,one)

   !work=p  ! RHS
   !p=0  ! initial guess
   !tol=1e-10
   !call cgsolver(p,work,alpha,beta,tol,work2,helmholtz_periodic,.false.)
else
   ! might try divfree_ghost
   stop 'divfree_gridspace: only supports periodic/reflection case'
endif




! compute u=u-grad(p)
do n=1,3
   call der(p,work,dummy,work2,1,n)
   if (n==3) work=work/Lz
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      u(i,j,k,n)=u(i,j,k,n)-work(i,j,k)
   enddo
   enddo
   enddo
enddo

if (dealias>0) then
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

if (n_var<3) call abort("vorticity() requires n_var>2")

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
      if (n==2) vor(i,j,k,1) = vor(i,j,k,1) -d1(i,j,k)/Lz
      if (n==1) vor(i,j,k,2) = vor(i,j,k,2) +d1(i,j,k)/Lz
   enddo
   enddo
   enddo
   endif

enddo


end subroutine


subroutine potential_vorticity(pv,vor,u,d1,work,pv_type)
use params
use fft_interface
use transpose
implicit none
!
! For the moment we've stored theta in u(:,:,:,4)
!
!

real*8 u(nx,ny,nz,n_var)    ! input
real*8 vor(nx,ny,nz,3)      ! input
real*8 d1(nx,ny,nz)     !  work array used for derivatives
real*8 work(nx,ny,nz)   ! work array 
real*8 pv(nx,ny,nz)          ! output
integer pv_type
!bw  pv_type
!bw  full    = 1   q = omega_a dot grad theta + fcor d theta/dz -bous omega_3
!bw  i         2   q = omega dot grad theta
!bw  ii        3   q = omega_3 - fcor/bous d theta / dz
!bw  iii       4   q = d theta / d z
!bw  iv        5   q = omega_3
!bw  v         6   q = bous omega_3 + omega_i diff_i theta
!bw  vi        7   q = fcor d theta / d z + omega_i diff_i theta

! local variables
integer i,j,k,n
real*8 dummy(1)

if (n_var<3) call abort("potential vorticity() requires n_var>2")

call vorticity(Q,vor,work1,work2)
   
if(pv_type == 1) then
   !bw
   !bw Next, dot this into the gradiant of the density.
   !bw
   !bw Assume the density, (\tilde(rho)) is in u(:,:,:,4).
   !bw
   !compute theta_x
   call der(u(1,1,1,4),d1,dummy,work,DX_ONLY,1)
   pv = d1*vor(:,:,:,1)
   
   ! compute theta_y
   call der(u(1,1,1,4),d1,dummy,work,DX_ONLY,2)
   pv=pv + d1*vor(:,:,:,2)
   
   ! compute theta_z -- do not forget to add in coriolis here
   call der(u(1,1,1,4),d1,dummy,work,DX_ONLY,3)
   pv = pv + (d1/Lz)*(vor(:,:,:,3)+fcor)
   pv = pv - bous*vor(:,:,:,3)
   
elseif (pv_type == 2) then
!
! pv = omega dot grad theta
!
   !compute theta_x
   call der(u(1,1,1,4),d1,dummy,work,DX_ONLY,1)
   pv = d1*vor(:,:,:,1)
   call der(u(1,1,1,4),d1,dummy,work,DX_ONLY,2)
   pv=pv + d1*vor(:,:,:,2)

   ! compute theta_z -- do not forget to add in coriolis here
   call der(u(1,1,1,4),d1,dummy,work,DX_ONLY,3)
   pv = pv + (d1/Lz)*vor(:,:,:,3)
   
elseif (pv_type == 3) then
   !bw
   !bw pv = -bous * omega_3 +  fcor theta / d z
   !bw
   ! compute theta_z -- do not forget to add in coriolis here
   call der(u(1,1,1,4),d1,dummy,work,DX_ONLY,3)
   pv = - bous*vor(:,:,:,3)+fcor*d1/Lz

elseif (pv_type == 4) then
   !bw
   !bw pv = fcor * d theta / d z
   !bw
   pv = 0
   ! compute theta_z -- do not forget to add in coriolis here
   call der(u(1,1,1,4),d1,dummy,work,DX_ONLY,3)
   pv = fcor*d1/Lz
elseif (pv_type == 5) then
   !bw
   !bw pv = -bous * omega_3
   !bw
   pv = -bous*vor(:,:,:,3)

elseif (pv_type == 6) then
!
! pv = omega dot grad theta - bous * omega_3
!

   !compute theta_x
   call der(u(1,1,1,4),d1,dummy,work,DX_ONLY,1)
   pv = d1*vor(:,:,:,1)
   
   ! compute theta_y
   call der(u(1,1,1,4),d1,dummy,work,DX_ONLY,2)
   pv = pv + d1*vor(:,:,:,2)

   ! compute theta_z -- do not forget to add in coriolis here
   call der(u(1,1,1,4),d1,dummy,work,DX_ONLY,3)
   pv = pv + (d1/Lz)*vor(:,:,:,3)
   pv = pv - bous*vor(:,:,:,3)

elseif (pv_type == 7) then
!
! pv = omega dot grad theta + fcor  d theta / d z
!
   !compute theta_x
   call der(u(1,1,1,4),d1,dummy,work,DX_ONLY,1)
   pv = d1*vor(:,:,:,1)

   ! compute theta_y
   call der(u(1,1,1,4),d1,dummy,work,DX_ONLY,2)
   pv = pv + d1*vor(:,:,:,2)
   
   ! compute theta_z -- do not forget to add in coriolis here
   call der(u(1,1,1,4),d1,dummy,work,DX_ONLY,3)
   pv = pv + (d1/Lz)*(vor(:,:,:,3)+fcor)

endif


end subroutine



subroutine vorticity2d(vor,u,d1,work)
use params
use fft_interface
use transpose
implicit none
real*8 u(nx,ny,nz,2)    ! input
real*8 vor(nx,ny,nz)    ! output
real*8 d1(nx,ny,nz) 
real*8 work(nx,ny,nz) 

! local variables
integer i,j,k,n
real*8 dummy(1)

vor=0

   n=2
   call der(u(1,1,1,n),d1,dummy,work,DX_ONLY,1)
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      vor(i,j,k) = vor(i,j,k) + d1(i,j,k)
   enddo
   enddo
   enddo

   n=1
   ! compute u_y, u_yy
   call der(u(1,1,1,n),d1,dummy,work,DX_ONLY,2)
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      vor(i,j,k) = vor(i,j,k) -d1(i,j,k)
   enddo
   enddo
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
   div(i,j,k) = div(i,j,k)+work1(i,j,k)/Lz
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

implicit none
real*8 f(nx,ny,nz)    ! input
real*8 fout(g_nz2,nslabx,ny_2dz)  ! output
real*8 work(nx,ny,nz) ! work array1
integer n1,n1d,n2,n2d,n3,n3d


!
!  for the full spectral method, we will be working mostly in
!  the z-transform fourier space, 
!     pt(g_nz2,nslabx,ny_2dz)
!  so we also need that ny_2dz is even
if (mod(ny_2dz,2)/=0) then
   call abort("ny_2dz is not even.  cant use z-decomp FFTs")
endif




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

if (numerical_method==FOURTH_ORDER) then
if (Lz/=1) call abort("dzscale must be 1 for FD methods")
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
         xfac= alpha + beta*(-im*im -km*km/(Lz*Lz) - jm*jm)*pi2_squared      
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
!  but is usefull to simulate 4th order code in the spectral code
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
real*8 gradf(nx,ny,nz) ! work array
real*8 fxx(nx,ny,nz)     ! work array
real*8 dummy
integer n

! linear helmholtz operator
lf=alpha*f
if (beta==0) return
do n=1,3
   call der(f,gradf,fxx,work,DX_AND_DXX,n)
   if (n==3) fxx=fxx/Lz/Lz
   lf=lf+beta*fxx
enddo

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
!   gradf(:,:,:,n)=H0*gradf(:,:,:,n)
   gradf(:,:,:,n)=h*gradf(:,:,:,n)
enddo
call divergence(lf,gradf,fxx,work)
lf=alpha*f*h + beta*lf

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
real*8 :: dzscale


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
! input data:   (for example, if n1=5)  g_nx=o_nx=5
!   1 2 3 4 5      x(1)=x(5)=0
! odd extension:
!   1 2 3 4 5 -4 -3 -2  0 0      (fft does not include last periodic point)
!                                and we add 2 for padding
! fft:
!   0 -0 1  -1  2  -2  3  -3 4 -4
! repacked fft:
!   0 4 1  -1  2  -2  3  -3  
! removal of cosine terms:
!  4  -1  -2 -3              (we picked up a cos(4x) mode by mistake, instead
!                             of the sin(0x) mode.  But sin(0x) mode is always 0
!                             and cos(4x) mode is zero since our data is odd)
!
! if offset_bdy, then the trailing 0 is missing:
!  g_nx=5  o_nx=6    fft: 2*g_nx
!
!  input:            1 2 3 4 5     x(1)=x(6)=0  (x(6) not included)
!  odd extension:    1 2 3 4 5 0 -5 -4 -3 -2  0 0 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sinfft1(p,n1,n1d,n2,n2d,n3,n3d)
use params
use fft_interface
implicit none
integer n1,n1d,n2,n2d,n3,n3d
real*8 :: p(n1d,n2d,n3d)
real*8 :: w(2*n1+2,n2)

integer i,j,k
if (n1==1) return

if (1==offset_bdy) then
do k=1,n3
   do j=1,n2
      ! make an odd extension
      w(n1,j)=p(n1,j,k)
      w(n1+1,j)=0
      do i=1,n1-1
         w(i,j)=p(i,j,k)               !  1 2 3 4 5
         w(i+n1+1,j)=-p(n1-i+1,j,k)    !            * -5 -4 -3 -2
      enddo
   enddo
   call fft1(w,2*n1,2*n1+2,n2,n2,1,1)
   do j=1,n2
      ! save only the sine modes:
      do i=1,n1
         p(i,j,k)=w(2*i,j)
      enddo
   enddo
enddo
else
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
endif
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute 3D sin transform along first dimension
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine isinfft1(p,n1,n1d,n2,n2d,n3,n3d)
use fft_interface
use params
implicit none
integer n1,n1d,n2,n2d,n3,n3d
real*8 :: p(n1d,n2d,n3d)
real*8 :: w(2*n1+2,n2)

integer i,j,k
if (n1==1) return

if (1==offset_bdy) then
do k=1,n3
   do j=1,n2
      ! unpack the sine modes into a sine/cosine array:
      ! save only the sine modes:
      do i=1,n1
         w(2*i-1,j)=0
         w(2*i,j)=p(i,j,k)
      enddo
   enddo
   call ifft1(w,2*n1,2*n1+2,n2,n2,1,1)
   do j=1,n2
      do i=1,n1
         p(i,j,k)=w(i,j)
      enddo
   enddo
enddo
else
do k=1,n3
   do j=1,n2
      ! unpack the sine modes into a sine/cosine array:
      ! save only the sine modes:
      do i=1,n1-1
         w(2*i-1,j)=0
         w(2*i,j)=p(i,j,k)
      enddo
   enddo
   call ifft1(w,2*n1-2,2*n1+2,n2,n2,1,1)
   do j=1,n2
      do i=1,n1
         p(i,j,k)=w(i,j)
      enddo
   enddo
enddo
endif
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

            if (dealias_remove(im,jm,km)) then 
               p(i,j,k)=0
            endif
         enddo
      enddo
   enddo

end subroutine






subroutine global_min(p,mn)
use params
use mpi
implicit none
real*8 :: p(nx,ny,nz)
real*8 :: mn,mn2
integer :: ierr

mn=minval(p(nx1:nx2,ny1:ny2,nz1:nz2))
#ifdef USE_MPI
mn2=mn
call mpi_allreduce(mn2,mn,1,MPI_REAL8,MPI_MIN,comm_3d,ierr)
#endif

end subroutine




subroutine global_max(p,mx)
use params
use mpi
implicit none
real*8 :: p(nx,ny,nz)
real*8 :: mx,mx2
integer :: ierr

mx=maxval(p(nx1:nx2,ny1:ny2,nz1:nz2))
#ifdef USE_MPI
mx2=mx
call mpi_allreduce(mx2,mx,1,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif

end subroutine


subroutine ke_shell_z(Qhat,hscale,nvar2)
!
! compute the KE in shell   kstart2 < k**2 <= kstop2
! for z-decomposition data
!
!  Note: for Lz<>1, this routine assumes scaling domain by Lz 
!  
!
use params
use mpi
implicit none
integer :: kstart2,kstop2,nvar2
real*8 Qhat(g_nz2,nslabx,ny_2dz,nvar2)           ! Fourier data at time t
real*8 :: ke(3),ke2(3),xw,u2,xfac,ierr,hscale(3),xw2,cfl
integer :: im,jm,km,i,j,k,n,km_start,jm_start,im_start

! units of E = m**2/s**2
! units of E(k) = m**3/s**2
! hyper viscosity:  E(kmax)* k**8 / kmax**(8-1.5)  
! scaling:  E(kmax)/(kmax**2)(4-.75)

! du/dt = (sqrt(E(kmax))  [1/ (kmax**8 kmax**alpha) ] k**8  u  
! m/s**2  =  m**1.5/s  kmax**-alpha  m/s
!  1 = m**1.5 kmax**-alpha     = m**(1.5+alpha)   alpha = -1.5

if (dealias==1 .or. dealias==0) then
   ! dealias_remove = ( (km>g_nz/3)  .or.  (jm>g_ny/3)  .or. (im>g_nx/3) )
   ! take energy in band km such that:  km+1>g_nz/3 .and. km< g_nz/3
   km_start = g_nz/3
   if (km_start >= g_nz/3) km_start=km_start-1
   jm_start = g_ny/3
   if (jm_start >= g_ny/3.0) jm_start=jm_start-1
   im_start = g_nx/3.0
   if (im_start >= g_nx/3) im_start=im_start-1

else if (dealias==2) then
   kstart2=dealias_sphere_kmax2_1
   kstop2=dealias_sphere_kmax2
else
   call abort('ke_shell_z():  Error: bad dealiasing type')
endif



   ke=0
   
   do j=1,ny_2dz
      jm=z_jmcord(j)
      do i=1,nslabx
         im=z_imcord(i)
         do k=1,g_nz
            km=z_kmcord(k)
            xfac = 2*2*2
            if (km==0) xfac=xfac/2
            if (jm==0) xfac=xfac/2
            if (im==0) xfac=xfac/2

            u2=0
            do n=1,nvar2
               u2=u2+Qhat(k,i,j,n)*Qhat(k,i,j,n)
            enddo

            if (dealias==1) then
               if (abs(im)==im_start) ke(1)=ke(1) + .5*xfac*u2
               if (abs(jm)==jm_start) ke(2)=ke(2) + .5*xfac*u2
               if (abs(km)==km_start) ke(3)=ke(3) + .5*xfac*u2
            endif
            if (dealias==2) then
               xw = jm*jm + im*im + km*km 
               if (kstart2 < xw  .and. xw <= kstop2) then
                  ke(1) = ke(1) + .5*xfac*u2
               endif
            endif
         enddo
      enddo
   enddo
#ifdef USE_MPI
   ke2 = ke
   call mpi_allreduce(ke2,ke,3,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif
if (dealias==1) then
   hscale(1) = sqrt(ke(1)) * (pi2_squared*im_start**2)**(-(mu_hyper-.75))
   hscale(2) = sqrt(ke(2)) * (pi2_squared*jm_start**2)**(-(mu_hyper-.75))
   hscale(3) = sqrt(ke(3)) * (pi2_squared*(km_start/Lz)**2)**(-(mu_hyper-.75))

   im=g_nx/3
   jm=g_ny/3
   km=g_nz/3

!
!  make sure that hyper viscosity does not exceed a CFL condition 
!  d(u_k)/dt = x u_k      delt*x < cfl   x < cfl/delt
!
!   physical units:   du/dt = mu k'**mu_hyper u
!   code units:       du/dt  = mu (k/Lz) **mu_hyper u
!
!   mu (kmax/Lz)**mu_hyper < cfl/delt
!
   cfl= 1      
   if (delt>0) then
      cfl = cfl/(3*delt)  ! divide by 3 and apply speretaly to each term:
   endif

   xw2=mu_hyper_value*hscale(1)*(im*im*pi2_squared)**mu_hyper
   max_hyper(1)= xw2
   if (xw2>cfl) then
      print *,'x: warning: hyper viscosity CFL: ',xw2*delt
      hscale(1)=(cfl/xw2)*hscale(1)
      xw2=mu_hyper_value*hscale(1)*(im*im*pi2_squared)**mu_hyper
      max_hyper(1)= xw2*delt
   endif
   
   xw2=mu_hyper_value*hscale(2)*(jm*jm*pi2_squared)**mu_hyper
   max_hyper(2)= xw2
   if (xw2>cfl) then
      print *,'y: warning: hyper viscosity CFL: ',xw2*delt
      hscale(2)=(cfl/xw2)*hscale(2)
      xw2=mu_hyper_value*hscale(2)*(jm*jm*pi2_squared)**mu_hyper
      max_hyper(2)= xw2*delt
   endif
   
   xw2=mu_hyper_value*hscale(3)*(km*km*pi2_squared/(Lz*Lz))**mu_hyper
   max_hyper(3)= xw2
   if (xw2>cfl) then
      print *,'z: warning: hyper viscosity CFL: ',xw2*delt
      hscale(3)=(cfl/xw2)*hscale(3)
      xw2=mu_hyper_value*hscale(3)*(km*km*pi2_squared/(Lz*Lz))**mu_hyper
      max_hyper(3)= xw2*delt
   endif
endif


if (dealias==2) then
   hscale(1) = sqrt(ke(1)) * (pi2_squared*kstop2)**(-(mu_hyper-.75))
   hscale(2)=hscale(1)
   hscale(3)=hscale(1)
endif

end subroutine




subroutine ke_shell(Qhat,ke,numk,kstart2,kstop2)
!
! compute the KE in shell   kstart2 < k**2 <= kstop2
! for z-decomposition data
!
use params
use mpi
implicit none
real*8 Qhat(nx,ny,nz,n_var)        ! Fourier data at time t
integer :: numk,kstart2,kstop2
real*8 :: ke,xw,u2,xfac,ierr
integer :: im,jm,km,i,j,k,n


   ke=0
   numk=0
   
   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)
            xw = jm*jm + im*im + km*km 
            if (kstart2 < xw  .and. xw <= kstop2) then
               numk=numk+1
               xfac = 2*2*2
               if (km==0) xfac=xfac/2
               if (jm==0) xfac=xfac/2
               if (im==0) xfac=xfac/2
               
               u2=0
               do n=1,ndim
                  u2=u2+Qhat(i,j,k,n)*Qhat(i,j,k,n)
               enddo
               ke = ke + .5*xfac*u2
            endif
         enddo
      enddo
   enddo
#ifdef USE_MPI
   xw = ke
   call mpi_allreduce(xw,ke,1,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif
end subroutine




