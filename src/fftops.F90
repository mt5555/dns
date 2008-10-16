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
   call abortdns('der() can only handle periodic boundaries')
endif
if (g_bdy_y1/=PERIODIC) then
   call abortdns('der() can only handle periodic boundaries')
endif
if (g_bdy_z1/=PERIODIC ) then
   call abortdns('der() can only handle periodic boundaries')
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
   call helmholtz_periodic_inv(p,work,alpha,beta)

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

if (n_var<3) call abortdns("vorticity() requires n_var>2")

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

real*8 u(nx,ny,nz,n_var)    ! input state variable
real*8 vor(nx,ny,nz,3)      ! vorticty (computed and returned output)
real*8 d1(nx,ny,nz)         ! work array used for derivatives
real*8 work(nx,ny,nz)       ! work array 
real*8 pv(nx,ny,nz)         ! output
integer pv_type
!bw  pv_type
!bw  full    = 1   q = omega dot grad theta + fcor d theta/dz -bous omega_3
!bw  i         2   q = omega dot grad theta    (3D NSE with p-scalar)
!bw  ii        3   q =  -bous * omega_3 +  fcor theta / d (QG)
!bw  iii       4   q = fcor * d theta / d z	(Ro->0, Fr=1)
!bw  iv        5   q = -bous omega_3		(Fr->0, Ro=1)
!bw  v         6   q = bous omega_3 + omega_i diff_i theta 
!bw  vi        7   q = fcor d theta / d z + omega_i diff_i theta

! local variables
integer i,j,k,n
real*8 dummy(1)


if (n_var<3) call abortdns("potential vorticity() requires n_var>2")
!if(pv_type/=1) call abortdns("fftops: potential_vorticity")

call vorticity(vor,u,d1,work)

   
if(pv_type == 1) then
   !bw
   !bw Next, dot this into the gradiant of the density.
   !bw
   !bw Assume the density, (\tilde(rho)) is in u(:,:,:,4).
   !bw
   !compute theta_x
   call der(u(1,1,1,np1),d1,dummy,work,DX_ONLY,1)
   pv = d1*vor(:,:,:,1)
   
   ! compute theta_y
   call der(u(1,1,1,np1),d1,dummy,work,DX_ONLY,2)
   pv=pv + d1*vor(:,:,:,2) 
   ! compute theta_z -- do not forget to add in coriolis here
   call der(u(1,1,1,np1),d1,dummy,work,DX_ONLY,3)
   pv = pv + (d1/Lz)*(vor(:,:,:,3)+fcor)
   pv = pv - bous*vor(:,:,:,3)
   
elseif (pv_type == 2) then
!
! pv = omega dot grad theta
!
   !compute theta_x
   call der(u(1,1,1,np1),d1,dummy,work,DX_ONLY,1)
   pv = d1*vor(:,:,:,1)

   !compute theta_y
   call der(u(1,1,1,np1),d1,dummy,work,DX_ONLY,2)
   pv=pv + d1*vor(:,:,:,2)

   ! compute theta_z 
   call der(u(1,1,1,np1),d1,dummy,work,DX_ONLY,3)
   pv = pv + (d1/Lz)*vor(:,:,:,3)
   
elseif (pv_type == 3) then
   !bw
   !bw pv = -bous * omega_3 +  fcor theta / d z
   !bw
   ! compute theta_z -- do not forget to add in coriolis here
   call der(u(1,1,1,np1),d1,dummy,work,DX_ONLY,3)
   pv = - bous*vor(:,:,:,3)+fcor*d1/Lz

elseif (pv_type == 4) then
   !bw
   !bw pv = fcor * d theta / d z
   !bw
   pv = 0
   ! compute theta_z -- do not forget to add in coriolis here
   call der(u(1,1,1,np1),d1,dummy,work,DX_ONLY,3)
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
   call der(u(1,1,1,np1),d1,dummy,work,DX_ONLY,1)
   pv = d1*vor(:,:,:,1)
   
   ! compute theta_y
   call der(u(1,1,1,np1),d1,dummy,work,DX_ONLY,2)
   pv = pv + d1*vor(:,:,:,2)

   ! compute theta_z 
   call der(u(1,1,1,np1),d1,dummy,work,DX_ONLY,3)
   pv = pv + (d1/Lz)*vor(:,:,:,3)
   pv = pv - bous*vor(:,:,:,3)

elseif (pv_type == 7) then
!
! pv = omega dot grad theta + fcor  d theta / d z
!
   !compute theta_x
   call der(u(1,1,1,np1),d1,dummy,work,DX_ONLY,1)
   pv = d1*vor(:,:,:,1)

   ! compute theta_y
   call der(u(1,1,1,np1),d1,dummy,work,DX_ONLY,2)
   pv = pv + d1*vor(:,:,:,2)
   
   ! compute theta_z -- do not forget to add in coriolis here
   call der(u(1,1,1,np1),d1,dummy,work,DX_ONLY,3)
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


subroutine v_vorticity2d(vor,work)
!
!  compute   v_vorticity = (1-alpha^2 laplacian) u_vorticity
!
!  vor should contain the u vorticity.  It will be overwritten
!  with the v-vorticity
!
! note: highest mode tweaked so that laplacian = div grad
!
use params
use fft_interface
use transpose
implicit none
real*8 vor(nx,ny,nz)    ! input and output
real*8 work(nx,ny,nz) 

! local variables
integer i,j,k,im,jm,km

if (alpha_value>0) then
call fft3d(vor,work)
do k=nz1,nz2
   km=kmcord(k)
   if (km==g_nz/2) km=0
   do j=ny1,ny2
      jm=jmcord(j)
      if (jm==g_ny/2) jm=0
      do i=nx1,nx2
         im=imcord(i)
         if (im==g_nx/2) im=0
             if(infinite_alpha ==1) then
                vor(i,j,k) = vor(i,j,k)*((im*im +km*km/(Lz*Lz) + jm*jm)*pi2_squared)         
             else       
                 vor(i,j,k) = vor(i,j,k)*(1 + (alpha_value**2)*(im*im +km*km/(Lz*Lz) + jm*jm)*pi2_squared)
         endif
      enddo
   enddo
enddo
call ifft3d(vor,work)
endif

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
!  f stored using reference (nx,ny,nz) decompostion
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
!  f uses reference (nx,ny,nz) decompostion
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
!  f:     data stored with reference (nx,ny,nz) decomposition
!  fout:  data stored with z-decompostion 
!
!  f,fout can overlap in memory
!  data in f is ovewritten
!
use params
use fft_interface
use transpose

implicit none
real*8 f(nx,ny,nz)    ! input
real*8 fout(g_nz2,nx_2dz,ny_2dz)  ! output
real*8 work(nx,ny,nz) ! work array1
integer n1,n1d,n2,n2d,n3,n3d


!
!  for the full spectral method, we will be working mostly in
!  the z-transform fourier space, 
!     pt(g_nz2,nx_2dz,ny_2dz)
!  so we also need that ny_2dz is even
if (mod(ny_2dz,2)/=0) then
   call abortdns("ny_2dz is not even.  cant use z-decomp FFTs")
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





subroutine zx_fft3d_trashinput(f,fout,work)
!
!  compute fft of f, return in fout.
!  f:     data stored with x-pencil decompostion
!  fout:  data stored with z-pencil decompostion 
!
!  f,fout can overlap in memory
!  data in f is ovewritten
!
use params
use fft_interface
use transpose

implicit none
real*8 f(nx,ny,nz)    ! input
real*8 fout(g_nz2,nx_2dz,ny_2dz)  ! output
real*8 work(nx,ny,nz) ! work array1
integer n1,n1d,n2,n2d,n3,n3d


!
!  for the full spectral method, we will be working mostly in
!  the z-transform fourier space, 
!     pt(g_nz2,nx_2dz,ny_2dz)
!  so we also need that ny_2dz is even
if (mod(ny_2dz,2)/=0) then
   call abortdns("ny_2dz is not even.  cant use z-decomp FFTs")
endif


n1=g_nx
n1d=g_nx2
n2=nslabz
n2d=nslabz
n3=ny_2dx
n3d=ny_2dx
call fft1(f,n1,n1d,n2,n2d,n3,n3d)     
call transpose_from_x(f,work,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_y(work,f,n1,n1d,n2,n2d,n3,n3d)
call fft1(f,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(f,work,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_z(work,fout,n1,n1d,n2,n2d,n3,n3d)
call fft1(fout,n1,n1d,n2,n2d,n3,n3d)


end




subroutine z_ifft3d(fin,f,work)
!
!  compute inverse fft 3d of fin, return in f
!  fin:  data stored with z-pencil decompostion 
!  f:    data stored with reference (nx,ny,nz) decomposition
!
!  fin and f can overlap in memory
!
use params
use fft_interface
use transpose
implicit none
! true size of all arrays:  nx,ny,nz:
! work used as a work array and its shape must match fin
real*8 fin(g_nz2,nx_2dz,ny_2dz)  ! input
real*8 f(nx,ny,nz)               ! output
real*8 work(g_nz2,nx_2dz,ny_2dz) ! work array


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k

n1=g_nz
n1d=g_nz2   	
n2=nx_2dz
n2d=nx_2dz
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






subroutine zx_ifft3d(fin,f,work)
!
!  compute inverse fft 3d of fin, return in f
!  fin:  data stored with z-pencil decompostion 
!  f:    data stored with x-pencil decompostion
!
!  fin and f can overlap in memory
!
use params
use fft_interface
use transpose
implicit none
! true size of all arrays:  nx,ny,nz:
! f used as a work array and its shape must match fin
real*8 fin(g_nz2,nx_2dz,ny_2dz)  ! input
real*8 work(nx,ny,nz)    
real*8 f(g_nz2,nx_2dz,ny_2dz) 


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k

n1=g_nz
n1d=g_nz2   	
n2=nx_2dz
n2d=nx_2dz
n3=ny_2dz
n3d=ny_2dz

f=fin
call ifft1(f,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(f,work,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_y(work,f,n1,n1d,n2,n2d,n3,n3d)
call ifft1(f,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(f,work,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_x(work,f,n1,n1d,n2,n2d,n3,n3d)
call ifft1(f,n1,n1d,n2,n2d,n3,n3d)

end





subroutine zx_ifft3d_and_dx(fin,f,fx,work)
!
!  compute inverse fft 3d of fin, return in f
!  also return df/dx in fx                               
!
!  fin:   data stored with z-pencil decompostion 
!  f,fx:  data stored with x-pencil decompostion
!
!  fin and f can overlap in memory
!
use params
use fft_interface
use transpose
implicit none
! true size of all arrays:  nx,ny,nz:
! f used as a work array and its shape must match fin
real*8 fin(g_nz2,nx_2dz,ny_2dz)  ! input
real*8 f(g_nz2,nx_2dz,ny_2dz)    ! output
real*8 fx(nx,ny,nz)              ! output
real*8 work(nx,ny,nz)            ! work array


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k

n1=g_nz
n1d=g_nz2   	
n2=nx_2dz
n2d=nx_2dz
n3=ny_2dz
n3d=ny_2dz

f=fin   
call ifft1(f,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(f,work,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_y(work,f,n1,n1d,n2,n2d,n3,n3d)
call ifft1(f,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(f,work,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_x(work,f,n1,n1d,n2,n2d,n3,n3d)

call mult_by_ik(f,fx,n1,n1d,n2,n2d,n3,n3d)       ! compute fx code
call ifft1(fx,n1,n1d,n2,n2d,n3,n3d)                 ! compute fx code

call ifft1(f,n1,n1d,n2,n2d,n3,n3d)




end




subroutine zx_ifft3d_and_dy(fin,f,fy,work)
!
!  compute inverse fft 3d of fin, return in f
!  also return df/dy in fy                               
!
!  fin:   data stored with z-pencil decompostion 
!  f,fy:  data stored with x-pencil decompostion
!
!  fin and f can overlap in memory
!
use params
use fft_interface
use transpose
implicit none
! true size of all arrays:  nx,ny,nz:
! f used as a work array and its shape must match fin
real*8 fin(g_nz2,nx_2dz,ny_2dz)  ! input
real*8 f(g_nz2,nx_2dz,ny_2dz)  ! output
real*8 fy(nx,ny,nz)              ! output
real*8 work(nx,ny,nz)



!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k

n1=g_nz
n1d=g_nz2   	
n2=nx_2dz
n2d=nx_2dz
n3=ny_2dz
n3d=ny_2dz

f=fin
call ifft1(f,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(f,work,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_y(work,f,n1,n1d,n2,n2d,n3,n3d)

! f -> df/dy   
call mult_by_ik(f,fy,n1,n1d,n2,n2d,n3,n3d)        ! new code for fy
call ifft1(fy,n1,n1d,n2,n2d,n3,n3d)                  ! new code for fy
call transpose_from_y(fy,work,n1,n1d,n2,n2d,n3,n3d)    ! new code for fy
call transpose_to_x(work,fy,n1,n1d,n2,n2d,n3,n3d)      ! new code for fy
call ifft1(fy,n1,n1d,n2,n2d,n3,n3d)                  ! new code for fy

call ifft1(f,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(f,work,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_x(work,f,n1,n1d,n2,n2d,n3,n3d)
call ifft1(f,n1,n1d,n2,n2d,n3,n3d)



end







subroutine z_ifft3d_and_dx(fin,f,fx,work)
!
!  compute inverse fft 3d of fin, return in f
!  also return df/dx in fx                               
!
!  fin:   data stored with z-pencil decompostion 
!  f,fx:  data stored with reference (nx,ny,nz) decompostion
!
!  fin and f can overlap in memory
!
use params
use fft_interface
use transpose
implicit none
! true size of all arrays:  nx,ny,nz:
! work used as a work array and its shape must match fin
real*8 fin(g_nz2,nx_2dz,ny_2dz)  ! input
real*8 f(nx,ny,nz)               ! output
real*8 fx(nx,ny,nz)              ! output
real*8 work(g_nz2,nx_2dz,ny_2dz) ! work array


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k

n1=g_nz
n1d=g_nz2   	
n2=nx_2dz
n2d=nx_2dz
n3=ny_2dz
n3d=ny_2dz

work=fin
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(work,f,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_y(f,work,n1,n1d,n2,n2d,n3,n3d)
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(work,f,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d)

call mult_by_ik(work,f,n1,n1d,n2,n2d,n3,n3d)       ! compute fx code
call ifft1(f,n1,n1d,n2,n2d,n3,n3d)                 ! compute fx code
call transpose_from_x(f,fx,n1,n1d,n2,n2d,n3,n3d)   ! compute fx code

call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_x(work,f,n1,n1d,n2,n2d,n3,n3d)



end




subroutine z_ifft3d_and_dy(fin,f,fy,work)
!
!  compute inverse fft 3d of fin, return in f
!  also return df/dy in fy                               
!
!  fin:   data stored with z-pencil decompostion 
!  f,fy:  data stored with reference (nx,ny,nz) decompostion
!
!  fin and f can overlap in memory
!
use params
use fft_interface
use transpose
implicit none
! true size of all arrays:  nx,ny,nz:
! work used as a work array and its shape must match fin
real*8 fin(g_nz2,nx_2dz,ny_2dz)  ! input
real*8 f(nx,ny,nz)  ! output
real*8 fy(nx,ny,nz) ! output
real*8 work(g_nz2,nx_2dz,ny_2dz) ! work array


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k

n1=g_nz
n1d=g_nz2   	
n2=nx_2dz
n2d=nx_2dz
n3=ny_2dz
n3d=ny_2dz

work=fin
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(work,f,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_y(f,work,n1,n1d,n2,n2d,n3,n3d)

! work -> df/dy    (using f as a work array)
call mult_by_ik(work,f,n1,n1d,n2,n2d,n3,n3d)        ! new code for fy
call ifft1(f,n1,n1d,n2,n2d,n3,n3d)                  ! new code for fy
call transpose_from_y(f,fy,n1,n1d,n2,n2d,n3,n3d)    ! new code for fy
call transpose_to_x(fy,f,n1,n1d,n2,n2d,n3,n3d)      ! new code for fy
call ifft1(f,n1,n1d,n2,n2d,n3,n3d)                  ! new code for fy
call transpose_from_x(f,fy,n1,n1d,n2,n2d,n3,n3d)    ! new code for fy

! work -> f
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(work,f,n1,n1d,n2,n2d,n3,n3d)

call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d)
call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_x(work,f,n1,n1d,n2,n2d,n3,n3d)


end








subroutine mult_by_ik(p,px,n1,n1d,n2,n2d,n3,n3d)
!
!  take derivative in Fourier space
!
use params
implicit none
integer n1,n1d,n2,n2d,n3,n3d
real*8 p(n1d,n2d,n3d)
real*8 px(n1d,n2d,n3d)
integer i,j,k,m


do k=1,n3
   do j=1,n2
      ! note: for i=2, m=0, we are actually working with the cos(n1/2) mode
      ! but d/dx of this mode goes to sin(n1/2) = 0, so just take m=0 
      do i = 1, n1,2
         m=(i-1)/2
         px(i+1,j,k) = pi2* m * p(i,j,k)
         px(i,j,k) = -pi2 *m * p(i+1,j,k)
      enddo
   enddo
enddo
end subroutine mult_by_ik




subroutine z_fft3d_nvar(Q_grid,Q,work1,work2) 
! convinience function for routines which only know Qn
! with dimensions nx,ny,nz
!
! input:  Q_grid:   state vector in grid space (nx,ny,nz) decomposition
! output: Q         state vector in spectral space, z-pencil decomposition
use params
implicit none
real*8 :: Q_grid(nx,ny,nz,n_var)
real*8 :: Q(g_nz2,nx_2dz,ny_2dz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
integer n
do n=1,n_var
   work2=Q_grid(:,:,:,n)
   call z_fft3d_trashinput(work2,Q(1,1,1,n),work1)
enddo
end subroutine



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
if (Lz/=1) call abortdns("dzscale must be 1 for FD methods")
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Filter out spherical wavenumbers greater than some spec_max 
! (low-pass filter)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_filter_trunc(p)
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
            if ((im**2 + jm**2 + km**2 ) > spec_max**2) then 
               p(i,j,k)=0
            endif
         enddo
      enddo
   enddo

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Filter out spherical wavenumbers less than some spec_max 
!(high-pass filter) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_filter_hpass(p) 
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
            if ((im**2 + jm**2 + km**2 ) < spec_max**2) then                    
               p(i,j,k)=0                                                       
            endif                                                               
         enddo                                                                  
      enddo                                                                     
   enddo                                                                        
                                                                                
end subroutine        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Filter out spherical wavenumbers other than k_shell specified 
! (band-filter)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_filter_shell00(p,kshell)
use params
implicit none
real*8 p(nx,ny,nz)
integer kshell

integer i,j,k,im,jm,km,ks
real*8 xw2

   do k=nz1,nz2
      km=(kmcord(k))
      do j=ny1,ny2
         jm=(jmcord(j))
         do i=nx1,nx2
            im=(imcord(i))
            xw2 =  (im**2 + jm**2 + (km/Lz)**2 )
            ks = nint(sqrt(xw2))
            if (ks /= kshell) p(i,j,k)=0
         enddo
      enddo
   enddo

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Filter out spherical wavenumbers other than k_shell specified 
! (band-filter)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_filter_shell(p,kshell)
use params
implicit none
real*8 p(nx,ny,nz)
integer kshell

integer i,j,k,im,jm,km,ks
real*8 xw2,xw,xl,xr

   do k=nz1,nz2
      km=(kmcord(k))
      do j=ny1,ny2
         jm=(jmcord(j))
         do i=nx1,nx2
            im=(imcord(i))
            xw2 =  (im**2 + jm**2 + (km/Lz)**2 )
            ks = nint(sqrt(xw2))
!            if (ks /= kshell) p(i,j,k)=0
		xw=sqrt(xw2)
		xl=kshell-0.5
		xr=kshell+0.5
    if( xw.le.xl.or.xw.ge.xr ) p(i,j,k)=0

         enddo
      enddo
   enddo

end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Filter out spherical wavenumbers other than k_shell specified 
! (band-filter)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_filter_shell1(p,kshell)
use params
implicit none
real*8 p(nx,ny,nz)
integer kshell

integer i,j,k,im,jm,km,ks
real*8 xw2,xw,xl,xr

   do k=nz1,nz2
      km=(kmcord(k))
      do j=ny1,ny2
         jm=(jmcord(j))
         do i=nx1,nx2
            im=(imcord(i))
            xw2 =  (im**2 + jm**2 + (km/Lz)**2 )
            ks = nint(sqrt(xw2))
!            if (ks /= kshell) p(i,j,k)=0
		xw=sqrt(xw2)
		xl=kshell-0.01
		xr=kshell+0.01
    if( xw.le.xl.or.xw.ge.xr ) p(i,j,k)=0

         enddo
      enddo
   enddo

end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Filter out spherical wavenumbers other than k_shell specified 
! (band-filter)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_filter_shell0(p,kshell)
use params
implicit none
real*8 p(nx,ny,nz)
integer kshell

integer i,j,k,im,jm,km,ks
real*8 xw2

   do k=nz1,nz2
      km=(kmcord(k))
      do j=ny1,ny2
         jm=(jmcord(j))
         do i=nx1,nx2
            im=(imcord(i))
            if ( im**2+jm**2+km**2 .ge. kshell) then
               ! keep just the cosine mode (k,0,0)
            else
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

subroutine global_max_abs(p,mx)
use params
use mpi
implicit none
real*8 :: p(nx,ny,nz)
real*8 :: mx,mx2
integer :: ierr

mx=maxval(abs(p(nx1:nx2,ny1:ny2,nz1:nz2)))
#ifdef USE_MPI
mx2=mx
call mpi_allreduce(mx2,mx,1,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif

end subroutine




subroutine ke_shell_z(Qhat,hscale)
!
!  Note: for Lz<>1, this routine assumes scaling domain by Lz 
!
! compute a hyperviscosity scaling based on the KE in the last shell
! for spherical dealiasing, we follow Leslie Smith (ref?):
!
!    du/dt + ... =  (h k^2)^mu_hyper  u
!
!    KE = energy in the last spherical shell.
!    scalar hypervis coefficient chosen of the form: 
!        h = [ sqrt(KE)*  kmax^2P   ]   ^(1/mu_hyper) 
!
!
! Spherical shells my not exist when running large aspect ration cases,
! and do not make sense when using 2/3 dealiasing.  So we introduce
! a tensor viscosity h(1:3):
!
!    du/dt + ... =  div ( h dot grad u) 
!                = (h(1) kx^2 + h(2) ky^2 + h(3) kz^2 )^(mu_hyper) u
!    
!    KE(:) = KE in last slab kx_max,,ky_max and kz_max
!    scalar hypervis coefficient chosen of the form: 
!        h(1) = [ sqrt(KE(1)) kx_max^2P ]^(1/mu_hyper) 
!        h(2) = [ sqrt(KE(2)) ky_max^2P ]^(1/mu_hyper) 
!        h(3) = [ sqrt(KE(3)) kz_max^2P ]^(1/mu_hyper) 
!
! P chosen to be dimensionally correct:
!
! units of E = m**2/s**2
! units of E(k) = m**3/s**2
!
!    du/dt   =    sqrt(KE)   kx_max^2P  k^(2 mu_hyper)  u
!    m/s**2  =    m**1.5/s   m^(-2P)    m^(-2 mu_hyper)  m/s          
!      1     =    m**1.5     m^(-2P)    m^(-2 mu_hyper)
!      1.5 -2P - 2mu_hyper = 0   
!
!      P  = .75 - mu_hyper
!    
! 
! What happens if we want to scale by ENS = k^2 KE ?
!        h = [ ENS^? *  kmax^2P    ]   ^(1/mu_hyper) 
!        h = [ KE^?  *  kmax^(2P+2)]   ^(1/mu_hyper) 
!  for correct time units, ? must be 0.5, and thus this produces the same results.
!
!We also compute scaling for the scalars:
!  hscale(:,1)     viscosity in (x,y,z) directions used for momentum
!  hscale(:,np1:np2)  viscosity in (x,y,z) directions used for scalars
!
use params
use mpi
implicit none
integer :: kstart2,kstop2
real*8 Qhat(g_nz2,nx_2dz,ny_2dz,n_var)           ! Fourier data at time t
real*8 hscale(n_var,n_var)
real*8 :: ke(3),ke2(3),xw,u2,xfac,ierr,xw2,cfl
integer :: im,jm,km,i,j,k,n,km_start,jm_start,im_start,shell_type

!
!  shell_type = 1     use spherical shell of thickness 1 wave number   
!               2     use slabs of thickness 1 wave number
!
if (dealias==1 .or. dealias==0) then
   ! use SLABS
   shell_type=2 
   ! dealias_remove = ( (km>g_nz/3)  .or.  (jm>g_ny/3)  .or. (im>g_nx/3) )
   ! take energy in band km such that:  km+1>g_nz/3 .and. km< g_nz/3
   km_start = g_nz/3
   jm_start = g_ny/3
   im_start = g_nx/3

! uncomment next 4 lines use Energy scaling based on spherical shell
! contained inside 2/3 dealiased cube:
    shell_type=1
    kstart2=(g_nx/3 - 1 )**2
    kstop2=(g_nx/3)**2
!   print *,'using hyper viscosity energy scaling based on shell: ',sqrt(real(kstart2)),sqrt(real(kstop2))

else if (dealias==2) then
   ! use spherical shell
   shell_type=1   
   kstart2=dealias_sphere_kmax2_1
   kstop2=dealias_sphere_kmax2
else
   call abortdns('ke_shell_z():  Error: bad dealiasing type')
endif


!
!  compute KE in last shell, or last 3 slabs in wave number space:
!
ke=0
do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nx_2dz
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)
         xfac = 2*2*2
         if (km==0) xfac=xfac/2
         if (jm==0) xfac=xfac/2
         if (im==0) xfac=xfac/2
         
         u2=0
         do n=1,ndim
            u2=u2+Qhat(k,i,j,n)*Qhat(k,i,j,n)
         enddo
         
         if (shell_type==2) then
            if (abs(im)==im_start) ke(1)=ke(1) + .5*xfac*u2
            if (abs(jm)==jm_start) ke(2)=ke(2) + .5*xfac*u2
            if (abs(km)==km_start) ke(3)=ke(3) + .5*xfac*u2
         endif
         if (shell_type==1) then
            xw = jm*jm + im*im + km*km/(Lz*Lz) 
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



!
!  Compute hyperviscosity coefficient, and check for viscous CFL violations:
!
if (shell_type==1) then
   hscale(1,1) = sqrt(ke(1)) * (pi2_squared*kstop2)**(-(mu_hyper-.75))
   hscale(1,1) = hscale(1,1)**(1d0/mu_hyper)
   hscale(2,1)=hscale(1,1)
   hscale(3,1)=hscale(1,1)
endif
if (shell_type==2) then
   hscale(1,1) = sqrt(ke(1)) * (pi2_squared*im_start**2)**(-(mu_hyper-.75))
   hscale(2,1) = sqrt(ke(2)) * (pi2_squared*jm_start**2)**(-(mu_hyper-.75))
   if (km_start==0) then ! 2D problems
      hscale(3,1)=0
   else
      hscale(3,1) = sqrt(ke(3)) * (pi2_squared*(km_start/Lz)**2)**(-(mu_hyper-.75))
   endif
   hscale(1:3,1) = hscale(1:3,1)**(1d0/mu_hyper)
endif

!
! check for CFL violations
!
if (dealias==0) then  ! no dealiasing
   im=g_nx/2
   jm=g_ny/2
   km=g_nz/2
   cfl = .9
endif
if (dealias==1) then ! 2/3 rule
   im=im_start
   jm=jm_start
   km=km_start
   cfl = 2.0   ! hyper8:  2.9 blows up  2.8 seems stable
endif
if (dealias==2) then ! spherical
   ! im^2 + jm^2 + km^2 =  dealias_sphere_kmax2    
   im= sqrt(dealias_sphere_kmax2/3.0)    
   jm=im
   km=im
   cfl = .9   ! hyper8:  1.5 blows up, 1.25 seems stable, lets use .9
endif


!
!  make sure that hyper viscosity does not exceed a CFL condition 
!  d(u_k)/dt = x u_k      delt*x < cfl   x < cfl/delt
!
!   physical units:   du/dt = mu k'**mu_hyper u
!   code units:       du/dt  = mu (k/Lz) **mu_hyper u
!
!   mu (kmax/Lz)**mu_hyper < cfl/delt
!
xw2=hscale(1,1)*(im*im*pi2_squared)
xw2=xw2+hscale(2,1)*(jm*jm*pi2_squared)
xw2=xw2+hscale(3,1)*(km*km*pi2_squared/(Lz*Lz))
xw2=mu_hyper_value*xw2**mu_hyper


! only enforce viscous CFL for explict viscosity:
if (delt*xw2>cfl .and. hyper_implicit==0) then
   !      if (my_pe==io_pe) then
   !         print *,'delt = ', delt
   !         print *,'mu_hyper_value = ', mu_hyper_value
   !         print *,'hscale = ', hscale(1,1), hscale(2,1), hscale(3,1)
   !         print *,'im,jm,km = ',im,jm,km
   !         print *,'mu_hyper = ',mu_hyper
!   print *,'warning: velocity hyper viscosity CFL: ',delt*xw2
   !     endif
   !  scale 'hscale' by alpha so that 
   !     (alpha hscale wavenumber_stuff)**mu_hyper = cfl/delt
   !     alpha**mu_hyper xw2 = cfl/delt
   !     alpha**mu_hyper = cfl/(delt*xw2)
   hscale(1:3,1)= hscale(1:3,1)*( cfl/(delt*xw2) )**(1d0/mu_hyper)
   
   ! recompute:
   xw2=hscale(1,1)*(im*im*pi2_squared)
   xw2=xw2+hscale(2,1)*(jm*jm*pi2_squared)
   xw2=xw2+hscale(3,1)*(km*km*pi2_squared/(Lz*Lz))
   xw2=mu_hyper_value*xw2**mu_hyper
endif
max_hyper = delt*xw2







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  now repeat all of that, but this time for passive scalars
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do n=np1,np2
   ke=0
   do j=1,ny_2dz
      jm=z_jmcord(j)
      do i=1,nx_2dz
         im=z_imcord(i)
         do k=1,g_nz
            km=z_kmcord(k)
            xfac = 2*2*2
            if (km==0) xfac=xfac/2
            if (jm==0) xfac=xfac/2
            if (im==0) xfac=xfac/2
            
            u2=Qhat(k,i,j,n)*Qhat(k,i,j,n)
            
            if (shell_type==2) then
               if (abs(im)==im_start) ke(1)=ke(1) + .5*xfac*u2
               if (abs(jm)==jm_start) ke(2)=ke(2) + .5*xfac*u2
               if (abs(km)==km_start) ke(3)=ke(3) + .5*xfac*u2
            endif
            if (shell_type==1) then
               xw = jm*jm + im*im + km*km/(Lz*Lz)
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
   
   if (shell_type==1) then
      hscale(1,n) = sqrt(ke(1)) * (pi2_squared*kstop2)**(-(mu_hyper-.75))
      hscale(1,n) = hscale(1,n)**(1d0/mu_hyper)
      hscale(2,n)=hscale(1,n)
      hscale(3,n)=hscale(1,n)
   endif
   if (shell_type==2) then
      hscale(1,n) = sqrt(ke(1)) * (pi2_squared*im_start**2)**(-(mu_hyper-.75))
      hscale(2,n) = sqrt(ke(2)) * (pi2_squared*jm_start**2)**(-(mu_hyper-.75))
      if (km_start==0) then ! 2D problems
         hscale(3,n)=0
      else
         hscale(3,n) = sqrt(ke(3)) * (pi2_squared*(km_start/Lz)**2)**(-(mu_hyper-.75))
      endif
      hscale(1:3,n) = hscale(1:3,n)**(1d0/mu_hyper)
   endif

   ! now check CFL:
   if (dealias==0) then  ! no dealiasing
      im=g_nx/2
      jm=g_ny/2
      km=g_nz/2
   endif
   if (dealias==1) then ! 2/3 rule
      im=im_start
      jm=jm_start
      km=km_start
   endif
   if (dealias==2) then ! spherical
      ! im^2 + jm^2 + km^2 =  dealias_sphere_kmax2    
      im= sqrt(dealias_sphere_kmax2/3.0)    
      jm=im
      km=im
   endif
   
   xw2=hscale(1,n)*(im*im*pi2_squared)
   xw2=xw2+hscale(2,n)*(jm*jm*pi2_squared)
   xw2=xw2+hscale(3,n)*(km*km*pi2_squared/(Lz*Lz))
   xw2=mu_hyper_value*xw2**mu_hyper

   ! enforce CFL only for explicit hyper viscosity
   if (delt*xw2>cfl .and. hyper_implicit==0) then   
      !        if(my_pe==io_pe) then
      !            print *,'warning: scalar hyper viscosity CFL: ',delt*xw2
      !         endif
      !  scale 'hscale' by alpha so that 
      !     (alpha hscale wavenumber_stuff)**mu_hyper = cfl/delt
      !     alpha**mu_hyper xw2 = cfl/delt
      !     alpha**mu_hyper = cfl/(delt*xw2)
      hscale(1:3,n)=hscale(1:3,n)*( cfl/(delt*xw2) )**(1d0/mu_hyper)
   endif
enddo


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



subroutine hyper_filter(Qhat,dt)
use params
implicit none
! 
! apply hyper viscosity as a filter to Q
!
! INPUT:
!  Q = Q(t+1):   should be a valid state after an inviscid RK stage:   
!        Q(t+1) = Q(t) + dt*RHS  
!
! OUTPUT 
!    Q(t+1) =  Q(t+1) / (1+dt*k^8)  
!
! This is equivilent to an implicit viscosity term since:
!
!           Q(t+1) = Q(t) + dt*RHS  - dt*k^8 Q(t+1)
!   (1+dt*k^8) Q(t+1) =  Q(t) + dt*RHS
!    Q(t+1) = [ Q(t) + dt*RHS ] / (1+dt*k^8)
!
!
real*8 Qhat(g_nz2,nx_2dz,ny_2dz,n_var)           ! Fourier data at time t
real*8 dt

! local variables:
real*8 hyper_scale(n_var,n_var)
real*8 xw2
integer i,j,k,im,jm,km,n,xfac
real*8 ke0(n_var)
real*8 ke1(n_var)

if (mu_hyper<2) return
if (hyper_implicit /= 1 ) return

! compute hyper viscosity scaling based on energy in last shell:
! print *,'calling ke_shell  Q=',(qhat(1,1,1,1:3))
call ke_shell_z(Qhat,hyper_scale)


ke0=0
ke1=0

do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nx_2dz
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

         xfac = 2*2*2
         if (km==0) xfac=xfac/2
         if (jm==0) xfac=xfac/2
         if (im==0) xfac=xfac/2
         


         xw2=hyper_scale(1,1)*(im*im*pi2_squared)
         xw2=xw2+hyper_scale(2,1)*(jm*jm*pi2_squared)
         xw2=xw2+hyper_scale(3,1)*(km*km*pi2_squared/(Lz*Lz))
         do n=1,ndim
            ke0(n)=ke0(n)+.5*xfac*Qhat(k,i,j,n)**2
            Qhat(k,i,j,n)  = Qhat(k,i,j,n) / ( 1+dt*mu_hyper_value*xw2**mu_hyper )
            ke1(n)=ke1(n)+.5*xfac*Qhat(k,i,j,n)**2
         enddo

         do n=np1,np2
            xw2=hyper_scale(1,n)*(im*im*pi2_squared)
            xw2=xw2+hyper_scale(2,n)*(jm*jm*pi2_squared)
            xw2=xw2+hyper_scale(3,n)*(km*km*pi2_squared/(Lz*Lz))
            ke0(n)=ke0(n)+.5*xfac*Qhat(k,i,j,n)**2
            Qhat(k,i,j,n)  = Qhat(k,i,j,n) / ( 1+dt*mu_hyper_value*xw2**mu_hyper )
            ke1(n)=ke1(n)+.5*xfac*Qhat(k,i,j,n)**2
         enddo

      enddo
   enddo
enddo
do n=1,n_var
!   print *,n,ke0(n),ke1(n)
enddo



end subroutine






subroutine sincos_to_complex_field(p,cmodes_r,cmodes_i)
#if 0
  Convert set of sine and cosine FFT coefficients to complex coefficients.

  Note: the N/2 mode (for which the standard FFT allows only the cosine
  component) will be ignored (assumbed to be zero) by this routine.
  This mode is usually zero if any type of dealiasing is used.

  After calling this routine, here is a loop that will extract the modes:

do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   ! wave number (im,jm,km)   im positive or negative
   im=imcord_exp(i)
   jm=jmcord_exp(j)
   km=kmcord_exp(k)
   real_part = c_r(i,j,k)
   imag_part = c_i(i,j,k)
enddo
enddo
enddo

For details of the 8x8 transform that maps sin/cos to complex
coefficients, see the sincos_to_complex routine in sforcing.F90
#endif      
use params
implicit none
integer :: nmax
real*8 :: p(nx,ny,nz)
real*8 :: cmodes_r(nx,ny,nz)
real*8 :: cmodes_i(nx,ny,nz)
real*8 :: a,b,mx
integer :: i,j,k,im,jm,km,ii,jj,kk,sm,i0,i1,j0,j1,k0,k1,ck



cmodes_r=0
cmodes_i=0

do k=nz1,nz2,2
do j=ny1,ny2,2
do i=nx1,nx2,2

   i0=i
   i1=i+1
   j0=j
   j1=j+1
   k0=k
   k1=k+1
   
   ! verify that this processer has all 8 modes:
   if (  abs(imcord_exp(i0))/=abs(imcord_exp(i1)) .or. &
         abs(jmcord_exp(j0))/=abs(jmcord_exp(j1)) .or. & 
         abs(kmcord_exp(k0))/=abs(kmcord_exp(k1)) ) then
      print *,'not all modes are on this processor: '
      print *,i0,imcord_exp(i0),imcord_exp(i1)
      print *,j0,jmcord_exp(j0),jmcord_exp(j1)
      print *,k0,kmcord_exp(k0),kmcord_exp(k1)
      call abortdns(" ")
   endif

   ! make sure that i0 is the positive mode, i1 is the negative mode:
   if ( imcord_exp(i0)<0 .or. imcord_exp(i1)>0 .or.  &
        jmcord_exp(j0)<0 .or. jmcord_exp(j1)>0 .or.  & 
        kmcord_exp(k0)<0 .or. kmcord_exp(k1)>0 ) then
      print *,'we have the sign wrong?'
      print *,i0,imcord_exp(i0),imcord_exp(i1)
      print *,j0,imcord_exp(j0),imcord_exp(j1)
      print *,k0,imcord_exp(k0),imcord_exp(k1)
      call abortdns(" ")
   endif


   ! loop over the 8 sin/cos modes, mapping into 8 complex modes
   do ii=i0,i1
   do jj=j0,j1
   do kk=k0,k1      

      ! sin/cos mode 
      ! these arrays are the same as imcord(), except the last sin(N/2) wave
      ! number is zero instead of N/2.  This code sets all N/2 modes to zero
      im=imcord_exp(ii)  
      jm=jmcord_exp(jj)
      km=kmcord_exp(kk)

      a=0; b=0
      ! count the number if sin() terms:
      sm=0; if (im<0) sm=sm+1;  if (jm<0) sm=sm+1;  if (km<0) sm=sm+1
      if (sm==0) then
         a=p(ii,jj,kk)/8
      else if (sm==1) then
         b=-p(ii,jj,kk)/8
      else if (sm==2) then
         a=-p(ii,jj,kk)/8
      else if (sm==3) then
         b=p(ii,jj,kk)/8
      else
         call abortdns("this cant happen")
      endif

      ! cos(N/2) mode is stored with the cos(0) mode.
      ! we ignore cos(N/2) mode, and pretend it is the sin(0) mode
      ! which must be zero:
      if ( imcord(ii)==g_nx/2 .or. jmcord(jj)==g_ny/2 .or. &
          kmcord(kk)==g_nz/2 ) then
         a=0
         b=0
      endif



      cmodes_r(i0,j0,k0)=cmodes_r(i0,j0,k0) + a;    
      cmodes_i(i0,j0,k0)=cmodes_i(i0,j0,k0) + b

      cmodes_r(i0,j0,k1)=cmodes_r(i0,j0,k1) + a*sign(1,km)   
      cmodes_i(i0,j0,k1)=cmodes_i(i0,j0,k1) + b*sign(1,km)
      
      cmodes_r(i0,j1,k0)=cmodes_r(i0,j1,k0) + a*sign(1,jm)
      cmodes_i(i0,j1,k0)=cmodes_i(i0,j1,k0) + b*sign(1,jm)
      
      cmodes_r(i0,j1,k1)=cmodes_r(i0,j1,k1) + a*sign(1,jm)*sign(1,km)  
      cmodes_i(i0,j1,k1)=cmodes_i(i0,j1,k1) + b*sign(1,jm)*sign(1,km)  

      cmodes_r(i1,j0,k0)=cmodes_r(i1,j0,k0) + a*sign(1,im)
      cmodes_i(i1,j0,k0)=cmodes_i(i1,j0,k0) + b*sign(1,im)

      cmodes_r(i1,j0,k1)=cmodes_r(i1,j0,k1) + a*sign(1,im)*sign(1,km)
      cmodes_i(i1,j0,k1)=cmodes_i(i1,j0,k1) + b*sign(1,im)*sign(1,km)

      cmodes_r(i1,j1,k0)=cmodes_r(i1,j1,k0) + a*sign(1,im)*sign(1,jm)
      cmodes_i(i1,j1,k0)=cmodes_i(i1,j1,k0) + b*sign(1,im)*sign(1,jm)

      cmodes_r(i1,j1,k1)=cmodes_r(i1,j1,k1) + a*sign(1,im)*sign(1,jm)*sign(1,km)
      cmodes_i(i1,j1,k1)=cmodes_i(i1,j1,k1) + b*sign(1,im)*sign(1,jm)*sign(1,km)

   enddo
   enddo
   enddo
enddo
enddo
enddo
end subroutine






subroutine complex_to_sincos_field(p,cmodes_r,cmodes_i)
#if 0
  convert set of complex FFT coefficients to complex coefficients

  For details of the 8x8 transform that maps sin/cos to complex
  coefficients, see the complex_to_sincos routine in sforcing.F90
#endif      
use params
implicit none
integer :: nmax
real*8 :: p(nx,ny,nz)
real*8 :: cmodes_r(nx,ny,nz)
real*8 :: cmodes_i(nx,ny,nz)
real*8 :: Rr, Ri
integer :: i,j,k,im,jm,km,ii,jj,kk,sm,i0,i1,j0,j1,k0,k1,ck

integer :: zerosign
external :: zerosign


p=0

do k=nz1,nz2,2
do j=ny1,ny2,2
do i=nx1,nx2,2

   i0=i
   i1=i+1
   j0=j
   j1=j+1
   k0=k
   k1=k+1
   
   ! verify that this processer has all 8 modes:
   if (  abs(imcord_exp(i0))/=abs(imcord_exp(i1)) .or. &
         abs(jmcord_exp(j0))/=abs(jmcord_exp(j1)) .or. & 
         abs(kmcord_exp(k0))/=abs(kmcord_exp(k1)) ) then
      print *,'not all modes are on this processor: '
      print *,i0,imcord_exp(i0),imcord_exp(i1)
      print *,j0,jmcord_exp(j0),jmcord_exp(j1)
      print *,k0,kmcord_exp(k0),kmcord_exp(k1)
      call abortdns(" ")
   endif

   ! make sure that i0 is the positive mode, i1 is the negative mode:
   if ( imcord_exp(i0)<0 .or. imcord_exp(i1)>0 .or.  &
        jmcord_exp(j0)<0 .or. jmcord_exp(j1)>0 .or.  & 
        kmcord_exp(k0)<0 .or. kmcord_exp(k1)>0 ) then
      print *,'we have the sign wrong?'
      print *,i0,imcord_exp(i0),imcord_exp(i1)
      print *,j0,imcord_exp(j0),imcord_exp(j1)
      print *,k0,imcord_exp(k0),imcord_exp(k1)
      call abortdns(" ")
   endif


   ! loop over the 8 sin/cos modes, mapping into 8 complex modes
   do ii=i0,i1
   do jj=j0,j1
   do kk=k0,k1      

      ! sin/cos mode 
      ! these arrays are the same as imcord(), except the last sin(N/2) wave
      ! number is zero instead of N/2.  Our complex modes do not have a N/2
      ! mode and store the sin(0) coefficient (which is always zero) 
      ! where the sin(N/2) mode is usually stored
      im=imcord_exp(ii)
      jm=jmcord_exp(jj)
      km=kmcord_exp(kk)

      Rr = cmodes_r(ii,jj,kk)
      Ri = cmodes_i(ii,jj,kk)


      p( i0, j0, k0) = p( i0, j0, k0) + Rr 
      p( i0, j0, k1) = p( i0, j0, k1) - Ri*zerosign(km) 
      p( i0, j1, k0) = p( i0, j1, k0) - Ri*zerosign(jm) 
      p( i0, j1, k1) = p( i0, j1, k1) - Rr*zerosign(jm*km) 
      p( i1, j0, k0) = p( i1, j0, k0) - Ri*zerosign(im)
      p( i1, j0, k1) = p( i1, j0, k1) - Rr*zerosign(im*km) 
      p( i1, j1, k0) = p( i1, j1, k0) - Rr*zerosign(im*jm) 
      p( i1, j1, k1) = p( i1, j1, k1) + Ri*zerosign(im*jm*km) 
   enddo
   enddo
   enddo
      
enddo
enddo
enddo

end subroutine





#if 0
phase shift algorithm

[0,2pi] domain:
delta = pi/N
let hh(k) = k delta  = k pi/N

[0,1] domain:
delta = 1/2N
let hh(k) = 2pi k delta = k pi/N


(a+ib) exp(ik (x + delta)) = (a+ib) exp(i hh) exp(ikx)
                          = (a+ib) (cos(hh) + i sin(hh)) exp(ikx)

double angle formula:
cos( kx + hh) + i sin( kx + hh) = exp(i (kx+hh)) = exp(i kx) exp( i hh)
   = ( cos(kx)+i sin(kx) ) ( cos(hh)+i sin(hh) )
   =  cos(kx) cos(hh) - sin(kx) sin(hh) + i [ sin(kx)cos(hh) + cos(kx) sin(hh) ]

a cos(kx + hh) = a cos(hh) cos(kx) - a sin(hh) sin(kx)
b sin(kx + hh) = b sin(hh) cos(kx) + b cos(hh) sin(kx)  

so for a pair of modes:  

a cos(kx+hh) + b sin(kx+hh) ==   [a cos(hh)+b sin(hh)]  cos(kx)
                               + [b cos(hh)-a sin(hh)]  sin(kx)

and thus the phase shift:
a -> a cos(hh)+b sin(hh)
b -> b cos(hh)-a sin(hh)

2D:
4 modes.  apply x:
a cos(kx) cos(jy)   ->  acos+bsin cos(kx) cos(jy)    a1
b sin(kx) cos(jy)   ->  bcos-asin cos(kx) cos(jy)    b1

c cos(kx) sin(jy)   ->  ccos+dsin cos(kx) sin(jy)    c1
d sin(kx) sin(jy)   ->  dcos-csin sin(kx) sin(jy)    d1

apply y:
a1 cos + c1 sin
b1 cos + d1 sin
c1 cos - a1 sin
d1 cos - b1 sin                          



#endif
subroutine z_phaseshift(p,shift,work)
!
!   shift = 1   apply phaseshift of .5 delta_x
!   shift =-1   apply phaseshift of -.5 delta_x  (inverse operation)
!
use params
implicit none
real*8 p(g_nz2,nx_2dz,ny_2dz)         
real*8 work(g_nz2,nx_2dz,ny_2dz)    
integer :: shift,i,j,k,im,jm,km
real*8 a,b,hh

do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nx_2dz
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)

         ! apply x shift
         ! if z_imsign(i)==0   then 
         !    im==0        do nothing (constant mode), hh = 0
         !    im=N/2 mode  so we only have the cosine component
         !                 we require this mode to be zero
         !                 so a=b=0 and shift will have no effect.

         hh = shift*pi*im/g_nx
         ! z_imsign(i) = 1:     a = cosine mode  
         ! im>0                 b = sine mode
         !                      p = a cos(hh) + b sin(hh)
         ! z_imsign(i) = -1:    b = cosine mode
         ! im<0                 a = sine mode
         !                      p = a cos(hh) - b sin(hh)
         a = p(k,i,j)               
         b = p(k,i+z_imsign(i),j) 
         work(k,i,j) = a*cos(hh) + b*sin(hh)
      enddo
   enddo
enddo

do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nx_2dz
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)
         ! apply y shift
         hh= shift*pi*jm/g_ny
         a = work(k,i,j)               
         b = work(k,i,j+z_jmsign(j))
         p(k,i,j) = a*cos(hh) + b*sin(hh)
      enddo
   enddo
enddo

do j=1,ny_2dz
   jm=z_jmcord(j)
   do i=1,nx_2dz
      im=z_imcord(i)
      do k=1,g_nz
         km=z_kmcord(k)
         ! apply z shift
         hh= shift*pi*km/g_nz
         a = p(k,i,j)               
         b = p(k+z_kmsign(k),i,j)
         work(k,i,j) = a*cos(hh) + b*sin(hh)
      enddo
   enddo
enddo

p=work
end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  compute hyperviscous dissipation term
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine hyperder(Q,nLaplace)

use params

implicit none

! input
real*8 time
integer compute_ints,rkstage

! input, but data can be trashed if needed
real*8 work(nx,ny,nz)                ! Fourier data at time t
real*8 Q(nx,ny,nz,3)                 ! grid data at time t
real*8 nLaplace(nx,ny,nz,3)

!local
real*8 xw,xw2,xw_viss
integer i,j,k,im,jm,km,n


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dissipation term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! take FFT of q, store result in nlaplace
nLaplace=q
do n=1,3
   call fft3d(nLaplace(:,:,:,n),work)
enddo

! compute viscous term, store result in nlaplace
do k=nz1,nz2
   km=kmcord(k)
   if (km==g_nz/2) km=0
   do j=ny1,ny2
      jm=jmcord(j)
      if (jm==g_ny/2) jm=0
      do i=nx1,nx2
         im=imcord(i)
         if (im==g_nx/2) im=0
         
         xw=(im*im + jm*jm + km*km/Lz/Lz)*pi2_squared
         xw_viss=mu*xw
         if (mu_hyper>=2 ) then
            xw2=(im*im*pi2_squared)
            xw2=xw2+(jm*jm*pi2_squared)
            xw2=xw2+(km*km*pi2_squared/(Lz*Lz))
            xw_viss=xw_viss + mu_hyper_value*(xw2/pi2_squared/k_Gt**2)**mu_hyper
         endif
         if (mu_hyper==0) then
            xw_viss=xw_viss + mu_hyper_value
         endif
         if (mu_hypo==1 .and. xw>0) then
            xw_viss=xw_viss + mu_hypo_value/xw
         endif
         
         nLaplace(i,j,k,1)=xw_viss*nlaplace(i,j,k,1)
         nLaplace(i,j,k,2)=xw_viss*nlaplace(i,j,k,2)
         nLaplace(i,j,k,3)=xw_viss*nlaplace(i,j,k,3)
         
      enddo
   enddo
enddo

! fft nlaplace back to grid point values
do n=1,3
   call ifft3d(nLaplace(:,:,:,n),work)
enddo


end subroutine



