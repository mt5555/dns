#include "macros.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  compute derivative along index index=1,2 or 3
!  numder = 1  compute p_x, return in px.   (pxx is not accessed)
!  numder = 2  compute p_xx, return in pxx.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine der(p,px,pxx,pt,numder,index)
use params
use fft_interface
use transform
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



subroutine divfree(u)
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
real*8 d1(nx,ny,nz)
real*8 work(nx,ny,nz)
real*8 p(nx,ny,nz)
real*8 :: dummy(1)
real*8 :: alpha=0
real*8 :: beta=1

integer i,j,k
real*8 divu(nx,ny,nz)
real*8 px(nx,ny,nz)
real*8 pxx(nx,ny,nz)
real*8 lap(nx,ny,nz)
real*8 lap2(nx,ny,nz)

! compute p = div(u)
i=1
call der(u(1,1,1,i),d1,dummy,work,1,i)
p = d1
i=2
call der(u(1,1,1,i),d1,dummy,work,1,i)
p = p + d1
i=3
call der(u(1,1,1,i),d1,dummy,work,1,i)
p = p + d1

divu=p

! solve laplacian(p)=div(u)
call poisson(p,work,alpha,beta)


! compute u=u-grad(p)
do i=1,3
   call der(p,d1,dummy,work,1,i)
   u(:,:,:,i) = u(:,:,:,i) - d1
enddo


end subroutine








subroutine poisson(f,work,alpha,beta)
!
!  solve laplacian(p) = f
!  input:  f 
!  ouput:  f   will be overwritten with the solution p
!
use params
use fft_interface
use transform
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
use transform
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
use transform
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_laplace_inverse(p,alpha,beta)
use params
implicit none
real*8 p(nx,ny,nz)
real*8 alpha,beta

integer i,j,k,im,jm,km
real*8 xfac

   do k=nz1,nz2
      if (k==nz1+1) then 
         km=(1+nz2-nz1)/2
      else
         km=(k-nz1)/2
      endif
      do j=ny1,ny2
         if (j==ny1+1) then
            jm=(1+ny2-ny1)/2
         else
            jm=(j-ny1)/2
         endif
         do i=nx1,nx2
            if (i==nx1+1) then
               im=(1+nx2-nx1)/2
            else
               im=(i-nx1)/2
            endif
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
! This mode has been moved to index = 2, taking the place of
! sin(0x) = 0 mode
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_filter(p)
use params
implicit none
real*8 p(nx,ny,nz)
real*8 alpha,beta

integer i,j,k,im,jm,km
real*8 xfac

   do k=nz1,nz2
!      km=(k-nz1)/2
      do j=ny1,ny2
!         jm=(j-ny1)/2
         do i=nx1,nx2
!            im=(i-nx1)/2
            if (k == nz1+1) p(i,j,k)=0
            if (j == ny1+1) p(i,j,k)=0
            if (i == nx1+1) p(i,j,k)=0
         enddo
      enddo
   enddo

end subroutine

