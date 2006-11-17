#include "macros.h"

#define COMPACT


subroutine helmholtz_dirichlet_inv(f,work,alpha,beta)
!
!  solve [alpha + beta*laplacian](p) = f
!  input:  f 
!  ouput:  f   will be overwritten with the solution p
!  b.c. for f specified in f. 
!
! THIS SOLVER USES dirichlet boundary conditions, NO MATTER WHAT
! how the code's boundary condtions (g_bdy_x1...) are set.  
!
use params
use fft_interface
use ghost
use transpose
implicit none
real*8 f(nx,ny,nz)    ! input/output
!real*8 work(nx,ny,nz) ! work array
real*8 work(g_ny2,nslabz,nx_2dy)
real*8 :: alpha
real*8 :: beta


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k,ix1,ix2,iy1,iy2
real*8 xm,ym,zm,xfac,im,jm
real*8 :: axy,x_axy,y_axy,diag
real*8 :: xb(ny,2),yb(nx,2)

real*8,save :: y_cosy(g_ny)
real*8,save :: y_cosx(g_nx) ! only needs to be of size nx_2dy, but that is not a parameter
logical,save :: firstcall=.true.

#if 0
! for checking b.c. tweak:
real*8 f2(nx,ny,nz)    
real*8 phi(nx,ny,nz)    
real*8 lphi(nx,ny,nz)    

f2=f
! phi = f on boundary, zero inside
phi=f2
call zero_boundary(phi)
phi=f2-phi
call helmholtz_dirichlet(phi,lphi,alpha,beta,work)
f2=f2-lphi
call zero_boundary(f2)
#endif

if (firstcall) then
   firstcall=.false.
   do j=2,g_ny
   do i=1,nx_2dy
      im=y_imsine(i)
      jm=j-1
      y_cosy(j)=cos(pi*dely*jm/yscale)
      y_cosx(i)=cos(pi*delx*im/xscale)
   enddo
   enddo
endif





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
#ifdef COMPACT
if (ndim/=2) then
   call abortdns("helmholtz_dirichlet_inv() not coded for 3D compact")
endif
axy=(delx**2+dely**2)/(12*delx*delx*dely*dely)

! this is only needed to compute helmholtz() near the boundary
! when converting problem to 0 b.c. problem (f=f-lphi requires
! ghost cell data to compute lphi for the COMPACT case
call ghost_update_x(f,1)
call ghost_update_y(f,1)

#else
axy=0

#endif


x_axy=1/(delx*delx) - 2*axy
y_axy=1/(dely*dely) - 2*axy
diag=-2/(delx*delx)-2/(dely*dely)+4*axy   




if (beta==0) then
   f=f/alpha
   return
endif
if (alpha/=0) then
   call abortdns("helmholtz_dirichlet_inv() cant yet handle alpha=0")
endif
!
!  let phi = f on the boundary, 0 inside
!      lphi = Helmholtz(phi)
!
!  solve:  Helmholtz(psi) = b - lphi  with psi=0 on boundary
!          then set psi=psi+phi
!
!
!NOTE: we dont arrays for phi, lphi:  
! save boundary data in single 1D array
! set f = 0 on boundary
! using boundary data alone: compute f = f - lphi 
!      along points just inside boundary
!
! solve as before
! restor boundary conditions in f from 1D array.




!
! compute f-lphi just above x1 boundary (uses values of f along x1 boundary).
! then zero out f along x1 boundary. 
! We need to set f=0 on the boundary, and doing this now will
! avoid double counting the contribution of the corner points when we do
! the y1 and y2 boundaries.
!
! NOTE: for this routine, we are solving a dirichlet problem NO MATTER
! WHAT type of boundary is being used.  
! We cant rely on bx1:bx2, or intx1:intx1.
!
iy1=by1; if (my_y==0) iy1=by1+1
iy2=by2; if (my_y==ncpu_y-1) iy2=by2-1
ix1=bx1; if (my_x==0) ix1=bx1+1
ix2=bx2; if (my_x==ncpu_x-1) ix2=bx2-1


k=1
if (my_x==0) then
   i=bx1
   do j=iy1,iy2
      f(i+1,j,k)=f(i+1,j,k) - ( axy*(f(i,j+1,k) + f(i,j-1,k)) + x_axy*f(i,j,k) )
   enddo
   do j=by1,by2
      xb(j,1)=f(i,j,k)
      f(i,j,k)=0
   enddo
endif
if (my_x==ncpu_x-1) then
   i=bx2
   do j=iy1,iy2
      f(i-1,j,k)=f(i-1,j,k) - (axy*(f(i,j+1,k) + f(i,j-1,k)) + x_axy*f(i,j,k) )
   enddo
   do j=by1,by2
      xb(j,2)=f(i,j,k)
      f(i,j,k)=0
   enddo
endif
! now do the y1 boundary:
if (my_y==0) then
   j=by1
   do i=ix1,ix2
      f(i,j+1,k)=f(i,j+1,k) - (axy*(f(i+1,j,k) + f(i-1,j,k)) + y_axy*f(i,j,k))
   enddo
   do i=bx1,bx2
      yb(i,1)=f(i,j,k)
      f(i,j,k)=0
   enddo
endif
! now do the y2 boundary:
if (my_y==ncpu_y-1) then
   j=by2
   do i=ix1,ix2
      f(i,j-1,k)=f(i,j-1,k) - (axy*(f(i+1,j,k) + f(i-1,j,k)) + y_axy*f(i,j,k))
   enddo
   do i=bx1,bx2
      yb(i,2)=f(i,j,k)
      f(i,j,k)=0
   enddo
endif




#if 0
k=1
do j=by1,by2
do i=bx1,bx2
   if (abs(f(i,j,k)-f2(i,j,k))>1e-6) then
      print *,bx1,bx2
      print *,by1,by2
      print *,i,j,f(i,j,k),f2(i,j,k)
      stop
   endif
enddo
enddo
#endif


! solve Helm(f) = b with 0 b.c.
#if 0
call transpose_to_z(f,work,n1,n1d,n2,n2d,n3,n3d)  
call sinfft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(work,f,n1,n1d,n2,n2d,n3,n3d)  
#endif

call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d) 
call sinfft1(work,n1,n1d,n2,n2d,n3,n3d)     
call transpose_from_x(work,f,n1,n1d,n2,n2d,n3,n3d) 

call transpose_to_y(f,work,n1,n1d,n2,n2d,n3,n3d)  ! x,y,z -> y,x,z
call sinfft1(work,n1,n1d,n2,n2d,n3,n3d)


do k=1,nslabz
do j=2,g_ny
do i=1,nx_2dy
   ! mode imsine(i),jmsine(j),kmsine(k)
   ! sin(x+h)sin(y+h) + sin(x-h)sin(y+h) + 
   ! sin(x+h)sin(y-h) + sin(x-h)sin(y-h) =
   ! 
   !   [  2*cos(h)  ]  sin(x) [sin(y+h)+sin(y-h) ]  = 
   !
   !   [  4*cos(hx) cos(hy) ]  sin(x) sin(y)

   im=y_imsine(i)
   jm=j-1
   ! y_cosy(j)=cos(pi*dely*jm/yscale)
   ! y_cosx(i)=cos(pi*delx*im/xscale)
   
   xfac = x_axy*2*y_cosx(i)                         ! x term
   xfac = xfac + y_axy*2*y_cosy(j)                  ! y term
   ! diagonal terms:
   xfac = xfac + axy*4*y_cosx(i)*y_cosy(j)
   ! center term:
   xfac = xfac + diag


   xfac = alpha + beta*xfac
   if (im+jm==0 ) then
      work(j,k,i) = 0
   else
      work(j,k,i) = work(j,k,i)/xfac
   endif
     
 
enddo
enddo
enddo


call isinfft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(work,f,n1,n1d,n2,n2d,n3,n3d)         ! y,x,z -> x,y,z

call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d) 
call isinfft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_x(work,f,n1,n1d,n2,n2d,n3,n3d )

#if 0
call transpose_to_z(f,work,n1,n1d,n2,n2d,n3,n3d)       
call isinfft1(work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(work,f,n1,n1d,n2,n2d,n3,n3d)       
#endif

! restore boundary values
! restore in the opposite order as saved above, because
! the corner points are not correct in yb(:,1) and yb(:,2)
k=1
if (my_y==ncpu_y-1) then
   j=by2
   do i=bx1,bx2
      f(i,j,k)=yb(i,2)
   enddo
endif
if (my_y==0) then
   j=by1
   do i=bx1,bx2
      f(i,j,k)=yb(i,1)
   enddo
endif
if (my_x==ncpu_x-1) then
   i=bx2
   do j=by1,by2
      f(i,j,k)=xb(j,2)
   enddo
endif
if (my_x==0) then
   i=bx1
   do j=by1,by2
      f(i,j,k)=xb(j,1)
   enddo
endif

end subroutine










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
integer i,j,k,setbdy,ix1,ix2,iy1,iy2

!
! NOTE: for this routine, we are solving a dirichlet problem NO MATTER
! WHAT type of boundary is being used.  
! We cant rely on bx1:bx2, or intx1:intx1 to compute the interior
!
iy1=by1; if (my_y==0) iy1=by1+1
iy2=by2; if (my_y==ncpu_y-1) iy2=by2-1
ix1=bx1; if (my_x==0) ix1=bx1+1
ix2=bx2; if (my_x==ncpu_x-1) ix2=bx2-1


if (ndim==2) then
   k=1

#ifdef COMPACT
   do j=inty1,inty2
   do i=intx1,intx2

      work(i,j,k)= &
            (                                       &
               (f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k)) + &
               (f(i,j+1,k)-2*f(i,j,k)+f(i,j-1,k))   &
             ) / 12

   enddo
   enddo
   do j=iy1,iy2
   do i=ix1,ix2
      f(i,j,k)=f(i,j,k)+work(i,j,k)
   enddo
   enddo
#endif

      
   if (setbdy==1) then
   ! now overwrite boundary with same data as in p
   if (my_x==0) then
      do j=by1,by2
         f(bx1,j,k)=p(bx1,j,k)
      enddo
   endif
   if (my_x==ncpu_x-1) then
      do j=by1,by2
         f(bx2,j,k)=p(bx2,j,k)
      enddo
   endif
   if (my_y==0) then
      do i=bx1,bx2
         f(i,by1,k)=p(i,by1,k)
      enddo
   endif
   if (my_y==ncpu_y-1) then
      do i=bx1,bx2
         f(i,by2,k)=p(i,by2,k)
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

   do k=bz1,bz2
   do j=by1,by2
   do i=bx1,bx2
      if ( (my_x==0 .and. i==bx1)  .or. &
           (my_x==ncpu_x-1 .and. i==bx2)  .or. &
           (my_y==0 .and. j==by1) .or. &
           (my_y==ncpu_y-1 .and. j==by2) .or.&
           (my_z==0 .and. k==bz1 .and. ndim==3) .or. &
           (my_z==ncpu_z-1 .and. k==bz2 .and. ndim==3) ) then
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
integer n,i,j,k,ix1,ix2,iy1,iy2,iz1,iz2
!
! NOTE: for this routine, we are solving a dirichlet problem NO MATTER
! WHAT type of boundary is being used.  
! We cant rely on bx1:bx2, or intx1:intx1 to compute the interior
!
ix1=bx1; if (my_x==0) ix1=bx1+1
ix2=bx2; if (my_x==ncpu_x-1) ix2=bx2-1

iy1=by1; if (my_y==0) iy1=by1+1
iy2=by2; if (my_y==ncpu_y-1) iy2=by2-1

iz1=bz1; if (my_z==0) iz1=bz1+1
iz2=bz2; if (my_z==ncpu_z-1) iz2=bz2-1



if (ndim==2) then
   call ghost_update_x(f,1)

   k=1

#ifdef COMPACT
   ! work = Dx(f)
   do j=by1,by2
   do i=bx1,bx2
      work(i,j,k)=(f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k))/(delx*delx)
   enddo
   enddo
   call ghost_update_y(work,1)
#endif

   call ghost_update_y(f,1)

   do j=by1,by2
   do i=bx1,bx2
      if ( i<ix1 .or. i>ix2 .or. j<iy1 .or. j>iy2) then
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
   do k=bz1,bz2
   do j=by1,by2
   do i=bx1,bx2
      if ( i<ix1 .or. i>ix2 .or. j<iy1 .or. j>iy2 .or. &
             k<iz1 .or. k>iz2) then
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
   do j=by1,by2
   do i=bx1,bx2
      work(i,j,k)=(f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k))/(delx*delx)
   enddo
   enddo
   call ghost_update_y(work,1)
#endif

   call ghost_update_y(f,1)

   do j=by1,by2
   do i=bx1,bx2
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
   do k=bz1,bz2
   do j=by1,by2
   do i=bx1,bx2
      lf(i,j,k)=alpha*f(i,j,k)+ &
            beta*(f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k))/(delx*delx) +&
            beta*(f(i,j+1,k)-2*f(i,j,k)+f(i,j-1,k))/(dely*dely) +&
            beta*(f(i,j,k+1)-2*f(i,j,k)+f(i,j,k-1))/(delz*delz) 
   enddo
   enddo
   enddo

endif


end subroutine




