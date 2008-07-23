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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  IBM (Immersed Bounday Method) Module
!  This module contains all the variables and subroutines for the 
!  immersed boundary method. Note that only ns_ibm uses this module. 
!  11/1/02 Jamaludin Mohd Yusof
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ibm
use params
implicit none

integer :: surf_pt,surf_pt_max
integer :: rev_pt,rev_pt_max
real*8 :: surf_loc(5000,3)
real*8 :: rev_loc(5000,3,2)
real*8 :: cx(10,3)
real*8 :: jacobian(3,3),jacobinv(3,3)
real*8 :: cd(3,5000)
real*8 :: jac_for_surf(5000,3,3)
real*8 :: jac_inv_surf(5000,3,3)
real*8 :: jac_for_rev(5000,3,3)
real*8 :: jac_inv_rev(5000,3,3)
integer :: n_body,n_bnd,n_rpr

contains

subroutine force(rhs,Q,p,time)
!
! calculate the forcing terms for the immersed boundary method and 
! modify the rhs
!
use params
use fft_interface
implicit none

! input
real*8 Q(nx,ny,nz,n_var)
real*8 p(nx,ny,nz)
real*8 time

! output
real*8 rhs(nx,ny,nz,n_var)

!local
real*8 d(nx,ny,nz,ndim)
real*8 work(nx,ny,nz)
real*8 dummy(nx,ny,nz)
real*8 xl(ndim),xlo(ndim),xli(ndim)
real*8 ul(ndim),ulo(ndim),uli(ndim)
real*8 t1(ndim),t1o(ndim),t1i(ndim)
real*8 interp(2,2,2)
real*8 usph(ndim),ucar(ndim)
real*8 drag(ndim),vort(ndim)
real*8 dd(ndim),pl
real*8 ibf(nx,ny,nz,ndim)
real*8 da

integer i,j,k,l,m,n
integer ii,jj,kk
real*8 rad,r,rdelt,bdelt,tempr
integer i_surf(ndim),i_int(ndim),i_ext(ndim)

! index conventions to keep all the loops somewhat readable
! n,m : step over dimension indices
! i,j,k : step in x,y,z respectively
! ii,jj,kk : step in x,y,z respectively
! l : step over bodies

! compute forcing terms

write(*,*) surf_pt_max,rev_pt_max

! calculate pressure gradient for use later

do m=1,ndim
   call der(p,d(1,1,1,m),dummy,work,1,m)
enddo

!-----------------------------------------------------------------------
!       zero the force array, ibf(nx,ny,nz,ndim)
 
  do n=1,ndim
    do k=1,nz
      do j=1,ny
        do i=1,nx
          ibf(i,j,k,n)=0.0
        enddo
      enddo
    enddo
  enddo
 
!-----------------------------------------------------------------------
!       loop over all bodies
 
n_body=1

cx(1,1)=32
cx(1,2)=32
cx(1,3)=32

do l=1,n_body

  do m=1,ndim
    drag(m)=0.
  enddo

  da = 4.*pi*64*delx*dely/n_bnd

  do surf_pt=1,surf_pt_max
    do n=1,ndim
      xl(n)=cx(l,n)+surf_loc(surf_pt,n)
      i_surf(n)=int(xl(n))
      t1(n)=xl(n)-i_surf(n)
      if (t1(n).gt.1.) write(*,*) 't1 overflow'
      if (t1(n).lt.0.) write(*,*) 't1 underflow'
      ul(n)=0.
    enddo
 
    write(*,*) surf_pt,i_surf(1),i_surf(2),i_surf(3)
 
    pl=0.
 
!  interpolate local velocity and pressure at surface points.

    do k=0,1
      do j=0,1
        do i=0,1
          interp(i+1,j+1,k+1)=((1-i)+(2*i-1)*t1(1))* &
                          ((1-j)+(2*j-1)*t1(2))* &
                          ((1-k)+(2*k-1)*t1(3))
 
          do n=1,ndim
 
            ul(n)=ul(n)+interp(i+1,j+1,k+1)* &
             Q(i_surf(1)+i,i_surf(2)+j,i_surf(3)+k,n)
 
          enddo
 
            pl=pl+interp(i+1,j+1,k+1)* &
             p(i_surf(1)+i,i_surf(2)+j,i_surf(3)+k)
 
        enddo
      enddo
    enddo
 
    do n=1,ndim
      usph(n)=0.
      ucar(n)=0.
    enddo
 
    write(*,*) surf_pt,ul(1),ul(2),ul(3)

    do m=1,ndim
      do n=1,ndim
        usph(m)=usph(m)+jac_inv_surf(surf_pt,n,m)*ul(n)
      enddo
    enddo

!-----------------------------------------------------------------------
!       set no-slip or no-penetration boundary conditions
 
!  no-slip

   usph(1)=0.
   usph(2)=0.
   usph(3)=0.
 
!  no-penetration

   usph(1)=0.
 
    do m=1,ndim
      do n=1,ndim
        ucar(m)=ucar(m)+jac_for_surf(surf_pt,n,m)*usph(n)
      enddo
    enddo
 
    do k=0,1
      do j=0,1
        do i=0,1
          do n=1,ndim
            ibf(i_surf(1)+i,i_surf(2)+j,i_surf(3)+k,n)= &
             ibf(i_surf(1)+i,i_surf(2)+j,i_surf(3)+k,n)+ &
             interp(i+1,j+1,k+1)*ul(n)

          enddo
        enddo
      enddo
    enddo
 
!-----------------------------------------------------------------------
!       calculate surface drag
 
    vort(1)=(Q(i_surf(1),i_surf(2)+1,i_surf(3),3) &
           -Q(i_surf(1),i_surf(2),i_surf(3),3))/dely &
          -(Q(i_surf(1),i_surf(2),i_surf(3)+1,2) &
           -Q(i_surf(1),i_surf(2),i_surf(3),2))/delz
    vort(2)=(Q(i_surf(1),i_surf(2),i_surf(3)+1,1) &
           -Q(i_surf(1),i_surf(2),i_surf(3),1))/delz &
          -(Q(i_surf(1)+1,i_surf(2),i_surf(3),3) &
           -Q(i_surf(1),i_surf(2),i_surf(3),3))/delx
    vort(3)=(Q(i_surf(1)+1,i_surf(2),i_surf(3),2) &
           -Q(i_surf(1),i_surf(2),i_surf(3),2))/delx &
          -(Q(i_surf(1),i_surf(2)+1,i_surf(3),1) &
           -Q(i_surf(1),i_surf(2),i_surf(3),1))/dely
 
    dd(1)=vort(3)*cd(2,l)-vort(2)*cd(3,l)
 
    dd(2)=vort(1)*cd(3,l)-vort(3)*cd(1,l)
 
    dd(3)=vort(2)*cd(1,l)-vort(1)*cd(2,l)
 
    do n=1,ndim
      drag(n)=drag(n)-dd(n)*da*mu+pl*cd(n,l)*da
    enddo

  enddo

! close loop over surface points

!goto 200

!-----------------------------------------------------------------------
!       loop over all reversal pairs

  do rev_pt=1,rev_pt_max

    do n=1,ndim

!-----------------------------------------------------------------------
!       location of internal points

      xli(n)=cx(l,n)+rev_loc(rev_pt,n,1)
      i_int(n)=int(xli(n))
      uli(n)=0.

!-----------------------------------------------------------------------
!       location of external points

      xlo(n)=cx(l,n)+rev_loc(rev_pt,n,2)
      i_ext(n)=int(xlo(n))
      ulo(n)=0.

!-----------------------------------------------------------------------
!       calculate external velocity
!       decomposed into u_normal, u_phi, u_theta
!       stored as   usph(1), usph(2), usph(3)
!-----------------------------------------------------------------------

      t1o(n)=xlo(n)-i_ext(n)
      t1i(n)=xli(n)-i_int(n)

    enddo

!   write(*,*) rev_pt,xli,xlo

!-----------------------------------------------------------------------
!       calculate external local velocities

    do k=0,1
      do j=0,1
        do i=0,1
          interp(i+1,j+1,k+1)=((1-i)+(2*i-1)*t1o(1))* &
                        ((1-j)+(2*j-1)*t1o(2))* &
                        ((1-k)+(2*k-1)*t1o(3))
          do n=1,ndim

            ulo(n)=ulo(n)+interp(i+1,j+1,k+1)* &
              Q(i_ext(1)+i,i_ext(2)+j,i_ext(3)+k,n)

          enddo
        enddo
      enddo
    enddo

!-----------------------------------------------------------------------
!       calculate internal local velocities

    do k=0,1
      do j=0,1
        do i=0,1
          interp(i+1,j+1,k+1)=((1-i)+(2*i-1)*t1i(1))* &
                        ((1-j)+(2*j-1)*t1i(2))* &
                        ((1-k)+(2*k-1)*t1i(3))
          do n=1,ndim

            uli(n)=uli(n)+interp(i+1,j+1,k+1)* &
              Q(i_int(1)+i,i_int(2)+j,i_int(3)+k,n)

          enddo
        enddo
      enddo
    enddo

    do n=1,ndim
      usph(n)=0.
      ucar(n)=0.
    enddo

    do n=1,ndim
      do m=1,ndim
        usph(n)=usph(n)+jac_inv_rev(rev_pt,m,n)*ulo(m)
      enddo
    enddo
 
!-----------------------------------------------------------------------
!       set internal forcing to reverse external velocity
!       reverse both tangential velocity components
!       zero radial component
!       scale tangential velocities by ratio of radii
 
    usph(1)=-usph(1)
!    usph(1)=0.
    usph(2)=0.8*usph(2)
    usph(3)=0.8*usph(3)
 
    do n=1,ndim
      do m=1,ndim
        ucar(n)=ucar(n)+jac_for_rev(rev_pt,m,n)*usph(m)
      enddo
    enddo
 
    do k=0,1
      do j=0,1
        do i=0,1
          do n=1,ndim
            ibf(i_int(1)+i,i_int(2)+j,i_int(3)+k,n)= &
              ibf(i_int(1)+i,i_int(2)+j,i_int(3)+k,n)+ &
              interp(i+1,j+1,k+1)*ul(n)
          enddo
        enddo
      enddo
    enddo
 
  enddo

200	continue

    do k=1,nz
      do j=1,ny
        do i=1,nx
          do n=1,ndim
!           rhs(i,j,k,n) = rhs(i,j,k,n)-ibf(i,j,k,n) 
            Q(i,j,k,n)=Q(i,j,k,n)-ibf(i,j,k,n)

!     rhs(i,j,k,n) = 0.
!     rhs(i,j,k,n) = 0.-d(i,j,k,n)
!     Q(i,j,k,n)=0.
!     Q(i,j,k,n)=0.+d(i,j,k,n)
 
          enddo
        enddo
      enddo
    enddo
   
enddo

end subroutine

!*************************************************************************

subroutine integer_points(rhs,Q,p,time)

! input
real*8 Q(nx,ny,nz,n_var)
real*8 p(nx,ny,nz)
real*8 time
 
! output
real*8 rhs(nx,ny,nz,n_var)

integer :: i,j,k,l,m,n
integer :: ii,jj,kk
real*8 dummy(1)
real*8 work(nx,ny,nz)
real*8 :: usph(3),ucar(3)
real*8 d(nx,ny,nz,ndim)

! calculate pressure gradient for use later

do m=1,ndim
   call der(p,d(1,1,1,n),dummy,work,1,m)
enddo

do l=1,n_body
  do rev_pt=1,rev_pt_max

!	set internal velocity
!	no-slip: reverse tangential, preserve normal

    i=cx(l,1)+rev_loc(rev_pt,1,1)
    j=cx(l,2)+rev_loc(rev_pt,2,1)
    k=cx(l,3)+rev_loc(rev_pt,3,1)
 
    ii=cx(l,1)+rev_loc(rev_pt,1,2)
    jj=cx(l,2)+rev_loc(rev_pt,2,2)
    kk=cx(l,3)+rev_loc(rev_pt,3,2)
 
    do n=1,ndim
      usph(n)=0
      ucar(n)=0
    enddo

    do n=1,ndim
      do m=1,ndim
        usph(n)=usph(n)+jac_inv_rev(rev_pt,m,n)*Q(ii,jj,kk,m)
      enddo
    enddo

!   set boundary conditions for no-slip or free-slip

!    usph(1)=usph(1)
    usph(1)=-usph(1)
!    usph(1)=0.

    usph(2)=-usph(2)
    usph(3)=-usph(3)
 
    do n=1,ndim
      do m=1,ndim
        ucar(n)=ucar(n)+jac_for_rev(rev_pt,m,n)*usph(m)
      enddo
    enddo

    do n=1,ndim
 
      rhs(i,j,k,n) = 0.-d(i,j,k,n)
      Q(i,j,k,n)=ucar(n)
 
    enddo
      
  enddo

enddo

end subroutine

!	***************************************************************

subroutine jacfor(nx,ny,nz,jj)

implicit none

integer i,j
real*8 nx,ny,nz
real*8 jj(ndim,ndim)
real*8 st,ct,sp,cp

real*8 r1,r2


r1=dsqrt(dble(ny*ny+nz*nz))
r2=dsqrt(dble(nx*nx+ny*ny+nz*nz))

if(r1.ne.(0.d0)) then

        ct=dble(ny)/r1
        st=dble(nz)/r1

        cp=dble(nx)/r2
        sp=r1/r2

        jj(1,1)=cp
        jj(1,2)=-sp
        jj(1,3)=0.d0

        jj(2,1)=sp*ct
        jj(2,2)=cp*ct
        jj(2,3)=-sp*st

        jj(3,1)=sp*st
        jj(3,2)=cp*st
        jj(3,3)=sp*ct

else

!       special case for phi=0.
!       otherwise matrix is singular

        do i=1,ndim
        do j=1,ndim
                jj(i,j)=0.d0
        enddo
        enddo

!        jj(1,1)=dble(nx)/dble(iabs(nx))

!        write(*,*) nx,jj(1,1)

        jj(1,1)=1.d0
        jj(2,2)=1.d0
        jj(3,3)=1.d0
 
endif
 
return
end subroutine
 
!***************************************************************
 
subroutine jacinv(nx,ny,nz,jj)

implicit none
 
integer i,j
real*8 nx,ny,nz
real*8 jj(ndim,ndim)
real*8 st,ct,sp,cp
 
real*8 r1,r2
 
 
r1=dsqrt(dble(ny*ny+nz*nz))
r2=dsqrt(dble(nx*nx+ny*ny+nz*nz))
 
if(r1.ne.(0.d0)) then
 
        ct=dble(ny)/r1
        st=dble(nz)/r1

        cp=dble(nx)/r2
        sp=r1/r2
 
        jj(1,1)=cp
        jj(1,2)=sp*ct
        jj(1,3)=sp*st
 
        jj(2,1)=-sp
        jj(2,2)=cp*ct
        jj(2,3)=cp*st
 
        jj(3,1)=0.d0
        jj(3,2)=-st/sp
        jj(3,3)=ct/sp
 
else
 
!       special case for phi=0.
!       otherwise matrix is singular
 
        do i=1,ndim
        do j=1,ndim
                jj(i,j)=0.d0
        enddo
        enddo
 
!        jj(1,1)=dble(nx)/dble(iabs(nx))

!        write(*,*) nx,jj(1,1)

        jj(1,1)=1.d0
        jj(2,2)=1.d0
        jj(3,3)=1.d0
 
endif
 
return
 
end subroutine

!*************************************************************************************

subroutine node_set

integer i,j,k,l,m,n
integer ii,jj,kk
real*8 rad,r,rdelt,bdelt,tempr

! compute forcing terms

!set up forcing locations
!list of body centers 

cx(1,1)=32
cx(1,2)=32
cx(1,3)=32

cx(2,1)=48
cx(2,2)=16
cx(2,3)=48

cx(3,1)=48
cx(3,2)=48
cx(3,3)=16

r=12
rdelt=0.5
bdelt=2.0
n_body=1

!loop over all bodies to set surface point (do only once if stationary)

do l=1,n_body
   surf_pt=1
   rev_pt=1

   do k=cx(l,3)-r-1,cx(l,3)+r+1
   do j=cx(l,2)-r-1,cx(l,2)+r+1
   do i=cx(l,1)-r-1,cx(l,1)+r+1

      rad=real((i-cx(l,1))*(i-cx(l,1))+(j-cx(l,2))*(j-cx(l,2))+(k-cx(l,3))*(k-cx(l,3)))
      rad=sqrt(rad)

! calculate location of surface points

      if (rad.le.r+rdelt.and.rad.gt.r-rdelt) then 

         surf_loc(surf_pt,1)=i
         surf_loc(surf_pt,2)=j
         surf_loc(surf_pt,3)=k
         surf_pt=surf_pt+1

! calculate location of reversal points

      elseif (rad.le.r-rdelt.and.rad.gt.r-bdelt) then 

         tempr=r-rad
         rev_loc(rev_pt,1,1)=i
         rev_loc(rev_pt,2,1)=j
         rev_loc(rev_pt,3,1)=k

         rev_loc(rev_pt,1,2)=nint(cx(l,1)+(i-cx(l,1))*(r+tempr)/(r-tempr))
         rev_loc(rev_pt,2,2)=nint(cx(l,2)+(j-cx(l,2))*(r+tempr)/(r-tempr))
         rev_loc(rev_pt,3,2)=nint(cx(l,3)+(k-cx(l,3))*(r+tempr)/(r-tempr))

! calculate Jacobian for flow reversal points

	 call jacfor(i-cx(l,1),j-cx(l,2),k-cx(l,3),jacobian)
         call jacinv(i-cx(l,1),j-cx(l,2),k-cx(l,3),jacobinv)

	 do n=1,3
	   do m=1,3
	     jac_for_rev(rev_pt,n,m)=jacobian(n,m)
	     jac_inv_rev(rev_pt,n,m)=jacobinv(n,m)
	   enddo
	 enddo

         rev_pt=rev_pt+1

      endif

      if(surf_pt.ge.5000) then 
         write(*,*) 'Too many surface points',surf_pt
         return
      endif

      if(rev_pt.ge.5000) then
         write(*,*) 'Too many reversal points',rev_pt
         return
      endif

   enddo
   enddo
   enddo

   surf_pt_max=surf_pt-1
   rev_pt_max=rev_pt-1

enddo

end subroutine


end module
