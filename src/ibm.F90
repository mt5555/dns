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
real*8 :: surf_pt,surf_pt_max,surf_loc(5000,3)
real*8 :: rev_pt,rev_pt_max,rev_loc(5000,3,2)
real*8 :: cx(10),cy(10),cz(10)
real*8 :: jacobian(3,3),jacobinv(3,3)
real*8 :: jac_mat(5000,3,3)
real*8 :: jac_inv(5000,3,3)


contains

subroutine force(rhs,Q,p,time,compute_ints)
!
! calculate the forcing terms for the immersed boundary method and modify the rhs
!
! vor(1) = w_y - v_z
! vor(2) = u_z - w_x 
! vor(3) = v_x - u_y
!
! hel = u (w_y-v_z) + v (u_z - w_x)  + w (v_x - u_y)
!
use params
use fft_interface
implicit none

! input
real*8 Q(nx,ny,nz,n_var)
real*8 p(nx,ny,nz)
real*8 time
integer compute_ints

! output
real*8 rhs(nx,ny,nz,n_var)

!local
real*8 d(nx,ny,nz,3)
real*8 work(nx,ny,nz)
real*8 dummy(1)

integer n_body
integer i,j,k,l,m,n
integer ii,jj,kk
!integer surf_pt,surf_pt_max,surf_loc(5000,3)
!integer rev_pt,rev_pt_max,rev_loc(5000,3,2)
!real*8 cx(10),cy(10),cz(10)
real*8 rad,r,rdelt,bdelt,tempr
!real*8 jacobian(3,3),jacobinv(3,3)
!real*8 jac_mat(5000,3,3)
!real*8 jac_inv(5000,3,3)
real*8 usph(3),ucar(3)

! compute forcing terms

!set up forcing locations

cx(1)=32
cy(1)=32
cz(1)=32

cx(2)=48
cy(2)=16
cz(2)=48

cx(3)=48
cy(3)=48
cz(3)=16

cx(4)=16
cy(4)=48
cz(4)=48

!r=12
!rdelt=0.5
!bdelt=2.0
!n_body=1
!
!!loop over all bodies to set surface point (do only once if stationary)
!
!do n=1,3
!call der(p,d(1,1,1,n),dummy,work,1,n) 
!enddo
!
!do l=1,n_body
!   surf_pt=1
!   rev_pt=1
!
!   do k=cx(l)-r-1,cx(l)+r+1
!   do j=cy(l)-r-1,cy(l)+r+1
!   do i=cz(l)-r-1,cz(l)+r+1
!
!      rad=real((i-cx(l))*(i-cx(l))+(j-cy(l))*(j-cy(l))+(k-cz(l))*(k-cz(l)))
!      rad=sqrt(rad)
!
!! calculate location of surface points
!
!      if (rad.le.r+rdelt.and.rad.gt.r-rdelt) then 
!
!         surf_loc(l,surf_pt,1)=i
!         surf_loc(l,surf_pt,2)=j
!         surf_loc(l,surf_pt,3)=k
!         surf_pt=surf_pt+1
!
!! calculate location of reversal points
!
!      elseif (rad.le.r-rdelt.and.rad.gt.r-bdelt) then 
!
!         tempr=r-rad
!         rev_loc(l,rev_pt,1,1)=i
!         rev_loc(l,rev_pt,2,1)=j
!         rev_loc(l,rev_pt,3,1)=k
!
!         rev_loc(l,rev_pt,1,2)=nint(cx(l)+(i-cx(l))*(r+tempr)/(r-tempr))
!         rev_loc(l,rev_pt,2,2)=nint(cy(l)+(j-cy(l))*(r+tempr)/(r-tempr))
!         rev_loc(l,rev_pt,3,2)=nint(cz(l)+(k-cz(l))*(r+tempr)/(r-tempr))
!
!! calculate Jacobian for flow reversal points
!
!	 call jac(i-cx(l),j-cy(l),k-cz(l),jacobian)
!	 call jacinv(i-cx(l),j-cy(l),k-cz(l),jacobinv)
!
!	 do n=1,3
!	   do m=1,3
!	     jac_mat(l,rev_pt,n,m)=jacobian(n,m)
!	     jac_inv(l,rev_pt,n,m)=jacobinv(n,m)
!	   enddo
!	 enddo
!
!         rev_pt=rev_pt+1
!
!      endif
!
!      if(surf_pt.ge.5000) then 
!         write(*,*) 'Too many surface points',surf_pt
!         return
!      endif
!
!      if(rev_pt.ge.5000) then
!         write(*,*) 'Too many reversal points',rev_pt
!         return
!      endif
!
!   enddo
!   enddo
!   enddo
!
!   surf_pt_max=surf_pt-1
!   rev_pt_max=rev_pt-1
!
!enddo
!

write(*,*) surf_pt_max,rev_pt_max

! set surface velocity to zero

do l=1,n_body
   do surf_pt=1,surf_pt_max
 
      i=cx(l)+surf_loc(surf_pt,1)
      j=cy(l)+surf_loc(surf_pt,2)
      k=cz(l)+surf_loc(surf_pt,3)
 
      do n=1,3
 
!         rhs(i,j,k,n) = 0.
         rhs(i,j,k,n) = 0.-d(i,j,k,n)
         Q(i,j,k,n)=0.
!         Q(i,j,k,n)=0.+d(i,j,k,n)
 
      enddo
   
   enddo

!	set internal velocity
!	no-slip: reverse tangential, preserve normal

   do rev_pt=1,rev_pt_max
      i=cx(l)+rev_loc(l,rev_pt,1,1)
      j=cy(l)+rev_loc(l,rev_pt,2,1)
      k=cz(l)+rev_loc(l,rev_pt,3,1)
 
      ii=cx(l)+rev_loc(l,rev_pt,1,2)
      jj=cy(l)+rev_loc(l,rev_pt,2,2)
      kk=cz(l)+rev_loc(l,rev_pt,3,2)
 
      do n=1,3
        usph(n)=0
        ucar(n)=0
      enddo

      do n=1,3
        do m=1,3
          usph(n)=usph(n)+jac_inv(l,rev_pt,m,n)*Q(ii,jj,kk,m)
        enddo
      enddo

!      usph(1)=usph(1)
      usph(1)=-usph(1)
!      usph(1)=0.
      usph(2)=-usph(2)
      usph(3)=-usph(3)
 
      do n=1,3
        do m=1,3
          ucar(n)=ucar(n)+jac_mat(l,rev_pt,m,n)*usph(m)
        enddo
      enddo

      do n=1,3
 
!         rhs(i,j,k,n) = 0.
         rhs(i,j,k,n) = 0.-d(i,j,k,n)
         Q(i,j,k,n)=ucar(n)
!         Q(i,j,k,n)=-Q(ii,jj,kk,n)+d(i,j,k,n)
!         Q(i,j,k,n)=0
 
      enddo
      
   enddo

enddo

end subroutine

!	***************************************************************

       subroutine jac(nx,ny,nz,jj)

        implicit none

        integer nx,ny,nz
        integer i,j
        real*8 jj(3,3)
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

!       special case for phi=zero
!       otherwise matrix is singular

                do i=1,3
                do j=1,3
                        jj(i,j)=0.d0
                enddo
                enddo

!                jj(1,1)=dble(nx)/dble(iabs(nx))

!        write(*,*) nx,jj(1,1)

                jj(1,1)=1.d0
                jj(2,2)=1.d0
                jj(3,3)=1.d0
 
        endif
 
        return
        end subroutine

 
!       ***************************************************************
 
        subroutine jacinv(nx,ny,nz,jj)
 
        implicit none
 
        integer nx,ny,nz
        integer i,j
        real*8 jj(3,3)
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
 
!       special case for phi=zero
!       otherwise matrix is singular
 
                do i=1,3
                do j=1,3
                        jj(i,j)=0.d0
                enddo
                enddo
 
!                jj(1,1)=dble(nx)/dble(iabs(nx))

!        write(*,*) nx,jj(1,1)

                jj(1,1)=1.d0
                jj(2,2)=1.d0
                jj(3,3)=1.d0
 
        endif
 
        return
 
        end subroutine

end module
