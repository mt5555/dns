subroutine ns3d(Q,rhs)
!
! vor(1) = w_y - v_z
! vor(2) = u_z - w_x 
! vor(3) = v_x - u_y
!

implicit none
use params

! input
real*8 Q(nxd,nyd,nzd,nvard)
! output
real*8 rhs(nxd,nyd,nzd,nvard)
!local
real*8 d1(nxd,nyd,nzd)
real*8 d2(nxd,nyd,nzd)
real*8 vor(nxd,nyd,nzd,3)     ! could be removed and data accumulated in rhs
integer i,j,k


rhs=0
vor=0
! compute viscous terms (in rhs) and vorticity
do i=1,3
   ! compute u_x, u_xx
   der(Q(1,1,1,i),d1,d2,2,1)
   rhs(:,:,:,i) += d2
   if (i==3) vor(2) += -d1
   if (i==2) vor(3) += -d1

   ! compute u_y, u_yy
   der(Q(1,1,1,i),d1,d2,2,2)
   rhs(:,:,:,i) += d2
   if (i==3) vor(1) += d1
   if (i==1) vor(3) += -d1

   ! compute u_z, u_zz
   der(Q(1,1,1,i),d1,d2,2,3)
   rhs(:,:,:,i) += d2
   if (i==2) vor(:,:,:,1) += -d1
   if (i==1) vor(:,:,:,2) += d1
enddo



! build up the rhs 
vor = q cross vor
rhs = mu*rhs + vor

! vor = divergence (q cross vor)
d2=0
do i=1,3
   der(vor(i),d1,dummy,1,i)
   d2 += d1
enddo

! solve laplacian d1 = d2
poisson(d1,d2)
do i=1,3
   der(d1,d2,dummy,1,i)
   rhs(:,:,:,i) += d2
enddo

end











