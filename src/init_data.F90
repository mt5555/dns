subroutine init_data(Q)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: d1(nx,ny,nz)
real*8 :: d2(nx,ny,nz)
real*8 :: d3(nx,ny,nz)

! local variables
integer i,j,k,l
real*8 delta,pi2,delsq,delalf,delgam,yval,xval,dify,difx,uu,vv,denom
real*8 :: eps=.10
integer :: km=1
integer,parameter :: n=500
real*8 :: x(n),y(n)

! uniform flow to the right
Q=0


delta = .10
pi2 = 2*pi
delsq = delta**2

! Initialize vortex sheet
delalf = 2d0/n
delgam = delalf/km         !SO ALL PERIODS HAVE TOTAL CIRC=1

do i=0,n
enddo

do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2

   xval=2*(xcord(i)-.5)         ! x ranges from -1 .. 1
   if (ycord(j)<=.5) then
      yval=4*(ycord(j)-.25)  ! bottom 50% goes from -.5 .. .5
   else
      yval=-4*(ycord(j)-.75) ! top 50% goes from .5 .. -.5
   endif



   uu = 0
   vv = 0
   do l=0,n

      ! COMPUTE VELO AT (XCORD(I),YCORD(J)) INDUCED BY BLOB AT (X(L),Y(L))
!      x(l) = -1 + l*delalf            ! x ranges from -1 .. 1
      x(l) = -1 + l*delalf + xval       ! x ranges from xval-1 .. xval+1
      y(l) = eps*sin( km*pi*x(l) )

      difx =  xval - x(l) 
      dify =  yval - y(l) 
      denom = difx**2 + dify**2 + delsq
      uu = uu - dify/denom
      vv = vv + difx/denom
   enddo

   Q(i,j,k,1) = uu/pi2*delgam
   Q(i,j,k,2) = vv/pi2*delgam

enddo
enddo
enddo


! remove that pesky highest cosine mode
do i=1,3
   call fft3d(Q(1,1,1,i),d1)
   call fft_filter(Q(1,1,1,i))
   call ifft3d(Q(1,1,1,i),d1)
enddo
call divfree(Q,d1,d2,d3)

call bc_preloop


end subroutine









