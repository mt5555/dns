subroutine init_data(Q)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)

! local variables
integer i,j,k,l
real*8 delta,pi2,delsq,delalf,delgam,yval,xval,dify,difx,uu,vv,denom
real*8 xscale,yscale
real*8 :: eps=.10
integer :: km=1
integer,parameter :: n=500
real*8 :: x(0:n),y(0:n)

! uniform flow to the right
Q=0


delta = .05
pi2 = 2*pi
delsq = delta**2

! Initialize vortex sheet
delalf = 2d0/n
delgam = delalf/km         !SO ALL PERIODS HAVE TOTAL CIRC=1


xscale=2
yscale=4

do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
  
   xval=xscale*(xcord(i)-.5)         ! x ranges from -1 .. 1
   if (ycord(j)<=.5) then
      yval=yscale*(ycord(j)-.25)  ! bottom 50% goes from -.1 .. .1
   else
      yval=-yscale*(ycord(j)-.75) ! top 50% goes from .1 .. -.1
      xval=-xval	
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

   print *,'B'
   stop

   Q(i,j,k,1) = 5*uu*delgam/(pi2*xscale)
   Q(i,j,k,2) = 5*vv*delgam/(pi2*yscale)

   print *,'bottom of loop'

enddo
enddo
enddo

end subroutine







subroutine init_data_projection(Q)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)

! local variables
integer i
real*8 :: d1(nx,ny,nz)
real*8 :: d2(nx,ny,nz)
real*8 :: d3(nx,ny,nz)

call bc_preloop(Q)

! will also dealias if dealias=1
call divfree(Q,d1,d2,d3)


end subroutine





