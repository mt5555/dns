subroutine init_data_test(Q)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: work(nx,ny,nz)
integer km,jm,im,i,j,k,n
real*8 xw,R,theta

Q=0
do n=1,3
call fft3d(Q(1,1,1,n),work)

   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)
            xw=sqrt(real(km**2+jm**2+im**2))
            if (xw>0 .and. xw<=2.5) then
               R=xw**(-5/3.0)
               call random_number(theta)
               theta=pi2*theta
               !  a_lmn 
            endif
         enddo
      enddo
   enddo

call ifft3d(Q(1,1,1,n),work)
enddo

end subroutine






subroutine init_data_kh(Q)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)

! local variables
integer i,j,k,l
real*8 delta,delsq,delalf,delgam,yval,xval,dify,difx,uu,vv,denom
real*8 xscale,yscale
real*8 :: eps=.10
integer :: km=1
integer,parameter :: n=500
real*8 :: x(0:n),y(0:n)

Q=0
#undef NEWKH
#ifdef  NEWKH

k=1
eps=30
do j=ny1,ny2
do i=nx1,nx2
   if (ycord(j)<=.5) then
      Q(i,j,k,1)=tanh(eps*(ycord(j)-.25))
   else
      Q(i,j,k,1)=tanh(eps*(.75-ycord(j)))
   endif
   Q(i,j,k,2)=.05*sin(2*pi*xcord(i))
enddo
enddo



#else

delta = .05
delsq = delta**2

! Initialize vortex sheet
delalf = 2d0/n
delgam = delalf/km         !SO ALL PERIODS HAVE TOTAL CIRC=1


xscale=2
yscale=4

k=1
do j=ny1,ny2
do i=nx1,nx2
  
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
      x(l) = -1 + l*delalf + xval       ! x ranges from xval-1 .. xval+1
      y(l) = eps*sin( km*pi*x(l) )

      difx =  xval - x(l) 
      dify =  yval - y(l) 
      denom = difx**2 + dify**2 + delsq
      uu = uu - dify/denom
      vv = vv + difx/denom

   enddo

   Q(i,j,k,1) = 5*uu*delgam/(pi2*xscale)
   Q(i,j,k,2) = 5*vv*delgam/(pi2*yscale)
enddo
enddo

#endif

do k=nz1+1,nz2
do i=nx1,nx2
do j=ny1,ny2
   Q(i,j,k,1)=Q(i,j,nz1,1)	
   Q(i,j,k,2)=Q(i,j,nz1,2)	
enddo
enddo
enddo



end subroutine







subroutine init_data_projection(Q)
use params
use fft_interface
implicit none
real*8 :: Q(nx,ny,nz,n_var)

! local variables
integer i
real*8 :: d1(nx,ny,nz)
real*8 :: d2(nx,ny,nz)
real*8 :: d3(nx,ny,nz)

call bc_preloop(Q)

! will also dealias if dealias=1
call divfree(Q,d1)


end subroutine





