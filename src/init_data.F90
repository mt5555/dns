

subroutine init_data_lwisotropic(Q)
!
! low wave number, quasi isotropic initial condition
!
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: work(nx,ny,nz)
integer km,jm,im,i,j,k,n,wn
real*8 xw,ener_target,ener,xfac

type wnforcing_d
   real*8  :: wn
   integer :: n
   integer, allocatable :: index(:,:)
end type
type(wnforcing_d) :: wnforcing(2)

Q=0

do n=1,2
   wnforcing(n)%wn=n
   wnforcing(n)%n=0
enddo

! count the number of wavenumbers in each band.  
do k=nz1,nz2
   km=kmcord(k)
   do j=ny1,ny2
      jm=jmcord(j)
      do i=nx1,nx2
         im=imcord(i)
         xw=sqrt(real(km**2+jm**2+im**2))
         if (xw>=.5 .and. xw<1.5) then
            wnforcing(1)%n=wnforcing(1)%n+1
         endif
         if (xw>=1.5 .and. xw<2.5) then
            wnforcing(2)%n=wnforcing(2)%n+1
         endif
      enddo
   enddo
enddo

! allocate storage
do n=1,2
   i=wnforcing(n)%n
   allocate(wnforcing(n)%index(i,3))
   print *,'wave no=',wnforcing(n)%wn,' number of coeffs=',wnforcing(n)%n
   wnforcing(n)%n=0  ! reset counter to use again below
enddo

! store all the indexes
do k=nz1,nz2
   km=kmcord(k)
   do j=ny1,ny2
      jm=jmcord(j)
      do i=nx1,nx2
         im=imcord(i)
         xw=sqrt(real(km**2+jm**2+im**2))
         if (xw>=.5 .and. xw<1.5) then
            wnforcing(1)%n=wnforcing(1)%n+1
            wnforcing(1)%index(wnforcing(1)%n,1)=i
            wnforcing(1)%index(wnforcing(1)%n,2)=j
            wnforcing(1)%index(wnforcing(1)%n,3)=k
         endif

         if (xw>=1.5 .and. xw<2.5) then
            wnforcing(2)%n=wnforcing(2)%n+1
            wnforcing(2)%index(wnforcing(2)%n,1)=i
            wnforcing(2)%index(wnforcing(2)%n,2)=j
            wnforcing(2)%index(wnforcing(2)%n,3)=k
         endif
      enddo
   enddo
enddo


do wn=1,2
   ener_target=(wnforcing(wn)%wn)**(-5.0/3.0)
   do 
      !random initial condition
      !pick a random point in the n dimensional sphere
      ! 1. first pick a randome point in R^n
      ! 2. through away points outside the sphere.  
      ! 3. When we get a point inside, project onto the surface.    
      ener=0
      do n=1,6 !wnforcing(wn)%n
         i=wnforcing(1)%index(n,1)
         j=wnforcing(1)%index(n,2)
         k=wnforcing(1)%index(n,3)
         km=kmcord(k)
         jm=jmcord(j)
         im=imcord(i)
         call random_number(Q(i,j,k,1))
         call random_number(Q(i,j,k,2))
         call random_number(Q(i,j,k,3))
         Q(i,j,k,:)=2*Q(i,j,k,:)-1

         xfac = 2*2*2*(g_nx*g_ny*g_nz)
         xfac=1
!         if (km==0) xfac=xfac/2
!         if (jm==0) xfac=xfac/2
!         if (im==0) xfac=xfac/2
         ener=ener + xfac*(Q(i,j,k,1)**2 + Q(i,j,k,2)**2 + Q(i,j,k,3)**2)
      enddo
!MPI SUM
      ener=sqrt(ener)
      print *,ener,wnforcing(wn)%n
      if (ener<=1 .and. ener>.000001) exit
   enddo
   stop
   ! normalize
   do n=1,wnforcing(wn)%n
      i=wnforcing(1)%index(n,1)
      j=wnforcing(1)%index(n,2)
      k=wnforcing(1)%index(n,3)
      Q(i,j,k,:)=Q(i,j,k,:)*sqrt(ener_target/ener)
   enddo

enddo

do n=1,3
   call ifft3d(Q(1,1,1,n),work)  ! use rhs as a work array
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

k=nz1
eps=200
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


do k=nz1+1,nz2
do i=nx1,nx2
do j=ny1,ny2
   Q(i,j,k,1)=Q(i,j,nz1,1)	
   Q(i,j,k,2)=Q(i,j,nz1,2)	
enddo
enddo
enddo



end subroutine



subroutine init_data_khblob(Q)
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

call bc_preloop(Q)

! will also dealias if dealias=1
call divfree(Q,d1)


end subroutine





