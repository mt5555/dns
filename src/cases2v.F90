#include "macros.h"
!
! test cases that require n_var >= 2
!


subroutine init_data_zero(Q,q1,work1,work2)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
integer :: i,j,k,n

Q=0

if (init_cond_subtype==1) then
do k=1,nz
do j=1,ny
do i=1,nx
   Q(i,j,k,1)=1
   do n=2,n_var
      Q(i,j,k,n)=0
   enddo
enddo
enddo
enddo
endif


end subroutine




subroutine init_data_kh(Q,q1,work1,work2)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

! local variables
integer :: i,j,k,l,use3d=0
integer :: kinit(50)
real*8 :: eps
real*8 :: amp

Q=0


if (init_cond_subtype==0) then
   call print_message("Using thin shear layer initial condition")
   eps=200
   amp=.05
else if (init_cond_subtype==1) then
   ! E & Liu case:
   call print_message("Using E & Liu shear layer initial condition")
   eps=10*pi
   amp=.25
else if (init_cond_subtype==3) then
   use3d=1
   call print_message("Using 3d shear layer initial condition")
   eps=200
   amp=.05
else if (init_cond_subtype==4) then
   call print_message("Using 5 mode thin shear layer initial condition")
   eps=100
   amp=.10
endif

! thickness = 1/eps
! gridpoints per transition layer: nx/eps
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   if (ycord(j)<=.5) then
      Q(i,j,k,1)=tanh(eps*(ycord(j)-.25))
   else
      Q(i,j,k,1)=tanh(eps*(.75-ycord(j)))
   endif
   if (use3d==1) then
      Q(i,j,k,2)=amp*sin(2*pi*xcord(i))*cos(2*pi*zcord(k))
   else
      Q(i,j,k,2)=amp*sin(2*pi*xcord(i))
   endif
enddo
enddo
enddo

   


end subroutine



subroutine init_data_khblob(Q,q1,work1,work2)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

! local variables
integer i,j,k,l,i2
real*8 delta,delsq,delalf,delgam,yval,xval,dify,difx,uu,vv,denom
real*8 xsc,ysc
real*8 :: eps=.10
integer :: km=1
integer,parameter :: n=500
real*8 :: x(0:n),y(0:n),xtmp

Q=0

delta = .05
delsq = delta**2

! Initialize vortex sheet
delalf = 2d0/n
delgam = delalf/km         !SO ALL PERIODS HAVE TOTAL CIRC=1


xsc=2
ysc=4

k=1
do j=ny1,ny2
do i=nx1,nx2
  
   xval=xsc*(xcord(i)-.5)         ! x ranges from -1 .. 1
   if (ycord(j)<=.5) then
      yval=ysc*(ycord(j)-.25)  ! bottom 50% goes from -.1 .. .1
   else
      yval=-ysc*(ycord(j)-.75) ! top 50% goes from .1 .. -.1
      xval=-xval	
   endif

   uu = 0
   vv = 0
   do l=0,n
      ! COMPUTE VELO AT (XCORD(I),YCORD(J)) INDUCED BY BLOB AT (X(L),Y(L))
      x(l) = -1 + l*delalf + xval       ! x ranges from xval-1 .. xval+1
      y(l) = eps*sin( km*pi*x(l) )
      if (init_cond_subtype==1) then
         y(l)=0
         do i2=1,11,10
            y(l) = y(l) + (eps/i2)*sin( i2*pi*x(l) )
         enddo
      endif

      difx =  xval - x(l) 
      dify =  yval - y(l) 
      denom = difx**2 + dify**2 + delsq
      uu = uu - dify/denom
      vv = vv + difx/denom

   enddo

   Q(i,j,k,1) = 5*uu*delgam/(pi2*xsc)
   Q(i,j,k,2) = 5*vv*delgam/(pi2*ysc)
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









subroutine init_data_vxpair(Q,Qhat,work1,w,init)
!
!  init=0    setup parameters only, then exit
!  init=1    setup parameters and compute initial condition + b.c.
!                                 initialize tracers
!  init=2    setup parameters, initialize b.c. (not saved in restart file)
!                                 read tracer restart file
!
use params
use tracers
use ellipse
use bc
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: w(nx,ny,nz)
integer :: init

! local variables
integer :: i,j,k,l,use3d=0,n,initial_vor
real*8 :: eps
real*8 :: amp,vx,uy
real*8 :: delsq,hold,difx,dify,sumy,denom1,denom2,wsum,psisum,delalf,testmax
real*8 :: delta,xlocation
integer,parameter :: nd=200
real*8 :: xd(0:nd),yd(0:nd),wd(0:nd)
real*8 :: a0,a1,a2,a3,b0,b1,b2,b3,yalf,gamp


a0 =  -0.86887373086487
a1 =   22.19093784381875
a2 =  -52.31046126329605
a3 =   34.05681079318104 

b0 =  1.4
b1 = 0 
b2 = 20 
b3 =  -2000.0/45.0


offset_bdy=1
xlocation=-1   ! not yet set
initial_vor=0

if (init_cond_subtype ==0) then
!  my standard test case:
   biotsavart_cutoff=.001  
   biotsavart_ubar=.089
   biotsavart_apply=5  ! disabled
   yscale=3
   delta=.2

   if (g_nx==65) then
      print *,'disabling offset'
      offset_bdy=0
   endif
endif



if (init_cond_subtype ==1) then
! Monika's case
! 660x380   290min.  t=30  mu=1e-4
! [-1.7, 1.6]  [0 1.9]   delx=dely=.005
!OR:
! [ 0 , 3.3 ]  sheet at .1.7         xlocation=1.7/3.3
!
! scaled up to 720x400:   yscale=2
!                         xscale=3.6   
!
   biotsavart_cutoff=5e-4
   biotsavart_apply=-1  ! disabled
   delta=.2
   biotsavart_ubar=.089
   yscale=2
endif



if (init_cond_subtype ==2) then
! the standard initial condition for delta=.2
   biotsavart_cutoff=5e-3
   biotsavart_apply=50
   delta=.2
   biotsavart_ubar=.089
   yscale=2
endif


if (init_cond_subtype ==3) then
! special quick init case for performance
   biotsavart_cutoff=5e-2
   biotsavart_apply=-2  ! disabled, quick initialization
   delta=.2
   biotsavart_ubar=.089
   yscale=2
endif



if (init_cond_subtype ==4) then
   ! tweaked version of type 2, to run a little faster
   ! used for highest resolution run.
   biotsavart_cutoff=5e-3
   biotsavart_apply=100
   delta=.2
   biotsavart_ubar=.089
   yscale=1.9
endif

if (init_cond_subtype ==5) then
! the standard initial condition for delta=.1
   biotsavart_cutoff=5e-3
   biotsavart_apply=100
   delta=.1
   biotsavart_ubar=.097
   yscale=1.8
   xlocation=1.5
endif


if (init_cond_subtype ==100) then
   ! Kras iniital condition.  run to t=5, with nu=1e-6
   ! x:  -2 ... 0.5          5*32      25% more grid points in x
   ! y:  0.. 2               4*32
   ! [5N,4N] = [640,512]  
   biotsavart_cutoff=5e-3
   biotsavart_apply=10
   delta=.1
   biotsavart_ubar=.000
   yscale=2.0

   initial_vor = 1
endif


! set xscale so that delx=dely
xscale = yscale*(g_nx+offset_bdy-1)/(g_ny+offset_bdy-1)
call init_grid()   ! redo grid points since we changed scalings


! xlocation not set.  default value:
if (xlocation<0) xlocation=2.0



! initialize the crossing and vxline points:
vxline_count=9
do i=1,9
   vxline_y(i)=.1*i
enddo

if (init==0) return



Q=0



! INITIALIZES VORTEX SHEET (XS,YS,WS,NS)
if (initial_vor==0) then
   delalf = pi/(2*nd)
   do k=0,nd
      hold = k*delalf
      xd(k)  = xlocation
      yd(k)  = cos(hold)
      wd(k)  = cos(hold)*delalf
   enddo
   wd(nd) = wd(nd)/2
   wd(0) = wd(0)/2
   yd(nd) = 0
else if (initial_vor==1) then
   delalf = pi/(2*nd)
   do k=0,nd
      hold = k*delalf
      xd(k)  = xlocation
      yd(k)  = cos(hold)
      yalf = -sin(hold)
      if (yd(k)<.3) then
         gamp  = b1 + 2*b2*yd(k) + 3*b3*yd(k)*yd(k)
      else if (yd(k)<.7) then
         gamp  = a1 + 2*a2*yd(k) + 3*a3*yd(k)*yd(k)
      else 
         gamp  = yalf  ! -y/sqrt(1-y^2)
      endif
      wd(k)  = gamp*(-sin(hold))*delalf
   enddo
   if (io_pe==my_pe) then
      !do k=1,nd
      !   write(11,'(i4,2f20.7)') k,yd(k),wd(k)
      !enddo
      !close(11)
   endif
   wd(nd) = wd(nd)/2
   wd(0) = wd(0)/2
   yd(nd) = 0
else
   call abort("error vxpair(): bad value of initial_vor")
endif

delsq = delta**2
!   INITIALIZES VORTICITY ON THE GRID = VORTICITY INDUCED
!     BY A BLOB AT (0,1)
do j=by1,by2
   if (mod(j,100)==0 .and. my_pe==io_pe) then
      print *,'initial condition: ',inty2,j
   endif
do i=bx1,bx2
   wsum = 0
   psisum=0
   do k=0,nd
      difx = (xcord(i) - xd(k))
      dify = (ycord(j) - yd(k))
      sumy = (ycord(j) + yd(k))

      denom1 = difx**2 + dify**2 + delsq
      denom2 = difx**2 + sumy**2 + delsq
      wsum = wsum + wd(k)*(1/denom1**2-1/denom2**2)
      psisum = psisum - wd(k)*log(denom1/denom2)
   enddo
   w(i,j,1) = delsq*wsum/pi
   Qhat(i,j,nz1,2)= psisum*delx*dely/(4*pi)  - biotsavart_ubar*ycord(j)
enddo
enddo
call comp_ellipse_reshape(w,1,1)  ! use initial w to set max vorticity for ellipse contours

! We are generating an initial condition, so store w in Qhat(:,:,:,1)
! otherwise, initial condition was read in, and is already in Qhat(:,:,:,1)
if (init==1) Qhat(:,:,:,1)=w

! do this even for restart case, because boundary conditions were not
! saved in restart file if offset_bdy==1
call bcw_impose(Qhat(1,1,1,1))

! Set the values of PSI on the boundary coming from the initial condition,
! if needed.  
if (biotsavart_apply==-1) then
   ! psi on the boundary fixed for all time from the initial conditon:
   call print_message("Setting PSI boundary values (fixed for all time)")
   call bc_biotsavart(w,Qhat(1,1,1,2),1)  ! recompute PSI on boundary from w
else if (biotsavart_apply==-2) then
   ! short cut for PSI boundary values.  used for quicker testing:
   ! psi on the boundary fixed for all time from the initial conditon:
   call print_message("Setting PSI boundary values (fixed for all time)")
   call set_biotsavart(Qhat)    ! set based on PSI computed above 
else
   ! psi on the boundary always computed from w in time stepping loop.
endif




if (init==1) then
   numt_insert=500              ! 500 particles with insertion
   n=numt_insert+vxline_count   ! 10 extra particles with no insertion between them
   call allocate_tracers(n)
   delalf = pi/(2*(numt_insert-1))
   do k=1,numt_insert
      hold=(k-1)*delalf
      tracer(k,1)=xlocation
      tracer(k,2)=cos(hold)
      tracer(k,ndim+1)=k      ! "alf" particle label
   enddo
   tracer(numt_insert,2)=0    ! put last point on boundary

   ! add some points in the non-insert region:
   delalf = pi/(2*(n-(numt_insert+1)))
   do k=numt_insert+1,n
      ! hold=(k-(numt_insert+1))*delalf
      tracer(k,1)=xlocation
      tracer(k,2)=vxline_y(k-numt_insert)
      tracer(k,ndim+1)=k   ! "alf" particle lable 
   enddo
endif
end subroutine



subroutine comp_ellipse_reshape(w,setmax,center_only)
use params
use ellipse
implicit none
real*8 :: w(nx,ny)
integer :: setmax,center_only
call comp_ellipse(w,setmax,center_only)
end
