#include "macros.h"
subroutine init_data_restart(Q,work1,work2)
!
! low wave number, quasi isotropic initial condition
!
use params
use mpi
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

!local
character(len=80) message
character(len=80) fname

Q=0
if (equations==NS_UVW) then
   call print_message("Restarting from file restart.[uvw]")
   fname = rundir(1:len_trim(rundir)) // "restart.u"
   call singlefile_io(time_initial,Q(1,1,1,1),fname,work1,work2,1,io_pe)
   fname = rundir(1:len_trim(rundir)) // "restart.v"
   call singlefile_io(time_initial,Q(1,1,1,2),fname,work1,work2,1,io_pe)
   fname = rundir(1:len_trim(rundir)) // "restart.w"
   call singlefile_io(time_initial,Q(1,1,1,3),fname,work1,work2,1,io_pe)
else if (equations==SHALLOW) then
   call print_message("Restarting from file restart.[uvh]")
   fname = rundir(1:len_trim(rundir)) // "restart.u"
   call singlefile_io(time_initial,Q(1,1,1,1),fname,work1,work2,1,io_pe)
   fname = rundir(1:len_trim(rundir)) // "restart.v"
   call singlefile_io(time_initial,Q(1,1,1,2),fname,work1,work2,1,io_pe)
   fname = rundir(1:len_trim(rundir)) // "restart.h"
   call singlefile_io(time_initial,Q(1,1,1,3),fname,work1,work2,1,io_pe)
else if (equations==NS_PSIVOR) then
   call print_message("Restarting from file restart.vor")
   fname = rundir(1:len_trim(rundir)) // "restart.vor"
   call singlefile_io(time_initial,Q(1,1,1,1),fname,work1,work2,1,io_pe)
endif



write(message,'(a,f10.4)') "restart time=",time_initial
call print_message(message)
end subroutine




subroutine init_data_lwisotropic(Q,PSI,work,work2)
!
! low wave number, quasi isotropic initial condition
!
use params
use mpi
use transpose
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: PSI(nx,ny,nz,n_var)
real*8 :: work(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

! local variables
real*8 :: alpha,beta
integer km,jm,im,i,j,k,n,wn,ierr
integer,allocatable :: seed(:)
real*8 xw,ener1,ener2,ener1_target,ener2_target,ener,xfac
character(len=80) message
CPOINTER :: null

!
! 
!   g_xcord(i)=(i-1)*delx	

!random vorticity

! set the seed - otherwise it will be the same for all CPUs,
! producing a bad initial condition
call random_seed(size=k)
allocate(seed(k))
call random_seed(get=seed)
seed=seed+my_pe
call random_seed(put=seed)
deallocate(seed)

do n=1,3
   ! input from random number generator
   ! this gives same I.C independent of cpus
   call input1(PSI(1,1,1,n),work2,work,null,io_pe,.true.)  
enddo



alpha=0
beta=1
do n=1,3
   call helmholtz_periodic_inv(PSI(1,1,1,n),work,alpha,beta)
enddo
call vorticity(Q,PSI,work,work2)


ener1=0
ener2=0
ener1_target=1.0**(-5.0/3.0)
ener2_target=2.0**(-5.0/3.0)

do n=1,3
   call fft3d(Q(1,1,1,n),work) 
   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)
            xw=sqrt(real(km**2+jm**2+im**2))

            xfac = (2*2*2)
            if (km==0) xfac=xfac/2
            if (jm==0) xfac=xfac/2
            if (im==0) xfac=xfac/2
            if (xw>=.5 .and. xw<1.5) then
               ener1=ener1+.5*xfac*Q(i,j,k,n)**2
            else if (xw>=1.5 .and. xw<2.5) then
               ener2=ener2+.5*xfac*Q(i,j,k,n)**2
            else
               Q(i,j,k,n)=0
            endif
         enddo
      enddo
   enddo
enddo

#ifdef USE_MPI
   ener=ener1
   call MPI_allreduce(ener,ener1,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   ener=ener2
   call MPI_allreduce(ener,ener2,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif

if (ener1<=0) then
   call abort("ener1 must be positive")
endif
if (ener2<=0) then
   call abort("ener2 must be positive")
endif


do n=1,3
   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)
            xw=sqrt(real(km**2+jm**2+im**2))

            if (xw>=.5 .and. xw<1.5) then
               Q(i,j,k,n)=Q(i,j,k,n)*sqrt(ener1_target/(ener1))
            else if (xw>=1.5 .and. xw<2.5) then
               Q(i,j,k,n)=Q(i,j,k,n)*sqrt(ener2_target/(ener2))
            endif
         enddo
      enddo
   enddo
enddo

ener=0
ener1=0
ener2=0
do n=1,3
   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)
            xw=sqrt(real(km**2+jm**2+im**2))

            xfac = (2*2*2)
            if (km==0) xfac=xfac/2
            if (jm==0) xfac=xfac/2
            if (im==0) xfac=xfac/2
            if (xw>=.5 .and. xw<1.5) then
               ener1=ener1+.5*xfac*Q(i,j,k,n)**2
            else if (xw>=1.5 .and. xw<2.5) then
               ener2=ener2+.5*xfac*Q(i,j,k,n)**2
            endif
            ener=ener+.5*xfac*Q(i,j,k,n)**2


         enddo
      enddo
   enddo
enddo
do n=1,3
   call ifft3d(Q(1,1,1,n),work) 
enddo
#ifdef USE_MPI
   xfac=ener
   call MPI_allreduce(xfac,ener,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   xfac=ener1
   call MPI_allreduce(xfac,ener1,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   xfac=ener2
   call MPI_allreduce(xfac,ener2,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif

call print_message("Isotropic initial condition in wave numbers:");
write(message,'(a,f8.4,a,f8.4)') "0.5 <= wn < 1.5   E=",ener1,"  E target=",ener1_target
call print_message(message)
write(message,'(a,f8.4,a,f8.4)') "1.5 <= wn < 2.5   E=",ener2,"  E target=",ener2_target
call print_message(message)
write(message,'(a,f8.4,a,f8.4)') "Total E=",ener
call print_message(message)




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
integer i,j,k,l
real*8 delta,delsq,delalf,delgam,yval,xval,dify,difx,uu,vv,denom
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






subroutine init_data_vxpair(Q,q1,work1,w)
use params
use ghost
use bc
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: w(nx,ny,nz)

! local variables
integer :: i,j,k,l,use3d=0,n
real*8 :: eps
real*8 :: amp,vx,uy
real*8 :: delsq,hold,difx,dify,sumy,denom1,denom2,wsum,psisum,delalf,testmax
real*8 :: delta,xlocation
integer,parameter :: nd=200
real*8 :: xd(0:nd),yd(0:nd),wd(0:nd)

offset_bdy=1

if (init_cond_subtype ==0) then
!  my standard test case:
   biotsavart_cutoff=.001  
   delta=.2
   ubar=.089
   yscale=3

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
   biotsavart_cutoff=0
   if (g_bdy_x1==INFLOW0_ONESIDED) biotsavart_cutoff=.005
   delta=.2
   ubar=.089
   yscale=2.0
endif

! choose xscale so that delx*xscale = dely*yscale
xscale = yscale*(o_nx-1)/(o_ny-1)

call init_grid()   ! redo grid points since we changed scalings


xlocation=xscale/2
ubar = ubar  
equations=NS_PSIVOR
Q=0


! INITIALIZES VORTEX SHEET (XS,YS,WS,NS)
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


delsq = delta**2
!   INITIALIZES VORTICITY ON THE GRID = VORTICITY INDUCED
!     BY A BLOB AT (0,1)
do i=intx1,intx2
do j=inty1,inty2
   wsum = 0
   psisum=0
   do k=0,nd
      difx = (xcord(i) - xd(k))
      dify = (ycord(j) - yd(k))
      sumy = (ycord(j) + yd(k))

      denom1 = difx**2 + dify**2 + delsq
      denom2 = difx**2 + sumy**2 + delsq
      wsum = wsum + wd(k)*(1/denom1**2-1/denom2**2)
      !psisum = psisum - wd(k)*log(denom1/denom2)
   enddo
   w(i,j,1) = delsq*wsum/pi
   !Q(i,j,1,1)= (psisum  - ubar*ycord(j))/(4*pi)  scaling is wrong
enddo
enddo

call bcw_impose(w)
Q(:,:,:,3)=w








end subroutine







subroutine init_data_sht(Q,PSI,work1,work2)
use params
use transpose
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: PSI(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

! local variables
CPOINTER :: null
integer i,j,k,l
integer :: m
integer :: n,im,jm,km,ixw,ierr
real*8 :: k_0,xw,xfac ,alpha,beta,dummy
real*8 :: E_target,ke,pe,ke2
real*8 :: E_k(0:max(nx,ny,nz))
real*8 :: E_k2(0:max(nx,ny,nz))
real*8 :: Len,U,R,F
character(len=80) :: message

equations = SHALLOW   ! shallow water equations, not the default NS equations

k_0=14  
m=25
U=1.0        ! velocity scale  U normalized below so that KE=.5 U**2 = .5
Len=1.0/(2*pi*k_0)  ! length scale, determined by mode k_0 above



! set 
! R = Rossby number
! F = Froude number
if (init_cond_subtype==0) then        ! run A
   R=.01
   F=.04
else if (init_cond_subtype==1) then   ! run B
   R=.05
   F=.05
else if (init_cond_subtype==2) then   ! run C
   R=.05
   F=.075
else if (init_cond_subtype==3) then   ! run D
   R=.25
   F=.05
else if (init_cond_subtype==4) then   ! run E
   R=.25
   F=.20
else if (init_cond_subtype==5) then   ! run F
   R=1.00
   F=0.05
else if (init_cond_subtype==6) then   ! run G
   R=1.00
   F=0.30
else if (init_cond_subtype==7) then   ! run H
   R=5.00
   F=0.05
else if (init_cond_subtype==8) then   ! run I
   R=5.00
   F=0.30
else if (init_cond_subtype==9) then   ! run J
   R=20.00
   F=00.30
else if (init_cond_subtype==10) then   ! run K
   R=10.00
   F=00.10
else if (init_cond_subtype==11) then   ! run L
   R=20.00
   F=00.05
else if (init_cond_subtype==12) then   ! run M
   R=2.00
   F=.1
else if (init_cond_subtype==13) then   ! run N
   R=.4
   F=.1
else
   call abort("init_data_swt(): init_cond_subtype set to unsupported value")
endif
   

! set dimensionless constants used in code: f-plane and gravity:
H0 = 1.0
fcor=U/(Len*R)
grav=U**2/(H0*F**2)              

write(message,'(a,2f6.2)') "R, F parameters: ",R,F
call print_message(message)


Q=0
! random vorticity
call input1(PSI,work1,work2,null,io_pe,.true.)  

! invert to get PSI
alpha=0
beta=1
call helmholtz_periodic_inv(PSI,work1,alpha,beta)

! initialize amplitude to satisfy:  E(k)/.5k^2
! E(k) = k^(m/2) / (k + k_0)^m
!
!


n=1
call fft3d(PSI,work1) 
E_k=0
do k=nz1,nz2
   km=kmcord(k)
   do j=ny1,ny2
      jm=jmcord(j)
      do i=nx1,nx2
         im=imcord(i)
         xw=sqrt(real(km**2+jm**2+im**2))
         ixw = nint(xw)
         
         xfac = (2*2*2)
         if (km==0) xfac=xfac/2
         if (jm==0) xfac=xfac/2
         if (im==0) xfac=xfac/2

         E_k(ixw)=E_k(ixw)+.5*pi2_squared*xw*xw*xfac*PSI(i,j,k,n)**2
      enddo
   enddo
enddo

#ifdef USE_MPI
   E_k2=E_k
   call MPI_allreduce(E_k2,E_k,max(nx,ny,nz),MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif


do k=nz1,nz2
   km=kmcord(k)
   do j=ny1,ny2
      jm=jmcord(j)
      do i=nx1,nx2
         im=imcord(i)
         xw=sqrt(real(km**2+jm**2+im**2))
         ixw = nint(xw)
         
         xfac = (2*2*2)
         if (km==0) xfac=xfac/2
         if (jm==0) xfac=xfac/2
         if (im==0) xfac=xfac/2

         if (xw>0) then
            E_target=(xw**(m/2.0) / (xw+k_0)**m) 

            PSI(i,j,k,n)=PSI(i,j,k,n)/sqrt(E_k(ixw))
            PSI(i,j,k,n)=PSI(i,j,k,n)*sqrt(E_target)
         else
            PSI(i,j,k,n)=0
         endif
      enddo
   enddo
enddo


call ifft3d(PSI,work1) 




! compute velocity = curl PSI
Q=0
call der(PSI,Q,dummy,work1,DX_ONLY,2)
Q(:,:,:,1)=-Q(:,:,:,1)
call der(PSI,Q(1,1,1,2),dummy,work1,DX_ONLY,1)

if (dealias) call dealias_gridspace(Q,work1)


! normalize so <u,u>=U**2
ke = .5*sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,1)**2) +  &
     .5*sum(Q(nx1:nx2,ny1:ny2,nz1:nz2,2)**2)  

#ifdef USE_MPI
   ke2=ke
   call MPI_allreduce(ke2,ke,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
ke=ke/g_nx/g_ny
Q=sqrt(.5*U**2) * Q/sqrt(ke)



if (dealias) call dealias_gridspace(Q,work1)

! compute a balenced hight field:
! div( u grad u + f cross u  + g grad h ) = 0
! 
! PSI = u grad u + f cross u 
PSI=0
PSI(:,:,:,1)= fcor*Q(:,:,:,2)
PSI(:,:,:,2)=-fcor*Q(:,:,:,1)

! advection and viscous terms
k=1
do n=1,2

   ! compute u_x
   call der(Q(1,1,1,n),work1,dummy,work2,DX_ONLY,1)
   do j=ny1,ny2
   do i=nx1,nx2
      PSI(i,j,k,n) = PSI(i,j,k,n) - Q(i,j,k,1)*work1(i,j,k) 
   enddo
   enddo


   ! compute u_y
   call der(Q(1,1,1,n),work1,dummy,work2,DX_ONLY,2)
   do j=ny1,ny2
   do i=nx1,nx2
      PSI(i,j,k,n) = PSI(i,j,k,n) - Q(i,j,k,2)*work1(i,j,k) 
   enddo
   enddo

enddo

call divfree_gridspace(PSI,Q(1,1,1,3),work1,work2)
Q(:,:,:,3)= H0 + Q(:,:,:,3)/grav




end subroutine








