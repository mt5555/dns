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


call print_message("Restarting from file restart.[uvw]")
fname = rundir(1:len_trim(rundir)) // "restart.u"
call singlefile_io(time_initial,Q(1,1,1,1),fname,work1,work2,1,io_pe)
fname = rundir(1:len_trim(rundir)) // "restart.v"
call singlefile_io(time_initial,Q(1,1,1,2),fname,work1,work2,1,io_pe)
fname = rundir(1:len_trim(rundir)) // "restart.w"
call singlefile_io(time_initial,Q(1,1,1,3),fname,work1,work2,1,io_pe)

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
CPOINTER :: null=0;
character(len=80) message

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
   call input1(PSI(1,1,1,n),Q,work,null,io_pe)  
enddo



alpha=0
beta=1
do n=1,3
   call poisson(PSI(1,1,1,n),work,alpha,beta)
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






subroutine init_data_kh(Q)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)

! local variables
integer i,j,k,l
real*8 :: eps

Q=0

k=nz1
eps=200

do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   if (ycord(j)<=.5) then
      Q(i,j,k,1)=tanh(eps*(ycord(j)-.25))
   else
      Q(i,j,k,1)=tanh(eps*(.75-ycord(j)))
   endif
   if (nslabz>16) then
      Q(i,j,k,2)=.05*sin(2*pi*xcord(i))*cos(2*pi*zcord(k))
   else
      Q(i,j,k,2)=.05*sin(2*pi*xcord(i))
   endif
enddo
enddo
enddo

#if 0

do k=nz1+1,nz2
do i=nx1,nx2
do j=ny1,ny2
   Q(i,j,k,1)=Q(i,j,nz1,1)	
   Q(i,j,k,2)=Q(i,j,nz1,2)	
enddo
enddo
enddo
#endif


   


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







subroutine init_data_projection(Q,d1)
use params
use fft_interface
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: d1(nx,ny,nz)

! local variables
integer i


call bc_preloop(Q)

! will also dealias if dealias=1
call divfree(Q,d1)


end subroutine





