#include "macros.h"
!
!  compute structure functions for a range of deltas, and
!  for a range of directions.  
!  
!  isoave1:   computed structure functions for the most directions,
!             but must be run in parallel
!
!
!
!

module isoave
implicit none

integer,parameter :: ndir=49
integer,parameter :: ndelta_max=72
integer :: dir(3,ndir)
integer :: ndelta
integer :: delta_val(ndelta_max)   ! delta values (in terms of grid points)
real*8  :: r_val(ndelta_max,ndir)      ! actual distance for each delta value


logical :: firstcall=.true.
!
!   r = r_val(idel,idir)
!   rhat = dir(:,idir) /  || dir(:,idir) ||
!   rperp1,2  two directions orthongonal to rhat
!
!   u_l  = (u(x+r)-u(x))  dot rhat
!   u_t1 = (u(x+r)-u(x))  dot rperp1
!   u_t2 = (u(x+r)-u(x))  dot rperp2
!
!   D_ll  = u_l**2
!   D_lll = u_l**3
!
!   D_tt  = u_t1**2
!   D_ltt = u_l * u_t1**2
!
!
!
real*8,allocatable  :: D_ll(:,:)         ! D_ll(ndelta,ndir)
real*8,allocatable  :: D_lll(:,:)        ! D_lll(ndelta,ndir)
real*8,allocatable  :: D_tt(:,:,:)       ! D_tt(ndelta,ndir,2)    
real*8,allocatable  :: D_ltt(:,:,:)      ! D_ltt(ndelta,ndir,2)    

! signed (positive part) of above
real*8,allocatable  :: SP_lll(:,:)        ! D_lll(ndelta,ndir)
real*8,allocatable  :: SP_ltt(:,:,:)      ! D_ltt(ndelta,ndir,2)    

! signed (positive part) of above
real*8,allocatable  :: SN_lll(:,:)        ! D_lll(ndelta,ndir)
real*8,allocatable  :: SN_ltt(:,:,:)      ! D_ltt(ndelta,ndir,2)    

! also added to the file for completeness:
real*8,private :: epsilon,mu,ke

private init

contains


subroutine writeisoave(fid,time)
use params
implicit none

CPOINTER fid
real*8 :: time
!local
integer :: i,idir
real*8 :: x


x=ndelta; call cwrite8(fid,x,1)   
x=ndir;   call cwrite8(fid,x,1)   
x=4;      call cwrite8(fid,x,1)   ! number of longitudinal (1 per direction)
x=4;      call cwrite8(fid,x,1)   ! number of transverse (2 per direction)
x=7;      call cwrite8(fid,x,1)   ! number of scalars
x=0;      call cwrite8(fid,x,1)   ! number of future type2


! write out the r values
do idir=1,ndir
   call cwrite8(fid,r_val(1,idir),ndelta)
enddo

! longitudinal
do idir=1,ndir
   call cwrite8(fid,D_ll(1,idir),ndelta)
enddo
do idir=1,ndir
   call cwrite8(fid,D_lll(1,idir),ndelta)
enddo
do idir=1,ndir
   call cwrite8(fid,SP_lll(1,idir),ndelta)
enddo
do idir=1,ndir
   call cwrite8(fid,SN_lll(1,idir),ndelta)
enddo

! transverse
do i=1,2
do idir=1,ndir
   call cwrite8(fid,D_tt(1,idir,i),ndelta)
enddo
enddo
do i=1,2
do idir=1,ndir
   call cwrite8(fid,D_ltt(1,idir,i),ndelta)
enddo
enddo
do i=1,2
do idir=1,ndir
   call cwrite8(fid,SP_ltt(1,idir,i),ndelta)
enddo
enddo
do i=1,2
do idir=1,ndir
   call cwrite8(fid,SN_ltt(1,idir,i),ndelta)
enddo
enddo

call cwrite8(fid,time,1)
x=g_nx; call cwrite8(fid,x,1)
x=g_ny; call cwrite8(fid,x,1)
x=g_nz; call cwrite8(fid,x,1)
call cwrite8(fid,mu,1)
call cwrite8(fid,ke,1)
call cwrite8(fid,epsilon,1)


end subroutine







subroutine isoave1(Q,d1,work)
use params

!input
real*8 :: Q(nx,ny,nz,ndim)
real*8 :: d1(nx,ny,nz)
real*8 :: work(nx,ny,nz)

!local
real*8 :: rhat(3),rvec(3),rperp1(3),rperp2(3),delu(3)
real*8 :: u_l,u_t1,u_t2,rnorm
real*8 :: eta,lambda,r_lambda,ke_diss
real*8 :: dummy
integer :: idir,idel,i2,j2,k2,i,j,k,n,m

if (firstcall) then
   firstcall=.false.
   if (ncpu_x*ncpu_y*ncpu_z>1) then
      call abort("isoave: can only be called if running on 1 cpu")
   endif
   call init
endif


D_ll=0
D_tt=0
D_ltt=0
D_lll=0
SP_ltt=0
SP_lll=0
SN_ltt=0
SN_lll=0


ke_diss=0
ke=0
do n=1,ndim
   do m=1,ndim
      call der(Q(1,1,1,n),d1,dummy,work,1,m)
      do k=nz1,nz2
      do j=ny1,ny2
      do i=nx1,nx2
         if (m==1) ke = ke + .5*Q(i,j,k,n)**2
         ke_diss=ke_diss + d1(i,j,k)*d1(i,j,k)
      enddo
      enddo
      enddo
   enddo
enddo
epsilon=mu*ke_diss/g_nx/g_ny/g_nz   
ke=ke/g_nx/g_ny/g_nz

eta = (mu**3 / epsilon)**.25
lambda=10*ke*mu/epsilon       ! single direction lambda
R_lambda = lambda*sqrt(2*ke/3)/mu 


print *,'KE:      ',ke
print *,'epsilon: ',epsilon
print *,'mu       ',mu
print *
print *,'eta      ',eta
print *,'delx/eta ',1/(g_nmin*eta)
print *,'lambda   ',lambda
print *,'R_l      ',R_lambda



      



!$XXX parallel do private(rhat,rperp1,rperp2,rvec,i2,j2,k2,n,delu,u_l,u_t1,u_t2)
do idir=1,ndir

   write(*,'(a,i3,a,i3,a,3i3,a)') 'direction: ',idir,'/',ndir,'  (',&
           dir(:,idir),')'

      rhat = dir(:,idir)*delta_val(1)
      rhat=rhat/sqrt(rhat(1)**2+rhat(2)**2+rhat(3)**2)
      call compute_perp(rhat,rperp1,rperp2)

#if 0
      ! check orthoginality
      print *,'norms: ',sqrt(rhat(1)**2+rhat(2)**2+rhat(3)**2), &
           sqrt(rperp1(1)**2+rperp1(2)**2+rperp1(3)**2), &
           sqrt(rperp2(1)**2+rperp2(2)**2+rperp2(3)**2)
      
      print *,'ortho: ',&
           rhat(1)*rperp1(1)+rhat(2)*rperp1(2)+rhat(3)*rperp1(3), &
           rhat(1)*rperp2(1)+rhat(2)*rperp2(2)+rhat(3)*rperp2(3), &
           rperp2(1)*rperp1(1)+rperp2(2)*rperp1(2)+rperp2(3)*rperp1(3)
      
#endif

!$omp parallel do private(rvec,i2,j2,k2,n,delu,u_l,u_t1,u_t2)
   do idel=1,ndelta

      rvec = dir(:,idir)*delta_val(idel)
      ! dont bother computing deltas above 50% domain size
      if ( (rvec(1)**2+rvec(2)**2+rvec(3)**2) < g_nmin**2/4) then
      
      if (rvec(1)<0) rvec(1)=rvec(1)+nslabx
      if (rvec(2)<0) rvec(2)=rvec(2)+nslaby
      if (rvec(3)<0) rvec(3)=rvec(3)+nslabz
      
      
      do k=nz1,nz2
      do j=ny1,ny2
      do i=nx1,nx2

         i2 = i + rvec(1)
         do
            if (i2<nx1) call abort("isoave: i2<nx1")
            if (i2<=nx2) exit
            i2=i2-nslabx
         enddo
         j2 = j + rvec(2)
         do
            if (j2<ny1) call abort("isoave: j2<ny1")
            if (j2<=ny2) exit
            j2=j2-nslaby
         enddo
         k2 = k + rvec(3)
         do
            if (k2<nz1) call abort("isoave: k2<nz1")
            if (k2<=nz2) exit
            k2=k2-nslabz
         enddo
         do n=1,3
            delu(n)=Q(i2,j2,k2,n)-Q(i,j,k,n)
         enddo
         u_l  = delu(1)*rhat(1)+delu(2)*rhat(2)+delu(3)*rhat(3)
         u_t1 = delu(1)*rperp1(1)+delu(2)*rperp1(2)+delu(3)*rperp1(3)
         u_t2 = delu(1)*rperp2(1)+delu(2)*rperp2(2)+delu(3)*rperp2(3)
         
         
         D_ll(idel,idir)=D_ll(idel,idir) + u_l**2
         D_tt(idel,idir,1)=D_tt(idel,idir,1) + u_t1**2
         D_tt(idel,idir,2)=D_tt(idel,idir,2) + u_t2**2
         
         D_lll(idel,idir)=D_lll(idel,idir) + u_l**3
         D_ltt(idel,idir,1)=D_ltt(idel,idir,1) + u_l*u_t1**2
         D_ltt(idel,idir,2)=D_ltt(idel,idir,2) + u_l*u_t2**2

         if (u_l>=0) then
            SP_lll(idel,idir)  =SP_lll(idel,idir) + u_l**3
            SP_ltt(idel,idir,1)=SP_ltt(idel,idir,1) + u_l*u_t1**2
            SP_ltt(idel,idir,2)=SP_ltt(idel,idir,2) + u_l*u_t2**2
         else
            SN_lll(idel,idir)  =SN_lll(idel,idir) - u_l**3
            SN_ltt(idel,idir,1)=SN_ltt(idel,idir,1) - u_l*u_t1**2
            SN_ltt(idel,idir,2)=SN_ltt(idel,idir,2) - u_l*u_t2**2
         endif
         
      enddo
      enddo
      enddo
      endif
enddo
!$omp end parallel do
enddo


D_ll=D_ll/g_nx/g_ny/g_nz
D_tt=D_tt/g_nx/g_ny/g_nz
D_ltt=D_ltt/g_nx/g_ny/g_nz
D_lll=D_lll/g_nx/g_ny/g_nz

SP_ltt=SP_ltt/g_nx/g_ny/g_nz
SP_lll=SP_lll/g_nx/g_ny/g_nz
SN_ltt=SN_ltt/g_nx/g_ny/g_nz
SN_lll=SN_lll/g_nx/g_ny/g_nz


end subroutine



!
! For each direction, find its intersection with the unit sphere
! (2 points) and output in lat-lon coordinates
!
!
subroutine writepoints()
use params

!input
real*8 :: Q(nx,ny,nz,ndim)

!local
real*8 :: rhat(3),rvec(3),rperp1(3),rperp2(3),delu(3)
real*8 :: lat,lon
integer :: idir,idel,i2,j2,k2,i,j,k,n

if (firstcall) then
   firstcall=.false.
   if (ncpu_x*ncpu_y*ncpu_z>1) then
      call abort("isoave: can only be called if running on 1 cpu")
   endif
   call init
endif

open(97,file="isoave.coords",form='formatted')
write(97,*) 2*ndir

do idir=1,ndir

      rhat = dir(:,idir)*delta_val(1)
      rhat=rhat/sqrt(rhat(1)**2+rhat(2)**2+rhat(3)**2)

      ! convert to lat-lon:
      !call latlon(rhat,lat,lon)
      write(97,'(3f25.20)') rhat
      !call latlon(rhat,lat,lon)
      write(97,'(3f25.20)') -rhat
enddo

close(97)
end subroutine



subroutine latlon(rhat,lat,lon)
use params  ! just to get pi
implicit none
real*8 rhat(3),lat,lon,cotheta

      if (abs(rhat(1)).lt.1e-9 .and. abs(rhat(2)).lt.1e-9) then 
         lon=0
      else
         lon=atan2(rhat(2),rhat(1))
      endif
      !if (lon.lt.0) lon=lon+2*pi
      cotheta=acos(rhat(3))
      lat=pi/2-cotheta
end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! initialize values of delta to use, and all directions
! to use
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init
use params
integer :: i,j,idel,idir
real*8 :: rvec(3)




! compute all deltas up to ndelta_max
! (but we only use deltas up to ndelta)
! 1..16
do i=1,16
   delta_val(i)=i
enddo
j=16

! 2..32  (start with 2*9)
do i=9,16
   j=j+1
   delta_val(j)=2*i
enddo

! 4..64  (start with 4*9)
do i=9,16
   j=j+1
   delta_val(j)=4*i
enddo

! 8..128 (start with 8*9)
do i=9,16
   j=j+1
   delta_val(j)=8*i
enddo

! 16..256  (start with 16*9)
do i=9,16
   j=j+1
   delta_val(j)=16*i
enddo

! 32..512  (start with 32*9)
do i=9,16
   j=j+1
   delta_val(j)=32*i
enddo

! 64..1024  (start with 64*9)
do i=9,16
   j=j+1
   delta_val(j)=64*i
enddo

! 128..2048
do i=9,16
   j=j+1
   delta_val(j)=128*i
enddo
! determine maximum value of delta to use for this grid
do idel=1,ndelta_max
   if (delta_val(idel) >= g_nmin/2) exit
   ndelta=idel
enddo





! x,y,z directions
dir(:,1)=(/1,0,0/)
dir(:,2)=(/0,1,0/)
dir(:,3)=(/0,0,1/)

! face diagonals
dir(:,4)=(/1,1,0/)
dir(:,5)=(/1,0,1/)
dir(:,6)=(/0,1,1/)
dir(:,7)=(/-1,1,0/)
dir(:,8)=(/-1,0,1/)
dir(:,9)=(/0,-1,1/)


! body diagonals
dir(:,10)=(/1,1,1/)
dir(:,11)=(/-1,1,1/)
dir(:,12)=(/1,-1,1/)
dir(:,13)=(/1,1,-1/)


! face 1,2 directions
dir(:,14)=(/0,1,2/)
dir(:,15)=(/0,2,1/)
dir(:,16)=(/0,-1,2/)
dir(:,17)=(/0,-2,1/)

dir(:,18)=(/1,2,0/)
dir(:,19)=(/1,0,2/)
dir(:,20)=(/-1,2,0/)
dir(:,21)=(/-1,0,2/)

dir(:,22)=(/2,1,0/)
dir(:,23)=(/2,0,1/)
dir(:,24)=(/-2,1,0/)
dir(:,25)=(/-2,0,1/)

! body 1,1,2 directions
dir(:,26)=(/1,1,2/)
dir(:,27)=(/1,2,1/)
dir(:,28)=(/2,1,1/)

dir(:,29)=(/-1,1,2/)
dir(:,30)=(/-1,2,1/)
dir(:,31)=(/-2,1,1/)

dir(:,32)=(/1,-1,2/)
dir(:,33)=(/1,-2,1/)
dir(:,34)=(/2,-1,1/)

dir(:,35)=(/1,1,-2/)
dir(:,36)=(/1,2,-1/)
dir(:,37)=(/2,1,-1/)

! 2,2,1 directions
dir(:,38)=(/1,2,2/)
dir(:,39)=(/2,1,2/)
dir(:,40)=(/2,2,1/)

dir(:,41)=(/-1,2,2/)
dir(:,42)=(/-2,1,2/)
dir(:,43)=(/-2,2,1/)

dir(:,44)=(/1,-2,2/)
dir(:,45)=(/2,-1,2/)
dir(:,46)=(/2,-2,1/)

dir(:,47)=(/1,2,-2/)
dir(:,48)=(/2,1,-2/)
dir(:,49)=(/2,2,-1/)
if (ndir/=49) then
   call abort("isoave: ndir /= 49")
endif


allocate(D_ll(ndelta,ndir))
allocate(D_lll(ndelta,ndir))
allocate(D_tt(ndelta,ndir,2))
allocate(D_ltt(ndelta,ndir,2))

allocate(SP_lll(ndelta,ndir))
allocate(SP_ltt(ndelta,ndir,2))

allocate(SN_lll(ndelta,ndir))
allocate(SN_ltt(ndelta,ndir,2))


do idel=1,ndelta
do idir=1,ndir
   rvec = dir(:,idir)*delta_val(idel)
   r_val(idel,idir) =(rvec(1)**2+rvec(2)**2+rvec(3)**2)
   r_val(idel,idir) = sqrt(r_val(idel,idir))
enddo
enddo


end subroutine



subroutine compute_perp(r,rp1,rp2)
real*8 r(3),rp1(3),rp2(3)

if (r(3)==0) then ! either r(1) or r(2) <> 0
   rp1(1)=-r(2)
   rp1(2)=r(1)
   rp1(3)=0
   rp1 = rp1/ sqrt(rp1(1)**2+rp1(2)**2)
else 
   rp1(1)=0
   rp1(2)=-r(3)
   rp1(3)=r(2)
   rp1 = rp1/ sqrt(rp1(2)**2+rp1(3)**2)
endif


! take cross product for remaining vector:
! rp2 = r cross rp1
rp2(1) =  r(2)*rp1(3) - r(3)*rp1(2)
rp2(2) =-(r(1)*rp1(3) - r(3)*rp1(1))
rp2(3) =  r(1)*rp1(2) - r(2)*rp1(1)


end subroutine


end module 
