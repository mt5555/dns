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

integer,parameter :: ndir=37
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



private init

contains


subroutine writeisoave(fid)
implicit none

CPOINTER fid
!local
integer :: i,idir
real*8 :: x


x=ndelta; call cwrite8(fid,x,1)   
x=ndir;   call cwrite8(fid,x,1)   
x=2;      call cwrite8(fid,x,1)   ! number of longitudinal (1 per direction)
x=2;      call cwrite8(fid,x,1)   ! number of transverse (2 per direction)
x=0;      call cwrite8(fid,x,1)   ! number of future type1
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
end subroutine







subroutine isoave1(Q)
use params

!input
real*8 :: Q(nx,ny,nz,ndim)

!local
real*8 :: rhat(3),rvec(3),rperp1(3),rperp2(3),delu(3)
real*8 :: u_l,u_t1,u_t2
integer :: idir,idel,i2,j2,k2,i,j,k,n

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


do idir=1,ndir

   write(*,'(a,i3,a,i3,a,3i3,a)') 'direction: ',idir,'/',ndir,'  (',&
           dir(:,idir),')'

   do idel=1,ndelta
      
      rvec = dir(:,idir)*delta_val(idel)
      rhat = rvec/sqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)
      call compute_perp(rhat,rperp1,rperp2)
#if 0
      ! check orthoginality
      print *,'idel,idir',idel,idir
      print *,'norms: ',sqrt(rhat(1)**2+rhat(2)**2+rhat(3)**2), &
           sqrt(rperp1(1)**2+rperp1(2)**2+rperp1(3)**2), &
           sqrt(rperp2(1)**2+rperp2(2)**2+rperp2(3)**2)
      
      print *,'ortho: ',&
           rhat(1)*rperp1(1)+rhat(2)*rperp1(2)+rhat(3)*rperp1(3), &
           rhat(1)*rperp2(1)+rhat(2)*rperp2(2)+rhat(3)*rperp2(3), &
           rperp2(1)*rperp1(1)+rperp2(2)*rperp1(2)+rperp2(3)*rperp1(3)
      
#endif
      
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
         
         D_lll(idel,idir)=D_ll(idel,idir) + u_l**3
         D_ltt(idel,idir,1)=D_tt(idel,idir,1) + u_l*u_t1**2
         D_ltt(idel,idir,2)=D_tt(idel,idir,2) + u_l*u_t2**2
         
      enddo
      enddo
      enddo
enddo
enddo

D_ll=D_ll/g_nx/g_ny/g_nz
D_tt=D_tt/g_nx/g_ny/g_nz
D_ltt=D_ltt/g_nx/g_ny/g_nz
D_lll=D_lll/g_nx/g_ny/g_nz


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

if (ndir/=37) then
   call abort("isoave: ndir /= 37")
endif
allocate(D_ll(ndelta,ndir))
allocate(D_lll(ndelta,ndir))
allocate(D_tt(ndelta,ndir,2))
allocate(D_ltt(ndelta,ndir,2))


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
rp2(1) =  r(2)*rp1(3) - r(3)*rp2(2)
rp2(2) =-(r(1)*rp1(3) - r(3)*rp2(1))
rp2(3) =  r(1)*rp1(2) - r(2)*rp2(1)


end subroutine


end module 
