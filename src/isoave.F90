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


integer,parameter :: ndir=73
integer,parameter :: ndelta_max=73
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

real*8,allocatable  :: dwork2(:,:)     
real*8,allocatable  :: dwork3(:,:,:)   

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

if (my_pe/=io_pe) return


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






!
! compute structure functions for many different directions
! 
!
! parallel x-y hyperslab:  averages over whole cube only
!       once this is debugged, put in flag so if running
!       on 1 cpu, code will do no transposes (so that the
!       subcube code will still work)
!
!
subroutine isoavep(Q,Qs,Qt,Qst,range)
use params
use transpose

!input
real*8 :: Q(nx,ny,nz,ndim)               ! original data
real*8 :: Qt(g_nz2,nslabx,ny_2dz,ndim)   ! transpose
real*8 :: Qs(nx,ny,nz,ndim)              ! shifted original data
real*8 :: Qst(g_nz2,nslabx,ny_2dz,ndim)  ! transpose
real*8 :: range(3,2)

!local
real*8 :: rhat(3),rvec(3),rperp1(3),rperp2(3),delu(3),dir_shift(3)
real*8 :: u_l,u_t1,u_t2,rnorm
real*8 :: eta,lambda,r_lambda,ke_diss
real*8 :: dummy,xtmp
integer :: idir,idel,i2,j2,k2,i,j,k,n,m,ntot,ishift,k_g,j_g
integer :: n1,n1d,n2,n2d,n3,n3d,ierr

if (firstcall) then
   firstcall=.false.
   if (ncpu_x*ncpu_y>1) then
      call abort("isoave: requires x-y hyperslab parallel decomposition")
   endif   
   call init
endif



ntot=g_nx*g_ny*g_nz

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
ntot=0
do n=1,ndim
   do m=1,ndim
      call der(Q(1,1,1,n),Qs(1,1,1,1),dummy,Qs(1,1,1,2),1,m)
      do k=nz1,nz2
      if (zcord(k)>=range(3,1) .and. zcord(k)<range(3,2)) then
      do j=ny1,ny2
      if (ycord(j)>=range(2,1) .and. ycord(j)<range(2,2)) then
      do i=nx1,nx2
      if (xcord(i)>=range(1,1) .and. xcord(i)<range(1,2)) then

         ntot=ntot+1
         if (m==1) ke = ke + .5*Q(i,j,k,n)**2
         ke_diss=ke_diss + Qs(i,j,k,1)*Qs(i,j,k,1)

      endif
      enddo
      endif
      enddo
      endif
      enddo
   enddo
enddo

epsilon=mu*ke_diss/ntot
ke=ke/ntot


#ifdef USE_MPI
   xtmp=epsilon
   call MPI_allreduce(xtmp,epsilon,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   xtmp=ke
   call MPI_allreduce(xtmp,ke,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif



eta = (mu**3 / epsilon)**.25
lambda=sqrt(10*ke*mu/epsilon)       ! single direction lambda
R_lambda = lambda*sqrt(2*ke/3)/mu 

if (my_pe==io_pe) then
   print *,'KE:      ',ke
   print *,'epsilon: ',epsilon
   print *,'mu       ',mu
   print *
   print *,'eta      ',eta
   print *,'delx/eta ',1/(g_nmin*eta)
   print *,'lambda   ',lambda
   print *,'R_l      ',R_lambda
endif


do n=1,3
   call transpose_to_z(Q(1,1,1,n),Qt(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo      



do idir=1,ndir

   if (my_pe==io_pe) then
      write(*,'(a,i3,a,i3,a,3i3,a)') 'direction: ',idir,'/',ndir,'  (',dir(:,idir),')'
   endif

      rhat = dir(:,idir)
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

      ! data is in x-y hyperslabs.  
      ! z-direction=0
      !    then compute in x-y hyperslab without any transpose.
      ! z-direction/=0 and y-direction=0
      !    tranpose to x-z plane, compute there
      ! z-direction/=0 and y-direction/=0
      !    shift so y direction becomes 0
      !    tranpose to x-z plane, compute there
      !
      dir_shift=dir(:,idir)

!     if (dir_shift(3)==0 .or. ncpu_z==1) then
      if (dir_shift(3)==0) then  
         ! no shifts - compute directly 
         call comp_str_xy(Q,idir,rhat,rperp1,rperp2,dir_shift,range)
      else if (dir_shift(2)==0) then
         ! no need to shift, y-direction alread 0
         call comp_str_xz(Qt,idir,rhat,rperp1,rperp2,dir_shift,range)
      else if (mod(dir_shift(2),dir_shift(3))==0) then
         ! 
         ! shift in y by (y/z)*z-index:
         ! 
         do k=nz1,nz2
         do j=ny1,ny2
         do i=nx1,nx2
            k_g = k-nz1+1 + nslabz*my_z  ! global z index
            ishift=(dir_shift(2)/dir_shift(3))*(k_g-1)
            j2=j+ishift
            do
               if (ny1<=j2 .and. j2<=ny2) exit
               if (j2<ny1) j2=j2+nslaby
               if (j2>ny2) j2=j2-nslaby
            enddo
            do n=1,3
               Qs(i,j,k,n)=Q(i,j2,k,n)
            enddo
         enddo
         enddo
         enddo
         do n=1,3
            call transpose_to_z(Qs(1,1,1,n),Qst(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
         enddo
         call comp_str_xz(Qst,idir,rhat,rperp1,rperp2,dir_shift,range)
      else if (mod(dir_shift(3),dir_shift(2))==0) then
         ! 
         ! shift in z by (z/y)*y-index
         !
         ! shift  Qst = Qt-shifted
         do j=1,ny_2dz
         do i=1,nslabx
         do k=1,g_nz
            ! j_g = -1 + ny1 + j + my_z*ny_2dz     ! orig y index = [ny1,ny2]
            ! j_g = j_g - ny1 +1                   ! convert so starts at 1
            j_g = j + my_z*ny_2dz    
            ishift=(dir_shift(3)/dir_shift(2))*(j_g-1)
            k2=k+ishift
            do
               if (1<=k2 .and. k2<=g_nz) exit
               if (k2<1) k2=k2+g_nz
               if (k2>g_nz) k2=k2-g_nz
            enddo
            do n=1,3
               Qst(k,i,j,n)=Qt(k2,i,j,n)
            enddo
         enddo
         enddo
         enddo
         do n=1,3
            call transpose_from_z(Qst(1,1,1,n),Qs(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
         enddo
         call comp_str_xy(Qs,idir,rhat,rperp1,rperp2,dir_shift,range)
      else
         call abort("parallel computation of direction not supported")
      endif
   enddo


D_ll=D_ll/ntot
D_tt=D_tt/ntot
D_ltt=D_ltt/ntot
D_lll=D_lll/ntot

SP_ltt=SP_ltt/ntot
SP_lll=SP_lll/ntot
SN_ltt=SN_ltt/ntot
SN_lll=SN_lll/ntot

#ifdef USE_MPI
   dwork2=D_ll
   call MPI_allreduce(dwork2,D_ll,ndelta*ndir,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   dwork3=D_tt
   call MPI_allreduce(dwork3,D_tt,ndelta*ndir*2,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   dwork3=D_ltt
   call MPI_allreduce(dwork3,D_ltt,ndelta*ndir*2,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   dwork2=D_lll
   call MPI_allreduce(dwork2,D_lll,ndelta*ndir,MPI_REAL8,MPI_SUM,comm_3d,ierr)


   dwork3=SP_ltt
   call MPI_allreduce(dwork3,SP_ltt,ndelta*ndir*2,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   dwork2=SP_lll
   call MPI_allreduce(dwork2,SP_lll,ndelta*ndir,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   dwork3=SN_ltt
   call MPI_allreduce(dwork3,SN_ltt,ndelta*ndir*2,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   dwork2=SN_lll
   call MPI_allreduce(dwork2,SN_lll,ndelta*ndir,MPI_REAL8,MPI_SUM,comm_3d,ierr)



#endif

end subroutine





!
! same as above, but serial version
!
subroutine isoave1(Q,d1,work,lx1,lx2,ly1,ly2,lz1,lz2)
use params

!input
real*8 :: Q(nx,ny,nz,ndim)
real*8 :: d1(nx,ny,nz)
real*8 :: work(nx,ny,nz)
integer :: lx1,lx2,lz1,lz2,ly1,ly2

!local
real*8 :: rhat(3),rvec(3),rperp1(3),rperp2(3),delu(3)
real*8 :: u_l,u_t1,u_t2,rnorm
real*8 :: eta,lambda,r_lambda,ke_diss
real*8 :: dummy
integer :: idir,idel,i2,j2,k2,i,j,k,n,m,ntot

if (firstcall) then
   if (ncpu_x*ncpu_y*ncpu_z>1) then
      call abort("isoave1: can only run with 1 MPI cpu")
   endif   
   firstcall=.false.
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
      do k=lz1,lz2
      do j=ly1,ly2
      do i=lx1,lx2
         if (m==1) ke = ke + .5*Q(i,j,k,n)**2
         ke_diss=ke_diss + d1(i,j,k)*d1(i,j,k)
      enddo
      enddo
      enddo
   enddo
enddo
ntot=(lz2-lz1+1)*(ly2-ly1+1)*(lx2-lx1+1)
epsilon=mu*ke_diss/ntot
ke=ke/ntot

eta = (mu**3 / epsilon)**.25
lambda=sqrt(10*ke*mu/epsilon)       ! single direction lambda
R_lambda = lambda*sqrt(2*ke/3)/mu 


print *,'KE:      ',ke
print *,'epsilon: ',epsilon
print *,'mu       ',mu
print *
print *,'eta      ',eta
print *,'delx/eta ',1/(g_nmin*eta)
print *,'lambda   ',lambda
print *,'R_l      ',R_lambda



      



!$omp parallel do private(rhat,rperp1,rperp2,rvec,i2,j2,k2,n,delu,u_l,u_t1,u_t2)
do idir=1,ndir

   write(*,'(a,i3,a,i3,a,3i3,a)') 'direction: ',idir,'/',ndir,'  (',dir(:,idir),')'

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

!$XXX parallel do private(rvec,i2,j2,k2,n,delu,u_l,u_t1,u_t2)
   do idel=1,ndelta

      rvec = dir(:,idir)*delta_val(idel)
      ! dont bother computing deltas above 50% domain size
      if ( (rvec(1)**2+rvec(2)**2+rvec(3)**2) < g_nmin**2/4) then
      
      if (rvec(1)<0) rvec(1)=rvec(1)+nslabx
      if (rvec(2)<0) rvec(2)=rvec(2)+nslaby
      if (rvec(3)<0) rvec(3)=rvec(3)+nslabz
      
      
      do k=lz1,lz2
      do j=ly1,ly2
      do i=lx1,lx2

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
enddo
!$omp end parallel do


D_ll=D_ll/ntot
D_tt=D_tt/ntot
D_ltt=D_ltt/ntot
D_lll=D_lll/ntot

SP_ltt=SP_ltt/ntot
SP_lll=SP_lll/ntot
SN_ltt=SN_ltt/ntot
SN_lll=SN_lll/ntot


end subroutine




subroutine comp_str_xy(Q,idir,rhat,rperp1,rperp2,dir_base,range)
use params
implicit none
!input
real*8 :: Q(nx,ny,nz,ndim)       
real*8 :: rhat(3),rperp1(3),rperp2(3),dir_base(3)
real*8 :: range(3,2)

!local
real*8 :: rvec(3),delu(3)
real*8 :: u_l,u_t1,u_t2
integer :: idir,idel,i2,j2,k2,i,j,k,n,m


  do idel=1,ndelta

      rvec = dir_base*delta_val(idel)
      ! dont bother computing deltas above 50% domain size
      if ( (rvec(1)**2+rvec(2)**2+rvec(3)**2) < g_nmin**2/4) then
      
      if (rvec(1)<0) rvec(1)=rvec(1)+nslabx
      if (rvec(2)<0) rvec(2)=rvec(2)+nslaby
      rvec(3)=0  ! data required to be in xy slab
      
      
      do k=nz1,nz2
      if (zcord(k)>=range(3,1) .and. zcord(k)<range(3,2)) then
      do j=ny1,ny2
      if (ycord(j)>=range(2,1) .and. ycord(j)<range(2,2)) then
      do i=nx1,nx2
      if (xcord(i)>=range(1,1) .and. xcord(i)<range(1,2)) then

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

         do n=1,3
            delu(n)=Q(i2,j2,k,n)-Q(i,j,k,n)
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
         
      endif
      enddo
      endif
      enddo
      endif
      enddo
      endif
   enddo

end subroutine









subroutine comp_str_xz(Q,idir,rhat,rperp1,rperp2,dir_base,range)
use params
implicit none
!input
real*8 :: Q(g_nz2,nslabx,ny_2dz,ndim)  
real*8 :: rhat(3),rperp1(3),rperp2(3),dir_base(3)
real*8 :: range(3,2)

!local
real*8 :: rvec(3),delu(3)
real*8 :: u_l,u_t1,u_t2
integer :: idir,idel,i2,j2,k2,i,j,k,n,m


  do idel=1,ndelta

      rvec = dir_base*delta_val(idel)
      ! dont bother computing deltas above 50% domain size
      if ( (rvec(1)**2+rvec(2)**2+rvec(3)**2) < g_nmin**2/4) then
      
      if (rvec(1)<0) rvec(1)=rvec(1)+nslabx
      rvec(2)=0  ! data required to be in xz slab
      if (rvec(3)<0) rvec(3)=rvec(3)+g_nz
      
      do j=1,ny_2dz
      if (z_jmcord(j)>=range(2,1) .and. z_jmcord(j)<range(2,2)) then
      do i=1,nslabx
      if (z_imcord(i)>=range(1,1) .and. z_imcord(i)<range(1,2)) then
      do k=1,g_nz
      if (z_kmcord(k)>=range(3,1) .and. z_kmcord(k)<range(3,2)) then


         i2 = i + rvec(1)
         do
            if (i2<1) call abort("isoave: i2<nx1")
            if (i2<=nslabx) exit
            i2=i2-nslabx
         enddo

         k2 = k + rvec(3)
         do
            if (k2<1) call abort("isoave: k2<nz1")
            if (k2<=g_nz) exit
            k2=k2-g_nz
         enddo
         do n=1,3
            delu(n)=Q(k2,i2,j,n)-Q(k,i,j,n)
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
         
      endif
      enddo
      endif
      enddo
      endif
      enddo
      endif
   enddo

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



delta_val=100000

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
if (j > ndelta_max) then
   stop "isoave init: j > ndelta_max"
endif
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


! face 1,3 directions
dir(:,50)=(/0,1,3/)
dir(:,51)=(/0,3,1/)
dir(:,52)=(/0,-1,3/)
dir(:,53)=(/0,-3,1/)

dir(:,54)=(/1,3,0/)
dir(:,55)=(/1,0,3/)
dir(:,56)=(/-1,3,0/)
dir(:,57)=(/-1,0,3/)

dir(:,58)=(/3,1,0/)
dir(:,59)=(/3,0,1/)
dir(:,60)=(/-3,1,0/)
dir(:,61)=(/-3,0,1/)

! body 1,1,3 directions
dir(:,62)=(/1,1,3/)
dir(:,63)=(/1,3,1/)
dir(:,64)=(/3,1,1/)

dir(:,65)=(/-1,1,3/)
dir(:,66)=(/-1,3,1/)
dir(:,67)=(/-3,1,1/)

dir(:,68)=(/1,-1,3/)
dir(:,69)=(/1,-3,1/)
dir(:,70)=(/3,-1,1/)

dir(:,71)=(/1,1,-3/)
dir(:,72)=(/1,3,-1/)
dir(:,73)=(/3,1,-1/)




if (ndir/=73) then
   call abort("isoave: ndir /= 49")
endif

allocate(dwork2(ndelta,ndir))
allocate(dwork3(ndelta,ndir,2))

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
