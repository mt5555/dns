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

! compute SK's helical structure function

module isoave
implicit none


integer,parameter :: ndir_max=73
integer,parameter :: ndelta_max=73
integer,parameter :: pmax=10       ! has to be 6 or greater
integer :: dir(3,ndir_max)
integer :: ndelta,ndir
integer :: delta_val(ndelta_max)   ! delta values (in terms of grid points)
real*8  :: r_val(ndelta_max,ndir_max)      ! actual distance for each delta value

integer :: ntranspose               ! number of data transpose operations performed

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
!   Dl  = u_l**p   p=2..10
!   Dt  = u_t**p   p=2..10
!
!   D_ltt = u_l * u_t1**2
!
!
!   
!
real*8  :: Dl_mean(2:pmax)
real*8,allocatable  :: Dl(:,:,:)         ! longitudinal, p=2..10
real*8,allocatable  :: Dt(:,:,:,:)       ! transverse, p=2..10
real*8,allocatable  :: D_ltt(:,:,:)      ! D_ltt(ndelta,ndir,2)    
real*8,allocatable  :: D_lltt(:,:)       ! 

real*8,allocatable  :: H_ltt(:,:)      ! H_ltt(ndelta,ndir)    
real*8,allocatable  :: H_tt(:,:)       ! 

! signed (positive part) of above
real*8,allocatable  :: SP_lll(:,:)        ! D_lll(ndelta,ndir)
real*8,allocatable  :: SP_ltt(:,:,:)      ! D_ltt(ndelta,ndir,2)    

! signed (positive part) of above
real*8,allocatable  :: SN_lll(:,:)        ! D_lll(ndelta,ndir)
real*8,allocatable  :: SN_ltt(:,:,:)      ! D_ltt(ndelta,ndir,2)    

! scalar structure functions:
!
!  w2s2(:,:,1)   <w2(x)*w2(x+r)>
!  w2s2(:,:,2)   <s2(x)*s2(x+r)>
!  w2s2(:,:,3)   <w2(x)*s2(x+r)>
!
real*8,allocatable :: w2s2(:,:,:)
real*8 :: w2s2_mean(3)


real*8,allocatable  :: dwork2(:,:)
real*8,allocatable  :: dwork3(:,:,:)




! also added to the file for completeness:
real*8,private :: epsilon,mu,ke

private init

contains




subroutine writeisoave_scalar(fid,time)
use params
implicit none

CPOINTER fid
real*8 :: time
!local
integer :: i,idir,p
real*8 :: x

if (my_pe/=io_pe) return

call cwrite8(fid,time,1)
x=ndelta; call cwrite8(fid,x,1)   
x=ndir;   call cwrite8(fid,x,1)   
x=pmax-2+1;      call cwrite8(fid,x,1)   ! number of str functions

call cwrite8(fid,Dl_mean,pmax-2+1)

! write out the r values
do idir=1,ndir
   call cwrite8(fid,r_val(1,idir),ndelta)
enddo

! longitudinal
do p=2,pmax
do idir=1,ndir
   call cwrite8(fid,Dl(1,idir,p),ndelta)
enddo
enddo
end subroutine




subroutine writeisoave_w2s2(fid,time)
use params
implicit none

CPOINTER fid
real*8 :: time
!local
integer :: i,idir,p
real*8 :: x

if (my_pe/=io_pe) return

call cwrite8(fid,time,1)
x=ndelta; call cwrite8(fid,x,1)   
x=ndir;   call cwrite8(fid,x,1)   
x=3;      call cwrite8(fid,x,1)   ! number of correlations

call cwrite8(fid,w2s2_mean,3)

! write out the r values
do idir=1,ndir
   call cwrite8(fid,r_val(1,idir),ndelta)
enddo

! longitudinal
do p=1,3
do idir=1,ndir
   call cwrite8(fid,w2s2(1,idir,p),ndelta)
enddo
enddo
end subroutine







subroutine writeisoave(fid,time)
use params
implicit none

CPOINTER fid
real*8 :: time
!local
integer :: i,idir,p
real*8 :: x

if (my_pe/=io_pe) return


x=ndelta; call cwrite8(fid,x,1)   
x=ndir;   call cwrite8(fid,x,1)   
x=pmax-1+5;   call cwrite8(fid,x,1)   ! number of longitudinal (1 per direction)
x=pmax-1+3;    call cwrite8(fid,x,1)   ! number of transverse (2 per direction)
x=7;      call cwrite8(fid,x,1)   ! number of scalars
x=0;      call cwrite8(fid,x,1)   ! number of future type2


! write out the r values
do idir=1,ndir
   call cwrite8(fid,r_val(1,idir),ndelta)
enddo

! longitudinal
do p=2,pmax
do idir=1,ndir
   call cwrite8(fid,Dl(1,idir,p),ndelta)
enddo
enddo

do idir=1,ndir
   call cwrite8(fid,SP_lll(1,idir),ndelta)
enddo
do idir=1,ndir
   call cwrite8(fid,SN_lll(1,idir),ndelta)
enddo
do idir=1,ndir
   call cwrite8(fid,H_ltt(1,idir),ndelta)
enddo
do idir=1,ndir
   call cwrite8(fid,H_tt(1,idir),ndelta)
enddo
do idir=1,ndir
   call cwrite8(fid,D_lltt(1,idir),ndelta)
enddo


! transverse
do p=2,pmax
do i=1,2
do idir=1,ndir
   call cwrite8(fid,Dt(1,idir,i,p),ndelta)
enddo
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
!  stype=3   assume Q is (u,v,w) data.  compute standard velocity
!                                    structure functions
!
!  stype=2   assume Q = (w,s), the vorticity and strain.  compute CJ's 
!            correlations, of the form <w(x)s(x+r>
!
!  stype=1   assume Q = scalar.  compute <(h(x+r)-h(x))^p>  p=2..6
!                   
!
subroutine isoavep(Q,Qs,Qt,Qst,stype,csig)
use params
use transpose

!input
integer :: stype
real*8 :: Q(nx,ny,nz,n_var)               ! original data
real*8 :: Qt(g_nz2,nslabx,ny_2dz,n_var)   ! transpose
real*8 :: Qs(nx,ny,nz,n_var)              ! shifted original data
real*8 :: Qst(g_nz2,nslabx,ny_2dz,n_var)  ! transpose


!local
real*8 :: rhat(3),rvec(3),rperp1(3),rperp2(3),delu(3),dir_shift(3)
real*8 :: u_l,u_t1,u_t2,rnorm
real*8 :: eta,lambda,r_lambda,ke_diss
real*8 :: dummy(3),xtmp,ntot
character(len=80) :: message
integer :: idir,idel,i2,j2,k2,i,j,k,n,m,ishift,k_g,j_g,nd
integer :: n1,n1d,n2,n2d,n3,n3d,ierr,p,csig

if (firstcall) then
   firstcall=.false.
   if (ncpu_x*ncpu_y>1) then
      call abort("isoave: requires x-y hyperslab parallel decomposition")
   endif   
   call init
endif

nd=3  ! number of variables stored in Q.
if (stype==3) nd=3
if (stype==2) nd=2
if (stype==1) nd=1

ntranspose=0
ntot=real(g_nx)*g_ny*g_nz


call zero_str

if (stype==3) then
ke_diss=0
ke=0
do n=1,ndim
   do m=1,ndim
      call der(Q(1,1,1,n),Qs(1,1,1,1),dummy,Qs(1,1,1,2),1,m)
      do k=nz1,nz2
         do j=ny1,ny2
            do i=nx1,nx2
               
               if (m==1) ke = ke + .5*Q(i,j,k,n)**2
               ke_diss=ke_diss + Qs(i,j,k,1)*Qs(i,j,k,1)
               
            enddo
         enddo
      enddo
   enddo
enddo


#ifdef USE_MPI
   xtmp=ke_diss
   call MPI_allreduce(xtmp,ke_diss,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   xtmp=ke
   call MPI_allreduce(xtmp,ke,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif


epsilon=mu*ke_diss/ntot
if (epsilon==0) epsilon=1e-20

ke=ke/ntot



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
endif

! compute the mean, r=0 values:
if (stype==2) then
   ! <w2,w2>
   w2s2_mean=0
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      w2s2_mean(1)=w2s2_mean(1)+Q(i,j,k,1)**2
      w2s2_mean(2)=w2s2_mean(2)+Q(i,j,k,2)**2
      w2s2_mean(3)=w2s2_mean(3)+Q(i,j,k,1)*Q(i,j,k,2)
   enddo
   enddo
   enddo
   w2s2_mean=w2s2_mean/g_nx/g_ny/g_nz
#ifdef USE_MPI
   dummy=w2s2_mean
   call MPI_allreduce(dummy,w2s2_mean,3,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif

endif

! compute the mean, r=0 values:
if (stype==1) then
   ! <w2,w2>
   Dl_mean=0
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      do p=2,pmax
         Dl_mean(p)=Dl_mean(p)+Q(i,j,k,1)**2
      enddo
   enddo
   enddo
   enddo
   Dl_mean=Dl_mean/g_nx/g_ny/g_nz
endif


do n=1,nd
   ntranspose=ntranspose+1
   call transpose_to_z(Q(1,1,1,n),Qt(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo   
   

do idir=1,ndir


!  check for SIGURG, telling us to stop because job will soon be killed
   call caught_sig(csig); 
#ifdef USE_MPI
   i=csig;
   call MPI_allreduce(i,csig,1,MPI_INTEGER,MPI_MAX,comm_3d,ierr)
#endif
   if (csig>0) exit;

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
         call comp_str_xy(Q,stype,idir,rhat,rperp1,rperp2,dir_shift)
      else if (dir_shift(2)==0) then
         ! no need to shift, y-direction already 0
         call comp_str_xz(Qt,stype,idir,rhat,rperp1,rperp2,dir_shift)
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
            do n=1,nd
               Qs(i,j,k,n)=Q(i,j2,k,n)
            enddo
         enddo
         enddo
         enddo
         do n=1,nd
            ntranspose=ntranspose+1
            call transpose_to_z(Qs(1,1,1,n),Qst(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
         enddo
         call comp_str_xz(Qst,stype,idir,rhat,rperp1,rperp2,dir_shift)
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
            do n=1,nd
               Qst(k,i,j,n)=Qt(k2,i,j,n)
            enddo
         enddo
         enddo
         enddo
         do n=1,nd
            ntranspose=ntranspose+1
            call transpose_from_z(Qst(1,1,1,n),Qs(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
         enddo
         call comp_str_xy(Qs,stype,idir,rhat,rperp1,rperp2,dir_shift)
      else
         call abort("parallel computation of direction not supported")
      endif
   enddo


   write(message,'(a,i5)') 'isoavep: number of calls to transpose_*() ',ntranspose
   call print_message(message)
   call normalize_and_reduce(ntot)
   call print_message('done with str normalize_and_reduce')

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
subroutine isoavep_subcube(Q,Qs,Qt,Qst,range,subcube,subcube_t,subcube_s,subcube_st)
use params
use transpose

!input
real*8 :: Q(nx,ny,nz,ndim)               ! original data
real*8 :: Qt(g_nz2,nslabx,ny_2dz,ndim)   ! transpose
real*8 :: Qs(nx,ny,nz,ndim)              ! shifted original data
real*8 :: Qst(g_nz2,nslabx,ny_2dz,ndim)  ! transpose
real*8 :: range(3,2)

! logical variable 1=included in subcube, 0=not included
! these have to be real*8 so we can use our tranpose() operators
real*8 :: subcube(nx,ny,nz)
real*8 :: subcube_s(nx,ny,nz)
real*8 :: subcube_t(g_nz2,nslabx,ny_2dz)
real*8 :: subcube_st(g_nz2,nslabx,ny_2dz)

!local
real*8 :: rhat(3),rvec(3),rperp1(3),rperp2(3),delu(3),dir_shift(3)
real*8 :: eta,lambda,r_lambda,ke_diss
real*8 :: dummy,xtmp,ntot
integer :: idir,idel,i2,j2,k2,i,j,k,n,m,ishift,k_g,j_g
integer :: n1,n1d,n2,n2d,n3,n3d,ierr,p

if (firstcall) then
   firstcall=.false.
   if (ncpu_x*ncpu_y>1) then
      call abort("isoave: requires x-y hyperslab parallel decomposition")
   endif   
   call init
endif



call zero_str

subcube=0

ke_diss=0
ke=0
ntot=0
do k=nz1,nz2
   if (zcord(k)>=range(3,1) .and. zcord(k)<range(3,2)) then
      do j=ny1,ny2
         if (ycord(j)>=range(2,1) .and. ycord(j)<range(2,2)) then
            do i=nx1,nx2
               if (xcord(i)>=range(1,1) .and. xcord(i)<range(1,2)) then
                  subcube(i,j,k)=1
                  ntot=ntot+1
               endif
            enddo
         endif
      enddo
   endif
enddo


do n=1,ndim
   do m=1,ndim
      call der(Q(1,1,1,n),Qs(1,1,1,1),dummy,Qs(1,1,1,2),1,m)
      do k=nz1,nz2
         do j=ny1,ny2
            do i=nx1,nx2
               
               if (subcube(i,j,k)/=0) then
                  if (m==1) ke = ke + .5*Q(i,j,k,n)**2
                  ke_diss=ke_diss + Qs(i,j,k,1)*Qs(i,j,k,1)
               endif
               
            enddo
         enddo
      enddo
   enddo
enddo


#ifdef USE_MPI
   xtmp=ke_diss
   call MPI_allreduce(xtmp,ke_diss,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   xtmp=ke
   call MPI_allreduce(xtmp,ke,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   xtmp=ntot
   call MPI_allreduce(xtmp,ntot,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif

epsilon=mu*ke_diss/ntot
ke=ke/ntot



eta = (mu**3 / epsilon)**.25
lambda=sqrt(10*ke*mu/epsilon)       ! single direction lambda
R_lambda = lambda*sqrt(2*ke/3)/mu 

if (my_pe==io_pe) then
   print *,'subcube range:'
   print *,'x: ',range(1,:)
   print *,'y: ',range(2,:)
   print *,'z: ',range(3,:)
   print *,'ntot:    ',ntot
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
call transpose_to_z(subcube,subcube_t,n1,n1d,n2,n2d,n3,n3d)
   



do idir=1,ndir

   if (my_pe==io_pe) then
      write(*,'(a,i3,a,i3,a,3i3,a)') 'direction: ',idir,'/',ndir,'  (',dir(:,idir),')'
   endif

      rhat = dir(:,idir)
      rhat=rhat/sqrt(rhat(1)**2+rhat(2)**2+rhat(3)**2)
      call compute_perp(rhat,rperp1,rperp2)


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
         call comp_str_subcube_xy(Q,idir,rhat,rperp1,rperp2,dir_shift,subcube)
      else if (dir_shift(2)==0) then
         ! no need to shift, y-direction alread 0
         call comp_str_subcube_xz(Qt,idir,rhat,rperp1,rperp2,dir_shift,subcube_t)
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
            subcube_s(i,j,k)=subcube(i,j2,k)
         enddo
         enddo
         enddo
         do n=1,3
            call transpose_to_z(Qs(1,1,1,n),Qst(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
         enddo
         call transpose_to_z(subcube_s,subcube_st,n1,n1d,n2,n2d,n3,n3d)
         call comp_str_subcube_xz(Qst,idir,rhat,rperp1,rperp2,dir_shift,subcube_st)
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
            subcube_st(k,i,j)=subcube_t(k2,i,j)
         enddo
         enddo
         enddo
         do n=1,3
            call transpose_from_z(Qst(1,1,1,n),Qs(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
         enddo
         call transpose_from_z(subcube_st,subcube_s,n1,n1d,n2,n2d,n3,n3d)
         call comp_str_subcube_xy(Qs,idir,rhat,rperp1,rperp2,dir_shift,subcube_s)
      else
         call abort("parallel computation of direction not supported")
      endif
   enddo


call normalize_and_reduce(ntot)

end subroutine





!
! same as below, but serial version
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
real*8 :: eta,lambda,r_lambda,ke_diss
real*8 :: dummy,ntot
integer :: idir,idel,i2,j2,k2,i,j,k,n,m,p

if (firstcall) then
   if (ncpu_x*ncpu_y*ncpu_z>1) then
      call abort("isoave1: can only run with 1 MPI cpu")
   endif   
   firstcall=.false.
   call init
endif

call zero_str

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
ntot=real(lz2-lz1+1)*(ly2-ly1+1)*(lx2-lx1+1)
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

         call accumulate_str(idir,idel, &
              Q(i,j,k,1),Q(i,j,k,2),Q(i,j,k,3),&
              Q(i2,j2,k2,1),Q(i2,j2,k2,2),Q(i2,j2,k2,3), &
              rhat,rperp1,rperp2)


         
      enddo
      enddo
      enddo
      endif
enddo
enddo
!$omp end parallel do

call normalize_and_reduce(ntot)


end subroutine




subroutine comp_str_xy(Q,stype,idir,rhat,rperp1,rperp2,dir_base)
use params
implicit none
!input
integer :: stype
real*8 :: Q(nx,ny,nz,*)       
real*8 :: rhat(3),rperp1(3),rperp2(3),dir_base(3)

!local
real*8 :: rvec(3),delu(3)
integer :: idir,idel,i2,j2,k2,i,j,k,n,m,p


  do idel=1,ndelta

      rvec = dir_base*delta_val(idel)
      ! dont bother computing deltas above 50% domain size
      if ( (rvec(1)**2+rvec(2)**2+rvec(3)**2) < g_nmin**2/4) then
      
      if (rvec(1)<0) rvec(1)=rvec(1)+nslabx
      if (rvec(2)<0) rvec(2)=rvec(2)+nslaby
      rvec(3)=0  ! data required to be in xy slab
      
      
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

         if (stype==2) then
            call accumulate_cj_str(idir,idel, &
              Q(i,j,k,1),Q(i,j,k,2),&
              Q(i2,j2,k,1),Q(i2,j2,k,2))
         else if (stype==1) then
            call accumulate_scalar_str(idir,idel, &
              Q(i,j,k,1),Q(i2,j2,k,1))
         else
            call accumulate_str(idir,idel, &
              Q(i,j,k,1),Q(i,j,k,2),Q(i,j,k,3),&
              Q(i2,j2,k,1),Q(i2,j2,k,2),Q(i2,j2,k,3), &
              rhat,rperp1,rperp2)
         endif

         
      enddo
      enddo
      enddo
      endif
   enddo

end subroutine









subroutine comp_str_xz(Q,stype,idir,rhat,rperp1,rperp2,dir_base)
use params
implicit none
!input
integer :: stype
real*8 :: Q(g_nz2,nslabx,ny_2dz,*)  
real*8 :: rhat(3),rperp1(3),rperp2(3),dir_base(3)

!local
real*8 :: rvec(3),delu(3)
integer :: idir,idel,i2,j2,k2,i,j,k,n,m,p


  do idel=1,ndelta

      rvec = dir_base*delta_val(idel)
      ! dont bother computing deltas above 50% domain size
      if ( (rvec(1)**2+rvec(2)**2+rvec(3)**2) < g_nmin**2/4) then
      
      if (rvec(1)<0) rvec(1)=rvec(1)+nslabx
      rvec(2)=0  ! data required to be in xz slab
      if (rvec(3)<0) rvec(3)=rvec(3)+g_nz
      
      do j=1,ny_2dz
      do i=1,nslabx
      do k=1,g_nz


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

         if (stype==2) then
            call accumulate_cj_str(idir,idel, &
              Q(k,i,j,1),Q(k,i,j,2),&
              Q(k2,i2,j,1),Q(k2,i2,j,2))
         else if (stype==2) then
            call accumulate_scalar_str(idir,idel, &
              Q(k,i,j,1),Q(k2,i2,j,1))
         else
            call accumulate_str(idir,idel, &
              Q(k,i,j,1),Q(k,i,j,2),Q(k,i,j,3),&
              Q(k2,i2,j,1),Q(k2,i2,j,2),Q(k2,i2,j,3), &
              rhat,rperp1,rperp2)
         endif

      enddo
      enddo
      enddo
      endif
   enddo

end subroutine





subroutine comp_str_subcube_xy(Q,idir,rhat,rperp1,rperp2,dir_base,subcube)
use params
implicit none
!input
real*8 :: Q(nx,ny,nz,ndim)       
real*8 :: subcube(nx,ny,nz)
real*8 :: rhat(3),rperp1(3),rperp2(3),dir_base(3)

!local
real*8 :: rvec(3),delu(3)
integer :: idir,idel,i2,j2,k2,i,j,k,n,m,p


  do idel=1,ndelta

      rvec = dir_base*delta_val(idel)
      ! dont bother computing deltas above 50% domain size
      if ( (rvec(1)**2+rvec(2)**2+rvec(3)**2) < g_nmin**2/4) then
      
      if (rvec(1)<0) rvec(1)=rvec(1)+nslabx
      if (rvec(2)<0) rvec(2)=rvec(2)+nslaby
      rvec(3)=0  ! data required to be in xy slab
      
      
      do k=nz1,nz2
      do j=ny1,ny2
      do i=nx1,nx2

         if (subcube(i,j,k)/=0) then  ! using an exit here exits *all* loops!

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

         call accumulate_str(idir,idel, &
              Q(i,j,k,1),Q(i,j,k,2),Q(i,j,k,3),&
              Q(i2,j2,k,1),Q(i2,j2,k,2),Q(i2,j2,k,3), &
              rhat,rperp1,rperp2)
      endif

      enddo
      enddo
      enddo
      endif
   enddo

end subroutine









subroutine comp_str_subcube_xz(Q,idir,rhat,rperp1,rperp2,dir_base,subcube)
use params
implicit none
!input
real*8 :: Q(g_nz2,nslabx,ny_2dz,ndim)  
real*8 :: subcube(g_nz2,nslabx,ny_2dz)
real*8 :: rhat(3),rperp1(3),rperp2(3),dir_base(3)


!local
real*8 :: rvec(3),delu(3)
integer :: idir,idel,i2,j2,k2,i,j,k,n,m,p

  do idel=1,ndelta

      rvec = dir_base*delta_val(idel)
      ! dont bother computing deltas above 50% domain size
      if ( (rvec(1)**2+rvec(2)**2+rvec(3)**2) < g_nmin**2/4) then
      
      if (rvec(1)<0) rvec(1)=rvec(1)+nslabx
      rvec(2)=0  ! data required to be in xz slab
      if (rvec(3)<0) rvec(3)=rvec(3)+g_nz


      
      do j=1,ny_2dz
      do i=1,nslabx
      do k=1,g_nz

         
         if (subcube(k,i,j)/=0) then   ! using an exit here exits *all* loops!

         
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
            

            call accumulate_str(idir,idel, &
                 Q(k,i,j,1),Q(k,i,j,2),Q(k,i,j,3),&
                 Q(k2,i2,j,1),Q(k2,i2,j,2),Q(k2,i2,j,3), &
                 rhat,rperp1,rperp2)

         endif
         
      enddo
      enddo
      enddo
      endif
   enddo

end subroutine



subroutine accumulate_str(idir,idel,u1,u2,u3,ur1,ur2,ur3,rhat,rperp1,rperp2)
!
!  (u1,u2,u3)      velocity at the point x
!  (ur1,ur2,ur3)   velocity at the point x+r
!
!
real*8 :: u1,u2,u3,ur1,ur2,ur3   
real*8 :: rhat(3),rperp1(3),rperp2(3)

real*8 :: delu1,delu2,delu3
real*8 :: u_l,u_t1,u_t2,ux_t1,ux_t2,ur_t1,ur_t2
real*8 :: u_t1_sq,u_t2_sq,u_l_sq
integer :: p,idel,idir

delu1=ur1-u1
delu2=ur2-u2
delu3=ur3-u3

u_l  = delu1*rhat(1)+delu2*rhat(2)+delu3*rhat(3)
u_t1 = delu1*rperp1(1)+delu2*rperp1(2)+delu3*rperp1(3)
u_t2 = delu1*rperp2(1)+delu2*rperp2(2)+delu3*rperp2(3)

u_t1_sq = u_t1*u_t1
u_t2_sq = u_t2*u_t2
u_l_sq  = u_l*u_l

! HELICAL structure function
ux_t1 = u1*rperp1(1)+u2*rperp1(2)+u3*rperp1(3)
ux_t2 = u1*rperp2(1)+u2*rperp2(2)+u3*rperp2(3)
ur_t1 = ur1*rperp1(1)+ur2*rperp1(2)+ur3*rperp1(3)
ur_t2 = ur1*rperp2(1)+ur2*rperp2(2)+ur3*rperp2(3)
H_ltt(idel,idir)=H_ltt(idel,idir) - u_l*(ux_t1*ur_t2-ux_t2*ur_t1)
H_tt(idel,idir)=H_tt(idel,idir) - (ux_t1*ur_t2-ux_t2*ur_t1)


Dl(idel,idir,2)  =  Dl(idel,idir,2) + u_l_sq
Dt(idel,idir,1,2)=Dt(idel,idir,1,2) + u_t1_sq
Dt(idel,idir,2,2)=Dt(idel,idir,2,2) + u_t2_sq

Dl(idel,idir,3)  =  Dl(idel,idir,3) + u_l*u_l_sq
Dt(idel,idir,1,3)=Dt(idel,idir,1,3) + u_t1*u_t1_sq
Dt(idel,idir,2,3)=Dt(idel,idir,2,3) + u_t2*u_t2_sq

Dl(idel,idir,4)  =  Dl(idel,idir,4) + u_l_sq*u_l_sq
Dt(idel,idir,1,4)=Dt(idel,idir,1,4) + u_t1_sq*u_t1_sq
Dt(idel,idir,2,4)=Dt(idel,idir,2,4) + u_t2_sq*u_t2_sq

Dl(idel,idir,5)  =  Dl(idel,idir,5) + u_l*u_l_sq*u_l_sq
Dt(idel,idir,1,5)=Dt(idel,idir,1,5) + u_t1*u_t1_sq*u_t1_sq
Dt(idel,idir,2,5)=Dt(idel,idir,2,5) + u_t2*u_t2_sq*u_t2_sq

Dl(idel,idir,6)  =  Dl(idel,idir,6) + u_l_sq*u_l_sq*u_l_sq
Dt(idel,idir,1,6)=Dt(idel,idir,1,6) + u_t1_sq*u_t1_sq*u_t1_sq
Dt(idel,idir,2,6)=Dt(idel,idir,2,6) + u_t2_sq*u_t2_sq*u_t2_sq

! this loop is very expensive.  
do p=7,pmax
   Dl(idel,idir,p)=Dl(idel,idir,p) + u_l**p
   Dt(idel,idir,1,p)=Dt(idel,idir,1,p) + u_t1**p
   Dt(idel,idir,2,p)=Dt(idel,idir,2,p) + u_t2**p
enddo

D_ltt(idel,idir,1)=D_ltt(idel,idir,1) + u_l*u_t1_sq
D_ltt(idel,idir,2)=D_ltt(idel,idir,2) + u_l*u_t2_sq
D_lltt(idel,idir)=D_lltt(idel,idir) + u_l_sq * .5*(u_t1_sq + u_t2_sq)


if (u_l>=0) then
   SP_lll(idel,idir)  =SP_lll(idel,idir) + u_l*u_l_sq
   SP_ltt(idel,idir,1)=SP_ltt(idel,idir,1) + u_l*u_t1_sq
   SP_ltt(idel,idir,2)=SP_ltt(idel,idir,2) + u_l*u_t2_sq
else
   SN_lll(idel,idir)  =SN_lll(idel,idir) - u_l*u_l_sq
   SN_ltt(idel,idir,1)=SN_ltt(idel,idir,1) - u_l*u_t1_sq
   SN_ltt(idel,idir,2)=SN_ltt(idel,idir,2) - u_l*u_t2_sq
endif
end subroutine 





subroutine accumulate_cj_str(idir,idel,u1,u2,ur1,ur2)
!
!  (u1,u2)         w(x), s(x)
!  (ur1,ur2)       w(x+r), s(x+r)
!
!
real*8 :: u1,ur1,u2,ur2
real*8 :: rhat(3),rperp1(3),rperp2(3)
integer :: idir,idel

! local
integer :: p


! compute:  <w^2w^2(r)>,  <s^2s^2(r)>, <w^2s^2(r)>  
w2s2(idel,idir,1)=w2s2(idel,idir,1)+u1*ur1
w2s2(idel,idir,2)=w2s2(idel,idir,2)+u2*ur2
w2s2(idel,idir,3)=w2s2(idel,idir,3)+u1*ur2


end subroutine 






subroutine accumulate_scalar_str(idir,idel,u1,ur1)
!
!  (u1,u2)         w(x), s(x)
!  (ur1,ur2)       w(x+r), s(x+r)
!
!
real*8 :: u1,ur1,u2,ur2
real*8 :: rhat(3),rperp1(3),rperp2(3)
integer :: idir,idel

! local
integer :: p

do p=2,pmax
   Dl(idel,idir,p)=(ur1-u1)**p
enddo







end subroutine 







!
! For each direction, find its intersection with the unit sphere
! (2 points) and output in lat-lon coordinates
!
!
subroutine writepoints()
use params


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
integer :: i,j,idel,idir,max_delta
real*8 :: rvec(3)



delta_val=100000

! compute all deltas up to ndelta_max
! (but we only use deltas up to ndelta)
! 1..16
do i=1,16
   delta_val(i)=i
enddo
j=16

! 18..32  (start with 2*9)
do i=9,16
   j=j+1
   delta_val(j)=2*i
enddo

! 36..64  (start with 4*9)
do i=9,16
   j=j+1
   delta_val(j)=4*i
enddo

! 72..128 (start with 8*9)
do i=9,16
   j=j+1
   delta_val(j)=8*i
enddo

! 144..256  (start with 16*9)
do i=9,16
   j=j+1
   delta_val(j)=16*i
enddo

! 288..512  (start with 32*9)
do i=9,16
   j=j+1
   delta_val(j)=32*i
enddo

! 576..1024  (start with 64*9)
do i=9,16
   j=j+1
   delta_val(j)=64*i
enddo

! 1152..2048
do i=9,16
   j=j+1
   delta_val(j)=128*i
enddo
! determine maximum value of delta to use for this grid
if (j > ndelta_max) then
   stop "isoave init: j > ndelta_max"
endif

max_delta = g_nmin/2
ndir=ndir_max

if (user_specified_isodel>0) then
   max_delta=min(max_delta,user_specified_isodel)
endif
if (user_specified_isodir>0) then
   ndir=min(user_specified_isodir,ndir_max)
endif



do idel=1,ndelta_max
   if (delta_val(idel) >= max_delta) exit
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



if (ndir_max<73) then
   call abort("isoave: ndir_max <  73")
endif



allocate(dwork2(ndelta,ndir))
allocate(dwork3(ndelta,ndir,2))

allocate(Dl(ndelta,ndir,2:pmax))
allocate(Dt(ndelta,ndir,2,2:pmax))
allocate(D_ltt(ndelta,ndir,2))
allocate(D_lltt(ndelta,ndir))

allocate(H_ltt(ndelta,ndir))
allocate(H_tt(ndelta,ndir))
allocate(SP_lll(ndelta,ndir))
allocate(SP_ltt(ndelta,ndir,2))

allocate(SN_lll(ndelta,ndir))
allocate(SN_ltt(ndelta,ndir,2))

allocate(w2s2(ndelta,ndir,3))


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






subroutine normalize_and_reduce(ntot)
!
! noramlize all structure functions by ntot, then
! to a global sum over all cpus.  
!
use params
use mpi
implicit none
real*8 :: ntot

integer :: ierr,p
Dl=Dl/ntot
Dt=Dt/ntot
D_lltt=D_lltt/ntot
D_ltt=D_ltt/ntot
H_ltt=H_ltt/ntot
H_tt=H_tt/ntot

SP_ltt=SP_ltt/ntot
SP_lll=SP_lll/ntot
SN_ltt=SN_ltt/ntot
SN_lll=SN_lll/ntot

w2s2=w2s2/ntot

#ifdef USE_MPI
   do p=2,pmax
   dwork2=Dl(:,:,p)
   call MPI_reduce(dwork2,Dl(1,1,p),ndelta*ndir,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   dwork3=Dt(:,:,:,p)
   call MPI_reduce(dwork3,Dt(1,1,1,p),ndelta*ndir*2,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   enddo

   dwork2=H_ltt
   call MPI_reduce(dwork2,H_ltt,ndelta*ndir,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   dwork2=H_tt
   call MPI_reduce(dwork2,H_tt,ndelta*ndir,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)


   dwork3=D_ltt
   call MPI_reduce(dwork3,D_ltt,ndelta*ndir*2,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   dwork2=D_lltt
   call MPI_reduce(dwork3,D_lltt,ndelta*ndir,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)


   dwork3=SP_ltt
   call MPI_reduce(dwork3,SP_ltt,ndelta*ndir*2,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   dwork2=SP_lll
   call MPI_reduce(dwork2,SP_lll,ndelta*ndir,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   dwork3=SN_ltt
   call MPI_reduce(dwork3,SN_ltt,ndelta*ndir*2,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   dwork2=SN_lll
   call MPI_reduce(dwork2,SN_lll,ndelta*ndir,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)

   do p=1,3
      dwork2=w2s2(:,:,p)
      call MPI_reduce(dwork2,w2s2(1,1,p),ndelta*ndir,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   enddo


#endif
end subroutine


subroutine zero_str
Dl=0
Dt=0
D_ltt=0
D_lltt=0
H_ltt=0
H_tt=0
SP_ltt=0
SP_lll=0
SN_ltt=0
SN_lll=0
w2s2=0
end subroutine

end module 
