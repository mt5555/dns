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


integer,parameter :: ndir_max=73
integer,parameter :: ndelta_max=100
integer,parameter :: pmax=10       ! has to be 6 or greater

!
!  str_type=0    standard structure functions (including helical)
!  str_type=1    fractional structure functions power < 1
!  str_type=2    fractional structure functions power > 1
!  str_type=3    fractional structure functions power (see subroutine)
!  str_type=4    anisotropic structure functions
!
integer :: str_type = 4

real*8 :: fractional_power(2:pmax) 
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

real*8 :: u_shear(3,3)    ! stress tensor
real*8 :: u_shear2(3,3)    ! copy of stress tensor

real*8,allocatable  :: dwork2(:,:)
real*8,allocatable  :: dwork3(:,:,:)


logical  :: use_max_shear_direction=.true.   ! requires str_type==4
integer :: idir_max
real*8  :: t1(3),t2(3)

! also added to the file for completeness:
real*8,private :: epsilon,mu,ke,h_epsilon=0

private init

contains




subroutine write_ux(fid,time)
use params
implicit none

CPOINTER fid
real*8 :: time
!local
integer :: i,idir,p
real*8 :: x
                                                               
if (my_pe/=io_pe) return

call cwrite8(fid,time,1)
call cwrite8(fid,u_shear,9)
x=idir_max
call cwrite8(fid,x,1)

end subroutine



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
x=8;      call cwrite8(fid,x,1)   ! number of scalars
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

! NOTE: if you add scalars here, be sure to update
! the number of scalars field in the header above.
call cwrite8(fid,time,1)
x=g_nx; call cwrite8(fid,x,1)
x=g_ny; call cwrite8(fid,x,1)
x=g_nz; call cwrite8(fid,x,1)
call cwrite8(fid,mu,1)
call cwrite8(fid,ke,1)
call cwrite8(fid,epsilon,1)
call cwrite8(fid,h_epsilon,1)


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
!  stype=4   assume Q is (u,v,w) and 1 scalar q.  compute standard velocity
!                structure functions and mixed velocity q structure functions
!
!  stype=3   assume Q is (u,v,w) data.  compute standard velocity
!                                    structure functions
!
!  stype=2   assume Q = (w,s), the vorticity and strain.  compute CJ's 
!            correlations, of the form <w(x)s(x+r)>
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
real*8 :: Qst(g_nz2,nslabx,ny_2dz,n_var)  ! transpose

! these two arrays can be overlapped in memory:
real*8 :: Qt(g_nz2,nslabx,ny_2dz,n_var)   ! transpose
real*8 :: Qs(nx,ny,nz,n_var)              ! shifted original data





!local
real*8 :: rhat(3),rvec(3),rperp1(3),rperp2(3),delu(3),dir_shift(3)
real*8 :: u_l,u_t1,u_t2,rnorm
real*8 :: eta,lambda,r_lambda,ke_diss,h_diss,xvec(3),xmin
real*8 :: dummy(pmax),xtmp,ntot,xfac,ux,uy,uz,vx,vy,vz,wx,wy,wz
character(len=80) :: message
integer :: idir,idel,i2,j2,k2,i,j,k,n,m,ishift,k_g,j_g,nd
integer :: n1,n1d,n2,n2d,n3,n3d,ierr,p,csig
integer :: im,jm,km
logical :: qt_uptodate

if (firstcall) then
   firstcall=.false.
   if (ncpu_x*ncpu_y>1) then
      call abort("isoave: requires x-y hyperslab parallel decomposition")
   endif   
   call init
endif

nd=3  ! number of variables stored in Q.
if (stype==4) nd=4
if (stype==2) nd=2
if (stype==1) nd=1

if (nd > n_var) then
   call abort("n_var insufficint for chosen stype/nd structure functions")
endif

ntranspose=0
ntot=real(g_nx)*g_ny*g_nz


call zero_str

if (stype==3) then
ke_diss=0
ke=0
u_shear=0
do n=1,ndim
   do m=1,ndim
      call der(Q(1,1,1,n),Qs(1,1,1,1),dummy,Qs(1,1,1,2),1,m)
      do k=nz1,nz2
         do j=ny1,ny2
            do i=nx1,nx2
               
               if (m==1) ke = ke + .5*Q(i,j,k,n)**2
               ke_diss=ke_diss + Qs(i,j,k,1)*Qs(i,j,k,1)
               u_shear(n,m) = u_shear(n,m)+ Qs(i,j,k,1)
            enddo
         enddo
      enddo
   enddo
enddo


Qs=Q
do n=1,3
   call fft3d(Qs(1,1,1,n),Qst)  ! use Qst as a work array
enddo

h_diss = 0
do k=nz1,nz2
   km=kmcord(k)
   do j=ny1,ny2
      jm=jmcord(j)
      do i=nx1,nx2
         im=imcord(i)
         
         xfac=-(im*im + jm*jm + km*km)*pi2_squared
         xfac = 2*2*2*xfac
         if (kmcord(k)==0) xfac=xfac/2
         if (jmcord(j)==0) xfac=xfac/2
         if (imcord(i)==0) xfac=xfac/2

         ! u_x term
         !ux = - pi2*im*Qs(i+imsign(i),j,k,1)
         vx = - pi2*im*Qs(i+imsign(i),j,k,2)
         wx = - pi2*im*Qs(i+imsign(i),j,k,3)
         
         uy = - pi2*jm*Qs(i,j+jmsign(j),k,1)
         !vy = - pi2*jm*Qs(i,j+jmsign(j),k,2)
         wy = - pi2*jm*Qs(i,j+jmsign(j),k,3)
         
         uz =  - pi2*km*Qs(i,j,k+kmsign(k),1)
         vz =  - pi2*km*Qs(i,j,k+kmsign(k),2)
         !wz =  - pi2*km*Qs(i,j,k+kmsign(k),3)
      
         ! vorticity: ( (wy - vz), (uz - wx), (vx - uy) )

         ! compute 2*k^2 u vor:
         h_diss = h_diss + 2*xfac*(Qs(i,j,k,1)*(wy-vz) + &
                                   Qs(i,j,k,2)*(uz-wx) + &
                                   Qs(i,j,k,3)*(vx-uy)) 

         ! ke_diss = ke_diss + xfac*(Qs(i,j,k,1)**2 + &
         !                            Qs(i,j,k,2)**2 + &
         !                            Qs(i,j,k,3)**2 )
         
      enddo
   enddo
enddo





#ifdef USE_MPI
   xtmp=ke_diss
   call mpi_allreduce(xtmp,ke_diss,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   xtmp=ke
   call mpi_allreduce(xtmp,ke,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   xtmp=h_diss
   call mpi_allreduce(xtmp,h_diss,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   u_shear2=u_shear
   call mpi_allreduce(u_shear2,u_shear,9,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif


epsilon=mu*ke_diss/ntot
if (epsilon==0) epsilon=1e-20
h_epsilon=mu*h_diss
if (h_epsilon==0) h_epsilon=1e-20
ke=ke/ntot
u_shear=u_shear/ntot
if (my_pe==io_pe) then
   write(6,*)'u_shear = ', u_shear
endif

eta = (mu**3 / epsilon)**.25
lambda=sqrt(10*ke*mu/epsilon)       ! single direction lambda
R_lambda = lambda*sqrt(2*ke/3)/mu 

if (my_pe==io_pe) then
   print *,'KE:      ',ke
   print *,'epsilon:   ',epsilon
   print *,'h-epsilon: ',h_epsilon
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
   w2s2_mean=w2s2_mean/ntot
#ifdef USE_MPI
   dummy(1:3)=w2s2_mean
   call mpi_allreduce(dummy,w2s2_mean,3,MPI_REAL8,MPI_SUM,comm_3d,ierr)
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
   Dl_mean=Dl_mean/ntot
#ifdef USE_MPI
   dummy(2:pmax)=Dl_mean(2:pmax)
   call mpi_allreduce(dummy(2),Dl_mean(2),pmax-1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
endif

idir_max=-1
if (use_max_shear_direction) then 
   call max_shear_coordinate_system(u_shear,idir_max,t1,t2)
endif

   
qt_uptodate=.false.
Do idir=1,ndir


!  check for SIGURG, telling us to stop because job will soon be killed
   call caught_sig(csig); 
#ifdef USE_MPI
   i=csig;
   call mpi_allreduce(i,csig,1,MPI_INTEGER,MPI_MAX,comm_3d,ierr)
#endif
   if (csig>0) exit;

      
   rhat = dir(:,idir)
   rhat=rhat/sqrt(rhat(1)**2+rhat(2)**2+rhat(3)**2)
   if (idir==idir_max) then
      rperp1=t1; rperp2=t2
   else
      call compute_perp(rhat,rperp1,rperp2)
   endif

   if (my_pe==io_pe) then
      xvec=rperp2;
      xmin=99;
      do i=1,3
         if (xvec(i)/=0 .and. abs(xvec(i))< xmin ) then
            xmin=abs(xvec(i))
         endif
      enddo
      xvec=xvec/xmin
      write(*,'(a,i3,a,i3,a,3i3,a,a,3f6.2,a,a,f7.2)') 'direction: ',idir,'/',ndir,&
           '  (',dir(:,idir),')' ,'  t2=(',xvec,')',' <t2,t2>=',sum(xvec*xvec)
   endif

#if 1
      ! check orthoginality
   if(my_pe==io_pe) then
      print *,'norms: ',sqrt(rhat(1)**2+rhat(2)**2+rhat(3)**2), &
           sqrt(rperp1(1)**2+rperp1(2)**2+rperp1(3)**2), &
           sqrt(rperp2(1)**2+rperp2(2)**2+rperp2(3)**2)
      
      print *,'ortho: ',&
           rhat(1)*rperp1(1)+rhat(2)*rperp1(2)+rhat(3)*rperp1(3), &
           rhat(1)*rperp2(1)+rhat(2)*rperp2(2)+rhat(3)*rperp2(3), &
           rperp2(1)*rperp1(1)+rperp2(2)*rperp1(2)+rperp2(3)*rperp1(3)
   endif
#endif
      xmin=abs(rhat(1)*rperp1(1)+rhat(2)*rperp1(2)+rhat(3)*rperp1(3))+ &
           abs(rhat(1)*rperp2(1)+rhat(2)*rperp2(2)+rhat(3)*rperp2(3)) + &
           abs(rperp2(1)*rperp1(1)+rperp2(2)*rperp1(2)+rperp2(3)*rperp1(3))
      if (xmin>1e-15) call abort("isaove.F90: orthogonality error")

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
         if (.not. qt_uptodate) then
            do n=1,nd
               ntranspose=ntranspose+1
               call transpose_to_z(Q(1,1,1,n),Qt(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
            enddo   
            qt_uptodate=.true.
         endif
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

         ! Qs is same array as Qt, so we just trashed Qt:
         qt_uptodate=.false.	

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
         if (.not. qt_uptodate) then
            do n=1,nd
               ntranspose=ntranspose+1
               call transpose_to_z(Q(1,1,1,n),Qt(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
            enddo   
            qt_uptodate=.true.
         endif

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
         ! Qs is the same array as Qt, so we just trashed Qt
         qt_uptodate=.false.
         call comp_str_xy(Qs,stype,idir,rhat,rperp1,rperp2,dir_shift)
      else
         call abort("parallel computation of direction not supported")
      endif

100 continue
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

#if 0
u_shear=0
u_shear(1,1)=-1
u_shear(2,2)=2
u_shear(3,3)=-1
u_shear(1,2)=1
u_shear(1,3)=.5
call max_shear_coordinate_system(u_shear,idir_max,t1,t2)
stop
#endif


call zero_str

subcube=0

ke_diss=0
ke=0
u_shear=0

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
                  u_shear(n,m) = u_shear(n,m)+ Qs(i,j,k,1)
               endif
               
            enddo
         enddo
      enddo
   enddo
enddo


#ifdef USE_MPI
   xtmp=ke_diss
   call mpi_allreduce(xtmp,ke_diss,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   xtmp=ke
   call mpi_allreduce(xtmp,ke,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   xtmp=ntot
   call mpi_allreduce(xtmp,ntot,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   u_shear2=u_shear
   call mpi_allreduce(u_shear2,u_shear,9,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif

epsilon=mu*ke_diss/ntot
ke=ke/ntot
u_shear=u_shear/ntot

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
   

idir_max=-1
if (use_max_shear_direction) then 
   call max_shear_coordinate_system(u_shear,idir_max,t1,t2)
endif


do idir=1,ndir

   if (my_pe==io_pe) then
      write(*,'(a,i3,a,i3,a,3i3,a,a,3i3,a)') 'direction: ',idir,'/',ndir,&
           '  (',dir(:,idir),')' 
   endif


   rhat = dir(:,idir)
   rhat=rhat/sqrt(rhat(1)**2+rhat(2)**2+rhat(3)**2)
   if (idir==idir_max) then
      rperp1=t1; rperp2=t2
   else
      call compute_perp(rhat,rperp1,rperp2)
   endif



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
u_shear=0
do n=1,ndim
   do m=1,ndim
      call der(Q(1,1,1,n),d1,dummy,work,1,m)
      do k=lz1,lz2
      do j=ly1,ly2
      do i=lx1,lx2
         if (m==1) ke = ke + .5*Q(i,j,k,n)**2
         ke_diss=ke_diss + d1(i,j,k)*d1(i,j,k)
         u_shear(n,m) = u_shear(n,m)+ d1(i,j,k)
      enddo
      enddo
      enddo
   enddo
enddo
ntot=real(lz2-lz1+1)*(ly2-ly1+1)*(lx2-lx1+1)
epsilon=mu*ke_diss/ntot
ke=ke/ntot
u_shear=u_shear/ntot

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



      

idir_max=-1
if (use_max_shear_direction) then 
   call max_shear_coordinate_system(u_shear,idir_max,t1,t2)
endif


do idir=1,ndir

   if (my_pe==io_pe) then
      write(*,'(a,i3,a,i3,a,3i3,a,a,3i3,a)') 'direction: ',idir,'/',ndir,&
           '  (',dir(:,idir),')' 
   endif


   rhat = dir(:,idir)
   rhat=rhat/sqrt(rhat(1)**2+rhat(2)**2+rhat(3)**2)
   if (idir==idir_max) then
      rperp1=t1; rperp2=t2
   else
      call compute_perp(rhat,rperp1,rperp2)
   endif


#if 1
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
         else if (stype==3) then
            call accumulate_str(idir,idel, &
              Q(i,j,k,1),Q(i,j,k,2),Q(i,j,k,3),&
              Q(i2,j2,k,1),Q(i2,j2,k,2),Q(i2,j2,k,3), &
              rhat,rperp1,rperp2)
         else if (stype==4) then
            call accumulate_str_uq(idir,idel, &
              Q(i,j,k,1),Q(i,j,k,2),Q(i,j,k,3),Q(i,j,k,4),&
              Q(i2,j2,k,1),Q(i2,j2,k,2),Q(i2,j2,k,3),Q(i2,j2,k,4), &
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
         else if (stype==1) then
            call accumulate_scalar_str(idir,idel, &
              Q(k,i,j,1),Q(k2,i2,j,1))
         else if (stype==3) then
            call accumulate_str(idir,idel, &
              Q(k,i,j,1),Q(k,i,j,2),Q(k,i,j,3),&
              Q(k2,i2,j,1),Q(k2,i2,j,2),Q(k2,i2,j,3), &
              rhat,rperp1,rperp2)
         else if (stype==4) then
            call accumulate_str_uq(idir,idel, &
              Q(k,i,j,1),Q(k,i,j,2),Q(k,i,j,3),Q(k,i,j,4),&
              Q(k2,i2,j,1),Q(k2,i2,j,2),Q(k2,i2,j,3),Q(k2,i2,j,4), &
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
real*8 :: u_t1_sq,u_t2_sq,u_l_sq,xp
real*8 :: u_l_3,u_t1_3,u_t2_3
integer :: p,idel,idir
real*8 :: y2_2, y22, y2_1, y21, y20, slt1, slt2, stt


delu1=ur1-u1
delu2=ur2-u2
delu3=ur3-u3

u_l  = delu1*rhat(1)+delu2*rhat(2)+delu3*rhat(3)
u_t1 = delu1*rperp1(1)+delu2*rperp1(2)+delu3*rperp1(3)
u_t2 = delu1*rperp2(1)+delu2*rperp2(2)+delu3*rperp2(3)


if (str_type==0 ) then
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
   do p=7,min(10,pmax)
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
   
endif

if (str_type==1) then
   xp=fractional_power(2)
   u_l=abs(u_l)**xp            ! ** .1
   u_t1=abs(u_t1)**xp
   u_t2=abs(u_t2)**xp
   
   u_l_sq=u_l*u_l              ! **.2
   u_t1_sq=u_t1*u_t1
   u_t2_sq=u_t2*u_t2
   
   Dl(idel,idir,2)  =  Dl(idel,idir,2) + u_l       ! ** .1
   Dt(idel,idir,1,2)=Dt(idel,idir,1,2) + u_t1
   Dt(idel,idir,2,2)=Dt(idel,idir,2,2) + u_t2
   
   Dl(idel,idir,3)  =  Dl(idel,idir,3) + u_l_sq    ! ** .2
   Dt(idel,idir,1,3)=Dt(idel,idir,1,3) + u_t1_sq
   Dt(idel,idir,2,3)=Dt(idel,idir,2,3) + u_t2_sq
   
   Dl(idel,idir,4)  =  Dl(idel,idir,4) + u_l*u_l_sq   ! ** .3
   Dt(idel,idir,1,4)=Dt(idel,idir,1,4) + u_t1*u_t1_sq
   Dt(idel,idir,2,4)=Dt(idel,idir,2,4) + u_t2*u_t2_sq
   
   Dl(idel,idir,5)  =  Dl(idel,idir,5) + u_l_sq*u_l_sq   ! ** .4
   Dt(idel,idir,1,5)=Dt(idel,idir,1,5) + u_t1_sq*u_t1_sq
   Dt(idel,idir,2,5)=Dt(idel,idir,2,5) + u_t2_sq*u_t2_sq
   
   Dl(idel,idir,6)  =  Dl(idel,idir,6) + u_l*u_l_sq*u_l_sq    ! ** .5
   Dt(idel,idir,1,6)=Dt(idel,idir,1,6) + u_t1*u_t1_sq*u_t1_sq
   Dt(idel,idir,2,6)=Dt(idel,idir,2,6) + u_t2*u_t2_sq*u_t2_sq
   
   Dl(idel,idir,7)  =  Dl(idel,idir,7) + u_l_sq*u_l_sq*u_l_sq   ! ** .6
   Dt(idel,idir,1,7)=Dt(idel,idir,1,7) + u_t1_sq*u_t1_sq*u_t1_sq
   Dt(idel,idir,2,7)=Dt(idel,idir,2,7) + u_t2_sq*u_t2_sq*u_t2_sq
   
   ! .1**7, .1**8, .1**9
   do p=8,min(10,pmax)   
      Dl(idel,idir,p)=Dl(idel,idir,p) + u_l**(p-1)             
      Dt(idel,idir,1,p)=Dt(idel,idir,1,p) + u_t1**(p-1)
      Dt(idel,idir,2,p)=Dt(idel,idir,2,p) + u_t2**(p-1)
   enddo
endif
if (str_type==2) then
   u_l=abs(u_l)
   u_t1=abs(u_t1)
   u_t2=abs(u_t2)
   
   u_l_sq=u_l**.25
   u_t1_sq=u_t1**.25
   u_t2_sq=u_t2**.25
   
   u_l_3=u_l_sq*u_l_sq                       ! .5
   u_t1_3=u_t1_sq*u_t1_sq
   u_t2_3=u_t2_sq*u_t2_sq
   
   Dl(idel,idir,6)  =  Dl(idel,idir,6) + u_l       ! ** 1
   Dt(idel,idir,1,6)=Dt(idel,idir,1,6) + u_t1
   Dt(idel,idir,2,6)=Dt(idel,idir,2,6) + u_t2
   
   Dl(idel,idir,7)  =  Dl(idel,idir,7) + u_l*u_l_sq       ! ** 1.25
   Dt(idel,idir,1,7)=Dt(idel,idir,1,7) + u_t1*u_t1_sq
   Dt(idel,idir,2,7)=Dt(idel,idir,2,7) + u_t2*u_t2_sq
   
   Dl(idel,idir,8)  =  Dl(idel,idir,8) + u_l*u_l_3       ! ** 1.50
   Dt(idel,idir,1,8)=Dt(idel,idir,1,8) + u_t1*u_t1_3
   Dt(idel,idir,2,8)=Dt(idel,idir,2,8) + u_t2*u_t2_3
   
   Dl(idel,idir,9)  =  Dl(idel,idir,9) + u_l*u_l_sq*u_l_3       ! ** 1.75
   Dt(idel,idir,1,9)=Dt(idel,idir,1,9) + u_t1*u_t1_sq*u_t1_3
   Dt(idel,idir,2,9)=Dt(idel,idir,2,9) + u_t2*u_t2_sq*u_t2_3
   
   
   Dl(idel,idir,10)  =  Dl(idel,idir,10) + u_l*u_l       ! ** 2
   Dt(idel,idir,1,10)=Dt(idel,idir,1,10) + u_t1*u_t1
   Dt(idel,idir,2,10)=Dt(idel,idir,2,10) + u_t2*u_t2
   
   
   u_l=(u_l)**(-.2)
   u_t1=(u_t1)**(-.2)
   u_t2=(u_t2)**(-.2)
   
   u_l_sq=u_l*u_l              ! **.4
   u_t1_sq=u_t1*u_t1
   u_t2_sq=u_t2*u_t2
   
   Dl(idel,idir,2)  =  Dl(idel,idir,2) + u_l_sq*u_l_sq   ! ** -.8
   Dt(idel,idir,1,2)=Dt(idel,idir,1,2) + u_t1_sq*u_t1_sq
   Dt(idel,idir,2,2)=Dt(idel,idir,2,2) + u_t2_sq*u_t2_sq
   
   Dl(idel,idir,3)  =  Dl(idel,idir,3) + u_l*u_l_sq   ! ** -.6
   Dt(idel,idir,1,3)=Dt(idel,idir,1,3) + u_t1*u_t1_sq
   Dt(idel,idir,2,3)=Dt(idel,idir,2,3) + u_t2*u_t2_sq
   
   Dl(idel,idir,4)  =  Dl(idel,idir,4) + u_l_sq    ! ** -.4
   Dt(idel,idir,1,4)=Dt(idel,idir,1,4) + u_t1_sq
   Dt(idel,idir,2,4)=Dt(idel,idir,2,4) + u_t2_sq
   
   Dl(idel,idir,5)  =  Dl(idel,idir,5) + u_l       ! ** -.2
   Dt(idel,idir,1,5)=Dt(idel,idir,1,5) + u_t1
   Dt(idel,idir,2,5)=Dt(idel,idir,2,5) + u_t2
endif
if (str_type==3) then
   ! 2.45 2.55 2.65 2.75 2.85 2.95 3.05 3.15 3.25 
   u_l=abs(u_l)
   
   u_l_sq=u_l**.1
   u_l_3=u_l**.2
   
   u_l=u_l**2.85
   
   Dl(idel,idir,6)  =  Dl(idel,idir,6)  + u_l                       ! 2.85
   Dl(idel,idir,7)  =  Dl(idel,idir,7)  + u_l*u_l_sq                ! 2.95
   Dl(idel,idir, 8) =  Dl(idel,idir, 8) + u_l*u_l_3                ! 3.05
   Dl(idel,idir, 9) =  Dl(idel,idir, 9) + u_l*u_l_sq*u_l_3         ! 3.15
   Dl(idel,idir,10) =  Dl(idel,idir,10) + u_l*u_l_3*u_l_3          ! 3.25
   
   
   u_l_sq=1/u_l_sq   ! -.1
   u_l_3=1/u_l_3     ! -.2
   
   Dl(idel,idir,5)  =  Dl(idel,idir,5)  + u_l*u_l_sq                    ! 2.75
   Dl(idel,idir,4)  =  Dl(idel,idir,4)  + u_l*u_l_3                     ! 2.65
   Dl(idel,idir,3)  =  Dl(idel,idir,3)  + u_l*u_l_sq*u_l_3              ! 2.55
   Dl(idel,idir,2)  =  Dl(idel,idir,2)  + u_l*u_l_3*u_l_3               ! 2.45
   
endif

! the mixed structure function $S_{ij}$, $i\neq j$

if (str_type==4) then
   slt1 = u_l*u_t1
   slt2 = u_l*u_t2	 
   stt = u_t1*u_t2
   
   !the polar angle is t and the azimuthal angle  is p
   !   y2_2 = 2*rhat(1)*rhat(2)              !sint*sint*sin2p
   !   y22 = rhat(1)**2 - rhat(2)**2         !sint*sint*cos2p
   !   y2_1 = rhat(1)*rhat(3)                !sint*cost*cosp
   !   y21 = rhat(2)*rhat(3)                 !sint*cost*sinp
   !   y20 = (3*rhat(3)**2 - 1)              ! 3*cost*cost - 1

   ! note: Dl() Dt() arrays are indexed 2:pmax - they start at 2
!   Dl(idel,idir,2) = Dl(idel,idir,2) + y2_2*(slt) 
!   Dl(idel,idir,3) = Dl(idel,idir,3) +  y22*(slt)
!   Dl(idel,idir,4) = Dl(idel,idir,4) + y2_1*(slt)
!   Dl(idel,idir,5) = Dl(idel,idir,5) + y21*(slt)
!   Dl(idel,idir,6) = Dl(idel,idir,6) + y20*(slt)
   Dl(idel,idir,2) = Dl(idel,idir,2) + slt1	!for matlab
   Dl(idel,idir,3) = Dl(idel,idir,3) + slt2	

!   Dt(idel,idir,1,2) = Dt(idel,idir,1,2) + y2_2*(stt)
!   Dt(idel,idir,1,3) = Dt(idel,idir,1,3) + y22*(stt)
!   Dt(idel,idir,1,4) = Dt(idel,idir,1,4) + y2_1*(stt)
!   Dt(idel,idir,1,5) = Dt(idel,idir,1,5) + y21*(stt)
!   Dt(idel,idir,1,6) = Dt(idel,idir,1,6) + y20*(stt)
   Dt(idel,idir,1,2) = Dt(idel,idir,1,2) + stt	!for matlab
	
endif

end subroutine 




subroutine accumulate_str_uq(idir,idel,u1,u2,u3,q1,ur1,ur2,ur3,qr,rhat,rperp1,rperp2)
!
!  (u1,u2,u3)      velocity at the point x
!  (ur1,ur2,ur3)   velocity at the point x+r
!      q           scalar at the point x
!      qr          scalar at the point x+r
!
!
real*8 :: u1,u2,u3,ur1,ur2,ur3,qr,q1   
real*8 :: rhat(3),rperp1(3),rperp2(3)

real*8 :: delu1,delu2,delu3
real*8 :: u_l,u_t1,u_t2,ux_t1,ux_t2,ur_t1,ur_t2
real*8 :: u_t1_sq,u_t2_sq,u_l_sq,xp
real*8 :: u_l_3,u_t1_3,u_t2_3
integer :: p,idel,idir
real*8 :: y2_2, y22, y2_1, y21, y20, slt1, slt2, stt

! compute the velocit increments:
delu1=ur1-u1
delu2=ur2-u2
delu3=ur3-u3

! compute longitudinal and transverse components
u_l  = delu1*rhat(1)+delu2*rhat(2)+delu3*rhat(3)
u_t1 = delu1*rperp1(1)+delu2*rperp1(2)+delu3*rperp1(3)
u_t2 = delu1*rperp2(1)+delu2*rperp2(2)+delu3*rperp2(3)

! note: Dl() Dt() arrays are indexed 2:pmax - they start at 2

!  Beth's third order correlation of q and velocity
!   Dl(idel,idir,4) = Dl(idel,idir,4) + q1*qr*u_l
Dl(idel,idir,2) = Dl(idel,idir,2) + q1*qr*u_l
	

end subroutine 


subroutine accumulate_cj_str(idir,idel,u1,u2,ur1,ur2)
!
!  (u1,u2)         w(x), s(x)
!  (ur1,ur2)       w(x+r), s(x+r)
!
!
implicit none
real*8 :: u1,ur1,u2,ur2
real*8 :: rhat(3),rperp1(3),rperp2(3)
integer :: idir,idel

! local
integer :: p


w2s2(idel,idir,1)=w2s2(idel,idir,1)+u1*ur1      ! w(x)*w(x+r)
w2s2(idel,idir,2)=w2s2(idel,idir,2)+u2*ur2      ! s(x)*s(x+r)
w2s2(idel,idir,3)=w2s2(idel,idir,3)+u1*ur2      ! w(x)*s(x+r)


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

if (str_type/=4 .and. use_max_shear_direction) then
   call abort("Error: isoave.F90 init():  use_max_shear_direction requires str_type==4")
endif

! fractional powers
if (pmax<10) then
   call abort("pmax set too small") 
endif

fractional_power(2)=.1
fractional_power(3)=.2
fractional_power(4)=.3
fractional_power(5)=.4
fractional_power(6)=.5
fractional_power(7)=.6
fractional_power(8)=.7
fractional_power(9)=.8
fractional_power(10)=.9


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
   call abort("isoave init: j > ndelta_max")
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

if (io_pe==my_pe) then
   print *,'isoave: ndir,ndelta: ',ndir,ndelta
endif



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

! 2,2,1 directions  norm^2 = 9
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


! face 1,3 directions  norm^2 = 10
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

! body 1,1,3 directions  norm^2 =  11
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

real*8 r(3),rp1(3),rp2(3),rmax
integer :: numzero,i,j,k,ii

#if 0
numzero=0
do i=1,3
   if (r(i)==0) numzero=numzero+1
enddo

if (numzero==2) then
   if (r(3)==0) then
       i=1  ! either 1 or 2 is non zero
       j=2
       k=3  ! this one is zero 
   else
       i=2  ! this one is zero
       j=3 
       k=1  ! this one is zero   
   endif
else if (numzero < 2 ) then
   ! two non zero entries.  set largest to k
   rmax=-1e99
   do ii=1,3
      if ( abs(r(ii))>rmax ) then
         k=ii
         i=ii+1; i=mod(i-1,3)+1
         j=ii+2; j=mod(j-1,3)+1
      endif
   enddo

else
   call abort("compute_perp(): this is not possible")
endif
rp1(i)=-r(j)
rp1(j)=r(i)
rp1(k)=0
#endif

#if 1
if (r(3)==0) then ! either r(1) or r(2) <> 0
   rp1(1)=-r(2)
   rp1(2)=r(1)
   rp1(3)=0
else 
   rp1(1)=0
   rp1(2)=-r(3)
   rp1(3)=r(2)
endif
#endif



rp1 = rp1/ sqrt(rp1(1)**2+rp1(2)**2+rp1(3)**2)


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
      call mpi_reduce(dwork2,Dl(1,1,p),ndelta*ndir,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)

      dwork3=Dt(:,:,:,p)
      call mpi_reduce(dwork3,Dt(1,1,1,p),ndelta*ndir*2,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   enddo

   dwork2=H_ltt
   call mpi_reduce(dwork2,H_ltt,ndelta*ndir,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   dwork2=H_tt
   call mpi_reduce(dwork2,H_tt,ndelta*ndir,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)


   dwork3=D_ltt
   call mpi_reduce(dwork3,D_ltt,ndelta*ndir*2,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   dwork2=D_lltt
   call mpi_reduce(dwork2,D_lltt,ndelta*ndir,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)


   dwork3=SP_ltt
   call mpi_reduce(dwork3,SP_ltt,ndelta*ndir*2,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   dwork2=SP_lll
   call mpi_reduce(dwork2,SP_lll,ndelta*ndir,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   dwork3=SN_ltt
   call mpi_reduce(dwork3,SN_ltt,ndelta*ndir*2,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   dwork2=SN_lll
   call mpi_reduce(dwork2,SN_lll,ndelta*ndir,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)

   do p=1,3
      dwork2=w2s2(:,:,p)
      call mpi_reduce(dwork2,w2s2(1,1,p),ndelta*ndir,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
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



subroutine max_shear_coordinate_system(u_shear,idir_max,t1,t2)
use params

real*8 :: u_shear(3,3),u_strain(3,3)
real*8 :: t1(3),t2(3),t1t(3),t2t(3)
integer :: idir_max,i,j,k,l
real*8 :: testmax, maxval, norm, sp12,sp13,sp12max,sp13max

! local variables
integer :: idir
real*8 :: rhat(3),rperp1(3),rperp2(3) ,err,tmp(3)

maxval=0.0d0
u_strain = (u_shear + transpose(u_shear))/2

do idir=1,ndir
	
   rhat = dir(:,idir)
   rhat=rhat/sqrt(rhat(1)**2+rhat(2)**2+rhat(3)**2)
   call compute_perp(rhat,rperp1,rperp2)

   sp12=0   ! sp12 = r u_shear rperp1
   do i=1,3
      tmp(i) = sum(u_strain(i,:)*rperp1)
   enddo
   sp12 = sum(tmp*rhat)
   
   sp13=0   ! sp13 = r u_shear rperp2
   do i=1,3
      tmp(i) = sum(u_strain(i,:)*rperp2)
   enddo
   sp13 = sum(tmp*rhat)
   
   
   ! find the (rhat,rperp1,rperp2) which maximizes S_12^2 + S_13^2
   testmax = sqrt(sp13**2 + sp12**2)
   
   if (my_pe==io_pe) then
      write(6,*)'dir, testmax: ',idir,testmax
   endif
   if (testmax .gt. maxval) then
      maxval = testmax
      idir_max = idir
      t1 = rperp1
      t2 = rperp2
      sp12max=sp12
      sp13max=sp13
   endif
enddo
if (my_pe==io_pe) then
write(6,*)'u_strain = '
write(6,*)u_strain
write(6,*)'Idir_max= ',idir_max	
endif



! rotate so that in the (rmax,t1,t2) coordinate system, S_13=0, (see notes rotate_gradu.tex) 
! return the new t1 and t2 as expressed in the ORIGINAL coordinate system

rhat = dir(:,idir_max)
rhat = rhat/sqrt(rhat(1)**2+rhat(2)**2+rhat(3)**2)

if (maxval==0) then
   ! in case we run this on a zero initial condition - dont crash
   ! and leave t1,t2 as original
else

   norm = 1/(maxval)
   
   t1t = (sp12max*t1 + sp13max*t2)*norm
  
   t2t(1) = rhat(2)*t1t(3) - t1t(2)*rhat(3)
   t2t(2) = -rhat(1)*t1t(3) + rhat(3)*t1t(1)
   t2t(3) = rhat(1)*t1t(2) - t1t(1)*rhat(2)

   t1 = t1t
   t2 = t2t


#if 1
   ! recompute sp12, sp13, using new t1,t2.  
   ! make sure sp12=maxval, sp13=0

   sp12=0   ! sp12 = r u_shear rperp1
   do i=1,3
      tmp(i) = sum(u_strain(i,:)*t1)
   enddo
   sp12 = sum(tmp*rhat)
   sp13=0   ! sp12 = r u_shear rperp1
   do i=1,3
      tmp(i) = sum(u_strain(i,:)*t2)
   enddo
   sp13 = sum(tmp*rhat)
   

   err = (sp12-maxval)/maxval
   if(my_pe==io_pe) then
      if (abs(err)>1e-30 ) then
         print *,'Error r S t1 /= maxval:  ',sp12,maxval, abs(err)
      endif
   endif
   
   if(my_pe==io_pe) then
      err = sp13/maxval
      if (abs(err)>1e-30 ) then
         print *,'Error r S t2 should be 0:  ',sp13, abs(err)
      endif
   endif
#endif
endif
end subroutine


end module 



