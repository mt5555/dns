#include "macros.h"
subroutine init_data(Q,Qhat,q1,work1,work2)
use params
use mpi
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
if (restart==1) then

   ! initialize some constants, if needed on restart runs:

   ! set grav, fcor
   if (init_cond==3) call init_data_sht(Q,Qhat,work1,work2,0)   

   ! set xscale, yscale, offset grid 
   if (init_cond==4) call init_data_vxpair(Q,Qhat,work1,work2,0)   

   !
   ! read input data from restart files
   !
   call init_data_restart(Q,Qhat,work1,work2)

   ! set boundary data using input data.
   if (init_cond==4) call init_data_vxpair(Q,Qhat,work1,work2,2)  

   ! rescale Energy spectrum
   if (init_cond==9) call init_data_decay(Q,Qhat,work1,work2,2,0,0)

else
   if (init_cond==0) call init_data_khblob(Q,Qhat,work1,work2)
   if (init_cond==1) call init_data_kh(Q,Qhat,work1,work2)
   if (init_cond==2) call init_data_lwisotropic(Q,Qhat,work1,work2,1,0)
   if (init_cond==3) call init_data_sht(Q,Qhat,work1,work2,1)
   if (init_cond==4) call init_data_vxpair(Q,Qhat,work1,work2,1)
   if (init_cond==5) call init_data_lwisotropic(Q,Qhat,work1,work2,1,1)
   if (init_cond==6) call init_data_zero(Q,Qhat,work1,work2)
   if (init_cond==7) call init_data_decay(Q,Qhat,work1,work2,1,0,0)
   if (init_cond==8) call init_data_decay(Q,Qhat,work1,work2,1,1,0)

   if (npassive>0) then
      call init_passive_scalars(Q,Qhat,work1,work2)
   endif

   if (equations==NS_UVW) then
      call print_message('Projecting initial data...')
      call divfree_gridspace(Q,work1,work2,q1) 
   else if (equations==SHALLOW .or. equations==NS_PSIVOR) then
      if (dealias>0)  then
         call print_message('Dealiasing initial data...')
         call dealias_gridspace(Q,work1)
      endif
   endif
endif
end subroutine









subroutine init_data_restart(Q,Qhat,work1,work2)
!
! low wave number, quasi isotropic initial condition
!
use params
use mpi
use fft_interface
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

!local
character(len=80) message
character(len=80) fname
integer :: n

Q=0
time_initial=-1
call input_uvw(time_initial,Q,Qhat,work1,work2,1)

write(message,'(a,f10.4)') "restart time=",time_initial
call print_message(message)
end subroutine










subroutine init_passive_scalars(Q,Qhat,work1,work2)
!
! low wave number, quasi isotropic initial condition
!
use params
use mpi
use fft_interface
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(nx,ny,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: schmidt_table(5)
integer :: n,k

schmidt=0
passive_type=0

schmidt_table(1)=.01
schmidt_table(2)=.05
schmidt_table(3)=.10
schmidt_table(4)=.5
schmidt_table(5)=1.0

k=0
call print_message("passive scalars:")
call print_message("    n   Schmidt   Type (0=Gaussian, 1=KE based)")
do n=np1,np2
   if (mod(n-np1,2)==0) then
      passive_type(n)=0   ! gaussian
      k=k+1
   else
      passive_type(n)=1   ! KE based
   endif
   if (k>5) then
      call abort("Error: passive_scalar_init(): schmidt_table too small.")
   endif
   schmidt(n)=schmidt_table(k)
   if (my_pe==io_pe) then
      write(*,'(i4,f11.3,i7)') n-np1+1,schmidt(n),passive_type(n)
   endif
enddo


do n=np1,np2
   if (passive_type(n)==0) call passive_gaussian_init(Q,work1,work2,n)
   if (passive_type(n)==1) call passive_KE_init(Q,work1,work2,n)
enddo


end subroutine








subroutine passive_KE_init(Q,work1,work2)
!
! low wave number, quasi isotropic initial condition
!
use params
use mpi
use fft_interface
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

real*8 :: ke2,ke,ke_thresh,check
integer :: i,j,k,n,ierr

character(len=80) ::  message

write(message,'(a,i3,a)') "Initializing KE correlated passive scalars n=",n,&
                        " ..."

call print_message(message)

ke=0
do n=1,ndim
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      ke = ke + .5*Q(i,j,k,n)**2
   enddo
   enddo
   enddo
enddo
ke=ke/g_nx/g_ny/g_nz

#ifdef USE_MPI
   ke2=ke
   call MPI_allreduce(ke2,ke,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif



ke_thresh=.82*ke

Q(:,:,:,n)=0
check=0
do k=nz1,nz2
   do j=ny1,ny2
      do i=nx1,nx2
         ke = .5*(Q(i,j,k,1)**2+Q(i,j,k,2)**2+Q(i,j,k,3)**2)
         if (ke>ke_thresh) then
            Q(:,:,:,n)=1
            check=check+1
         endif
      enddo
   enddo
enddo

check=check/g_nx/g_ny/g_nz
#ifdef USE_MPI
   ke2=check
   call MPI_allreduce(ke2,check,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif


write(message,'(a,f7.3)') "Mean density: ",check
call print_message(message)


end subroutine













subroutine passive_gaussian_init(Q,work1,work2)
!
! low wave number, quasi isotropic initial condition
!
use params
use mpi
use fft_interface
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

real*8 :: ke2,ke,ke_thresh,check
integer :: i,j,k,n,ierr

character(len=80) ::  message

write(message,'(a,i3,a)') "Initializing KE correlated passive scalars n=",n,&
                        " ..."

call print_message(message)

ke=0
do n=1,ndim
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      ke = ke + .5*Q(i,j,k,n)**2
   enddo
   enddo
   enddo
enddo
ke=ke/g_nx/g_ny/g_nz

#ifdef USE_MPI
   ke2=ke
   call MPI_allreduce(ke2,ke,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif



ke_thresh=.82*ke

Q(:,:,:,n)=0
check=0
do k=nz1,nz2
   do j=ny1,ny2
      do i=nx1,nx2
         ke = .5*(Q(i,j,k,1)**2+Q(i,j,k,2)**2+Q(i,j,k,3)**2)
         if (ke>ke_thresh) then
            Q(:,:,:,n)=1
            check=check+1
         endif
      enddo
   enddo
enddo

check=check/g_nx/g_ny/g_nz
#ifdef USE_MPI
   ke2=check
   call MPI_allreduce(ke2,check,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif


write(message,'(a,f7.3)') "Mean density: ",check
call print_message(message)


end subroutine











