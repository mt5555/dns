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
character(len=280) :: fname

if (restart==1) then

   ! initialize some constants, if needed on restart runs:

   ! set grav, fcor
   if (init_cond==3) call init_data_sht(Q,Qhat,work1,work2,0)   

   ! set xscale, yscale, offset grid 
   if (init_cond==4) call init_data_vxpair(Q,Qhat,work1,work2,0)   

   !
   ! read input uvw data from restart files.  sets 'time_initial'
   ! from restart data
   !
   call init_data_restart(Q,Qhat,work1,work2)

   ! set boundary data using input data.
   if (init_cond==4) call init_data_vxpair(Q,Qhat,work1,work2,2)  

   ! rescale Energy spectrum
   if (init_cond==9) then
      call init_data_decay(Q,Qhat,work1,work2,2,0,0)
      ! output rescaled data, but with different name:
      fname=runname(1:len_trim(runname)) // '-rescale'
      call output_uvw(fname,time_initial,Q,Qhat,work1,work2)		  
   endif

   if (npassive>0) then
   if (compute_passive_on_restart) then
      call init_passive_scalars(1,Q,Qhat,work1,work2)
      ! restart runs wont output data at t=0, for scalar output:
      call output_passive(runname,time_initial,Q,work1,work2)	
   else
      call init_passive_scalars(0,Q,Qhat,work1,work2)
      call input_passive(runname,time_initial,Q,work1,work2)
   endif
   endif

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
      call init_passive_scalars(1,Q,Qhat,work1,work2)
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










subroutine init_passive_scalars(init,Q,Qhat,work1,work2)
!
! low wave number, quasi isotropic initial condition
!
! init=0    initialize schmidt number only
! init=1    initialize schmidt number and passive scalar data
!
!
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(nx,ny,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: schmidt_table(5),xfac,mn,mx
integer :: n,k,i,j,im,jm,km,init,count,iter
character(len=80) :: message

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


if (init==0) return

count=0
do n=np1,np2
   count=count+1
   if (count<=2) then
      if (passive_type(n)==0) call passive_gaussian_init(Q,work1,work2,n)
      if (passive_type(n)==1) call passive_KE_init(Q,work1,work2,n)


      call global_min(Q(1,1,1,n),mn)
      call global_max(Q(1,1,1,n),mx)
      
      write(message,'(a,2f17.5)') 'initial passive scalar min/max: ',mn,mx
      call print_message(message)	
      call print_message("smothing passive scalar...")
      
      ! filter
      call fft3d(Q(1,1,1,n),work1)
      call fft_filter_dealias(Q(1,1,1,n))

      work2=Q(:,:,:,n)
      call ifft3d(work2,work1)
      call global_min(work2,mn)
      call global_max(work2,mx)
      
      ! laplacian smoothing:
      ! du/dt  =  laplacian(u)
      !  u_k = ((1 - k**2))**p  u_k    p = number of timesteps
      ! iterate until max and mn are within 5% of 0 and 1.  max iter=40
      iter=0
      do while ((mx>1.05 .or. mn<-.05) .and. (count<40) )
         iter=iter+1
         do k=nz1,nz2
            km=abs(kmcord(k))
            do j=ny1,ny2
               jm=abs(jmcord(j))
               do i=nx1,nx2
                  im=abs(imcord(i))
                  
                  xfac=(im*im+jm*jm+km*km)
                  xfac=xfac/(.25*g_nx*g_nx + .25*g_ny*g_ny + .25*g_nz*g_nz)
                  xfac=(1-.40*xfac)
                  
                  Q(i,j,k,n)=Q(i,j,k,n)*xfac
               enddo
            enddo
         enddo
         if (iter>=10) then  ! start checking 'mx'
            work2=Q(:,:,:,n)
            call ifft3d(work2,work1)
            call global_min(work2,mn)
            call global_max(work2,mx)
         endif
      enddo
      call ifft3d(Q(1,1,1,n),work1)
      write(message,'(a,2f17.5,a,i3)') 'after smoothing: min/max: ',mn,mx,'  iter=',iter
      call print_message(message)	
   else
      if (mod(count,2)==1) then
         Q(:,:,:,n)=Q(:,:,:,np1) 
         call print_message('Re-using i.c. from passive scalar 1')
      else
         Q(:,:,:,n)=Q(:,:,:,np1+1)          
         call print_message('Re-using i.c. from passive scalar 2')
      endif
   endif
enddo


end subroutine








subroutine passive_KE_init(Q,work1,work2,np)
!
! low wave number, quasi isotropic initial condition
!
use params
use mpi
implicit none
integer :: np
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

real*8 :: ke2,ke,ke_thresh,check
integer :: i,j,k,n,ierr

character(len=80) ::  message

write(message,'(a,i3,a)') "Initializing KE correlated double delta passive scalar n=",np
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



ke_thresh=.77*ke  ! .82 has too much 0, not enough 1 at 256^3
                  ! .77 has a touch too much at 1, not enought at 0

Q(:,:,:,np)=0
check=0
do k=nz1,nz2
   do j=ny1,ny2
      do i=nx1,nx2
         ke = .5*(Q(i,j,k,1)**2+Q(i,j,k,2)**2+Q(i,j,k,3)**2)
         if (ke>ke_thresh) then
            Q(i,j,k,np)=1
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













subroutine passive_gaussian_init(Q,work1,work2,np)
!
! low wave number, quasi isotropic initial condition
!
use params
use transpose
implicit none
integer :: np
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

real*8 :: ke2,ke,ke_thresh,check
integer :: i,j,k,n,ierr,NUMBANDS
real*8,allocatable :: enerb_target(:)
real*8,allocatable :: enerb(:)
real*8 :: ener
CPOINTER :: null
character(len=80) ::  message

write(message,'(a,i3,a)') "Initializing double delta passive scalars n=",np," ..."
call print_message(message)


NUMBANDS=.5 + sqrt(2.0)*g_nmin/3  ! round up
allocate(enerb_target(NUMBANDS))
allocate(enerb(NUMBANDS))

call livescu_spectrum(enerb_target,NUMBANDS)

call input1(Q(1,1,1,np),work1,work2,null,io_pe,.true.,-1)  
call rescale_e(Q(1,1,1,np),work1,ener,enerb,enerb_target,NUMBANDS,1)
! convert to 0,1:
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   if (Q(i,j,k,np)<0) then
      Q(i,j,k,np)=0
   else
      Q(i,j,k,np)=1
   endif
enddo
enddo
enddo

deallocate(enerb_target)
deallocate(enerb)
end subroutine










subroutine ranvor(Q,PSI,work,work2,rantype)
use params
use transpose
implicit none
real*8 :: Q(nx,ny,nz,3)
real*8 :: PSI(nx,ny,nz,3)
real*8 :: work(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
integer :: rantype,n,im,jm,km,i,j,k
real*8 :: xfac,alpha,beta
character(len=80) :: message
CPOINTER :: null


call print_message("computing random initial vorticity")
!random vorticity
if (rantype==0) then
do n=1,3
   ! input from random number generator
   ! this gives same I.C independent of cpus
   call input1(PSI(1,1,1,n),work2,work,null,io_pe,.true.,-1)  
enddo
else if (rantype==1) then
   do n=1,3
   write(message,*) 'random initial vorticity n=',n   
   call print_message(message)
   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         call gaussian( PSI(nx1,j,k,n), nslabx  )
         do i=nx1,nx2
            im=imcord(i)
            xfac = (2*2*2)
            if (km==0) xfac=xfac/2
            if (jm==0) xfac=xfac/2
            if (im==0) xfac=xfac/2
            PSI(i,j,k,n)=PSI(i,j,k,n)/xfac
         enddo
      enddo
   enddo
   enddo
else
   call abort("decay():  invalid 'rantype'")
endif
alpha=0
beta=1
do n=1,3
   write(message,*) 'solving for PSI n=',n   
   call print_message(message)
   call helmholtz_periodic_inv(PSI(1,1,1,n),work,alpha,beta)
enddo

call print_message("computing curl PSI")
if (ndim==2) then
   ! 2D case, treat PSI(:,:,:,1) = stream function, ignore other components
   Q=0
   ! u = PSI_y
   call der(PSI,Q,work2,work,DX_ONLY,2)     
   ! v = -PSI_x
   call der(PSI,Q(1,1,1,2),work2,work,DX_ONLY,1)
   Q(:,:,:,2)=-Q(:,:,:,2)
else
   call vorticity(Q,PSI,work,work2)
endif
call print_message("random initial condition complete")

end subroutine






subroutine rescale_e(Q,work,ener,enerb,enerb_target,NUMBANDS,nvec)
!
! rescale initial data
!
use params
use transpose
use mpi
implicit none
integer :: NUMBANDS,nvec
real*8 :: Q(nx,ny,nz,nvec)
real*8 :: work(nx,ny,nz)
real*8 :: enerb_target(NUMBANDS)
real*8 :: enerb(NUMBANDS)
real*8 :: enerb_work(NUMBANDS)
real*8 :: ener

integer :: n,im,jm,km,i,j,k,nb,ierr
real*8 :: xfac,xw
character(len=80) :: message


call print_message("Rescaling initial condition")
enerb=0
do n=1,nvec
   write(message,*) 'FFT to spectral space n=',n   
   call print_message(message)
   call fft3d(Q(1,1,1,n),work) 
   write(message,*) 'computing E(k) n=',n   
   call print_message(message)
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

            do nb=1,NUMBANDS
               if (xw>=nb-.5 .and. xw<nb+.5) &
                    enerb(nb)=enerb(nb)+.5*xfac*Q(i,j,k,n)**2
            enddo
            !remaining coefficients to 0:
            nb=NUMBANDS+1
            if (xw>=nb-.5) Q(i,j,k,n)=0

         enddo
      enddo
   enddo
enddo

#ifdef USE_MPI
   enerb_work=enerb
   call MPI_allreduce(enerb_work,enerb,NUMBANDS,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif


do n=1,nvec
   write(message,*) 'normalizing to E_target(k) n=',n   
   call print_message(message)

   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)
            xw=sqrt(real(km**2+jm**2+im**2))
            
            do nb=1,NUMBANDS
            if (xw>=nb-.5 .and. xw<nb+.5) then
               Q(i,j,k,n)=Q(i,j,k,n)*sqrt(enerb_target(nb)/(enerb(nb)))
            endif
            enddo
         enddo
      enddo
   enddo
enddo

ener=0
enerb=0
do n=1,nvec
   write(message,*) 're-computing E(k) n=',n   
   call print_message(message)
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

            do nb=1,NUMBANDS
            if (xw>=nb-.5 .and. xw<nb+.5) then
               enerb(nb)=enerb(nb)+.5*xfac*Q(i,j,k,n)**2
            endif
            enddo
            ener=ener+.5*xfac*Q(i,j,k,n)**2


         enddo
      enddo
   enddo
enddo

do n=1,nvec
   write(message,*) 'FFT back to grid space, n=',n   
   call print_message(message)
   call ifft3d(Q(1,1,1,n),work) 
enddo
#ifdef USE_MPI
   xfac=ener
   call MPI_allreduce(xfac,ener,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   enerb_work=enerb
   call MPI_allreduce(enerb_work,enerb,NUMBANDS,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
end subroutine
