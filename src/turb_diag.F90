#include "macros.h"
subroutine output_model(doit_model,time,Q,Qhat,q1,q2,q3,work1,work2)
use params
use pdf
use spectrum
use isoave
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(g_nz2,nslabx,ny_2dz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: q2(nx,ny,nz,n_var)
real*8 :: q3(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time
logical :: doit_model

! local variables
integer,parameter :: nints_e=28,npints_e=23	
real*8 :: ints_e(nints_e)
real*8 :: pints_e(npints_e,n_var)
real*8 :: x
integer i,j,k,n,ierr,csig
character(len=80) :: message
CPOINTER fid,fidj,fidS




if (compute_transfer) then
   compute_transfer=.false.
   ! spec_r computed last time step
   ! spec_diff, spec_f, spec_rhs were computed in RHS computation at the
   ! beginning of this flag (becuase compute_transfer flag was set)
   ! So they are all known at time_old. Now compute spec_r_new 
   ! (used to compute edot_r)
   call compute_Edotspec(time,Q,q1,work1,work2)
   ! output all the spectrum:
   call output_tran(time,Q,q1,q2,q3,work1,work2)
endif


! compute spectrum
! always compute at first timestep because transfer cannot be computed
! on last timestep.   
if (doit_model .or. time==time_initial) then
if ( g_bdy_x1==PERIODIC .and. &
     g_bdy_y1==PERIODIC .and. &
     g_bdy_z1==PERIODIC) then
   call compute_spec(time,Q,q1,work1,work2)
   call output_spec(time,Q,q1,q2,q3,work1,work2)
   
   !set this flag so that for next timestep, we will compute and save
   !spectral transfer functions:
   compute_transfer=.true.
endif
endif


if (.not.doit_model) return

!
! the "expensive" scalars
!

   call compute_expensive_scalars(Q,q1,q2,q3,work1,work2,ints_e,nints_e)
   if (npassive>0) then
      ! copy data computed above so that q3 = (ux,vy,wz)
      q3(:,:,:,1)=q1(:,:,:,1)
      q3(:,:,:,2)=q2(:,:,:,2)
      !q3(:,:,:,3)=q3(:,:,:,3)
      do n=np1,np2
         call compute_expensive_pscalars(Q,n,q1,q2,q3,work1,&
              pints_e(1,n),npints_e)
      enddo
   endif


   if (minval(ints_e(1:3))>0) then
      write(message,'(a,3f14.8)') 'skewness ux,vw,wz: ',&
           (ints_e(n+3)/ints_e(n)**1.5,n=1,3)
      call print_message(message)
      
      write(message,'(a,f14.8)') 'wSw: ',&
           ints_e(10)/ ( (ints_e(1)**2 + ints_e(2)**2 + ints_e(3)**2)/3 )
      call print_message(message)
   endif


   ! output turb scalars
   if (my_pe==io_pe) then
      write(message,'(f10.4)') 10000.0000 + time
      message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".scalars-turb"
      call copen(message,"w",fid,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "diag_output(): Error opening .scalars-turb file errno=",ierr
         call abort(message)
      endif
      x=nints_e; call cwrite8(fid,x,1)
      call cwrite8(fid,time,1)
      call cwrite8(fid,ints_e,nints_e)
      call cclose(fid,ierr)
   endif


   ! output turb passive scalars data
   if (my_pe==io_pe .and. npassive>0) then
      write(message,'(f10.4)') 10000.0000 + time
      message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".pscalars-turb"
      call copen(message,"w",fid,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "diag_output(): Error opening .pscalars-turb file errno=",ierr
         call abort(message)
      endif
      x=nints_e; call cwrite8(fid,x,1)
      x=npassive; call cwrite8(fid,x,1)
      call cwrite8(fid,time,1)
      x=mu; call cwrite8(fid,x,1)
      do n=np1,np2	
         call cwrite8(fid,pints_e(1,n),nints_e)
      enddo	 
      
      call cclose(fid,ierr)
   endif



!
! output structure functions and time averaged forcing
! 
if (diag_struct==1) then

   ! angle averaged functions:
   call isoavep(Q,q1,q2,q3,3,csig)
   ! if csig>0, isoavep did not complete - interrupted by SIGURG
   if (my_pe==io_pe .and. csig==0) then
      write(message,'(f10.4)') 10000.0000 + time
      message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".isostr"
      call copen(message,"w",fid,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "output_model(): Error opening .isostr file errno=",ierr
         call abort(message)
      endif
      call writeisoave(fid,time)
      call cclose(fid,ierr)
   endif
endif


if (diag_pdfs==1) then
   ! uscale=.01 / pi2
   ! epsscale=.01 / pi2_squared  
   call compute_all_pdfs(Q,q1)
   
   
   
   if (structf_init==1) then
   if (my_pe==io_pe) then
      write(message,'(f10.4)') 10000.0000 + time
      message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".sf"
      call copen(message,"w",fid,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "output_model(): Error opening .sf file errno=",ierr
         call abort(message)
      endif

      if (compute_uvw_jpdfs) then
      write(message,'(f10.4)') 10000.0000 + time
      message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".jpdf"
      call copen(message,"w",fidj,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "output_model(): Error opening .jpdf file errno=",ierr
         call abort(message)
      endif
      endif

      if (compute_passive_pdfs) then
      write(message,'(f10.4)') 10000.0000 + time
      message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".spdf"
      call copen(message,"w",fidS,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "output_model(): Error opening .spdf file errno=",ierr
         call abort(message)
      endif
      endif

   endif
   call output_pdf(time,fid,fidj,fidS)
   if (my_pe==io_pe) call cclose(fid,ierr)
   if (compute_uvw_jpdfs .and. my_pe==io_pe) call cclose(fidj,ierr)
   if (compute_passive_pdfs .and. my_pe==io_pe) call cclose(fidS,ierr)
   endif
endif

! time averaged dissapation and forcing:
!call compute_time_averages(Q,q1,q2,q3(1,1,1,1),q3(1,1,1,2),q3(1,1,1,3),time)

end subroutine




#if 0
subroutine compute_time_averages(Q,Qhat,f,wsum,work1,dxx,time)
use params
use sforcing
use fft_interface
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: wsum(nx,ny,nz)
real*8 :: work1(nx,ny,nz)
real*8 :: dxx(nx,ny,nz)
real*8 :: Qhat(g_nz2,nslabx,ny_2dz,n_var)
real*8 :: f(g_nz2,nslabx,ny_2dz,n_var)
real*8 :: time

! local
integer :: i,j,k,n,n1,n2,ierr
real*8,save,allocatable :: diss(:,:,:)
real*8,save,allocatable :: diss2(:,:,:)
real*8,save,allocatable :: uf(:,:,:)
real*8,save,allocatable :: uf2(:,:,:)
integer,save :: ntave=0
real*8 :: f_diss,x
character(len=80) message
character(len=280) fname



if (ntave==0) then
   allocate(diss(nx,ny,nz))
   allocate(diss2(nx,ny,nz))
   allocate(uf(nx,ny,nz))
   allocate(uf2(nx,ny,nz))
   diss=0
   diss2=0
   uf=0
   uf2=0
endif
ntave=ntave+1

wsum=0
do n1=1,3
do n2=1,3
   ! Q(:,:,:,n1)* d(Q)/dn2(:,:,:,n1)
   call der(Q(1,1,1,n1),f,dxx,work1,DX_AND_DXX,n2)
   wsum=wsum+mu*Q(:,:,:,n1)*dxx(:,:,:)
enddo
enddo
diss=(diss*(ntave-1) + wsum) / ntave
diss2=(diss2*(ntave-1) + wsum**2) / ntave


do n=1,3
   wsum=Q(:,:,:,n)
   call z_fft3d_trashinput(wsum,Qhat(1,1,1,n),work1)
enddo
f=0
call sforce(f,Qhat,f_diss)
wsum=0
do n=1,3
   call z_ifft3d(f(1,1,1,n),dxx,work1)
   wsum=wsum+dxx(:,:,:)*Q(:,:,:,n)
enddo
uf=(uf*(ntave-1) + wsum) / ntave
uf2=(uf2*(ntave-1) + wsum**2) / ntave



if (time>=time_final) then
   ! time to save the output
   write(message,'(f10.4)') 10000.0000 + time_initial
   fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".diss"
   x=ntave
   call singlefile_io(x,diss,fname,work1,dxx,0,io_pe)

   write(message,'(f10.4)') 10000.0000 + time_initial
   fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".diss2"
   x=ntave
   call singlefile_io(x,diss2,fname,work1,dxx,0,io_pe)

   write(message,'(f10.4)') 10000.0000 + time_initial
   fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".uf"
   x=ntave
   call singlefile_io(x,uf,fname,work1,dxx,0,io_pe)

   write(message,'(f10.4)') 10000.0000 + time_initial
   fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".uf2"
   x=ntave
   call singlefile_io(x,uf2,fname,work1,dxx,0,io_pe)


endif


end subroutine
#endif




subroutine iso_stats(Q,PSI,work,work2)
!
! compute random isotropic grid data, take the FFT and then
! compute stats for the coefficients.
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
integer km,jm,im,i,j,k,n,wn,ierr,nb
integer,allocatable :: seed(:)
integer,parameter :: NUMBANDS=100
real*8 xw,enerb(NUMBANDS),enerb_target(NUMBANDS),ener,xfac,theta
real*8 enerb_work(NUMBANDS)
character(len=80) message

integer :: trys,num_trys=100
integer,parameter :: nbin=100
real*8 :: bindel=.2e-5
real*8 :: count(-nbin:nbin,-3:3,-3:3,-3:3)

count=0

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

if (ncpus>1) then
   call abort("iso stats requires only 1 cpu in cartesian communicator")
endif

do trys=1,num_trys

do n=1,3
do j=ny1,ny2
do k=nz1,nz2
   call gaussian(PSI(nx1,j,k,n),nx2-nx1+1)
enddo
enddo
enddo

alpha=0
beta=1
do n=1,3
   call helmholtz_periodic_inv(PSI(1,1,1,n),work,alpha,beta)
enddo
call vorticity(Q,PSI,work,work2)


enerb=0
do n=1,3
   call fft3d(Q(1,1,1,n),work) 
   do k=nz1,nz2
      km=kmcord(k)
      do j=ny1,ny2
         jm=jmcord(j)
         do i=nx1,nx2
            im=imcord(i)
            xw=sqrt(real(km**2+jm**2+im**2))
            if (xw<3.5) then
               nb=nint(Q(i,j,k,n)/bindel)
               if (nb>100) nb=100
               if (nb<-100) nb=-100
               count(nb,im,jm,km)=count(nb,im,jm,km)+1

               xfac = (2*2*2)
               if (km==0) xfac=xfac/2
               if (jm==0) xfac=xfac/2
               if (im==0) xfac=xfac/2
               
               nb=nint(xw)
               enerb(nb)=enerb(nb)+.5*xfac*Q(i,j,k,n)**2
            endif
         enddo
      enddo
   enddo
enddo

print *,'try=',trys,num_trys
do nb=1,3
   print *,nb,enerb(nb)
enddo

enddo

do i=-nbin,nbin
   write(*,'(i10,3f10.5)') i,(count(i,1,1,1),n=1,3)
enddo

stop
end subroutine




















subroutine compute_expensive_scalars(Q,gradu,gradv,gradw,work,work2,scalars,ns)
!
!
use params
use fft_interface
use transpose
implicit none
integer :: ns
real*8 :: scalars(ns)
real*8 Q(nx,ny,nz,n_var)    
real*8 work(nx,ny,nz)
real*8 work2(nx,ny,nz)
real*8 gradu(nx,ny,nz,n_var)    
real*8 gradv(nx,ny,nz,n_var)    
real*8 gradw(nx,ny,nz,n_var)    

!local
real*8 :: scalars2(ns)
integer n1,n1d,n2,n2d,n3,n3d,ierr
integer i,j,k,n,m1,m2
real*8 :: vor(3),Sw(3),wS(3),Sww,ux2(3),ux3(3),ux4(3),uij,uji,u2(3)
real*8 :: vor2(3),vor3(3),vor4(3)
real*8 :: uxx2(3)
real*8 :: dummy(1),S2sum,ensave,S4sum,S2,S4,S2w2
real*8 :: tmx1,tmx2,xtmp

!
! compute derivatives
!
uxx2=0
do n=1,3
   if (n==1) then
      call der(Q(1,1,1,1),gradu(1,1,1,n),work2,work,DX_AND_DXX,n)
   else
      call der(Q(1,1,1,1),gradu(1,1,1,n),dummy,work,DX_ONLY,n)
   endif

   if (n==2) then
      call der(Q(1,1,1,2),gradv(1,1,1,n),work2,work,DX_AND_DXX,n)
   else
      call der(Q(1,1,1,2),gradv(1,1,1,n),dummy,work,DX_ONLY,n)
   endif

   if (n==3) then
      call der(Q(1,1,1,3),gradw(1,1,1,n),work2,work,DX_AND_DXX,n)
   else
      call der(Q(1,1,1,3),gradw(1,1,1,n),dummy,work,DX_ONLY,n)
   endif

   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      uxx2(n)=uxx2(n)+work2(i,j,k)*2
   enddo
   enddo
   enddo   
enddo






! scalars
S2sum=0
S4sum=0
S2w2=0
ensave=0
Sww=0
ux2=0
ux3=0
ux4=0
u2=0
vor2=0
vor3=0
vor4=0

do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   do n=1,3
      u2(n)=u2(n)+Q(i,j,k,n)**2
   enddo

   vor(1)=gradw(i,j,k,2)-gradv(i,j,k,3)
   vor(2)=gradu(i,j,k,3)-gradw(i,j,k,1)
   vor(3)=gradv(i,j,k,1)-gradu(i,j,k,2)

   ! compute Sw = Sij*wj
   Sw=0
   !wS=0
   S2=0
   S4=0
   do m1=1,3
      do m2=1,3
         if (m1==1) uij=gradu(i,j,k,m2)
         if (m1==2) uij=gradv(i,j,k,m2)
         if (m1==3) uij=gradw(i,j,k,m2)
         if (m2==1) uji=gradu(i,j,k,m1)
         if (m2==2) uji=gradv(i,j,k,m1)
         if (m2==3) uji=gradw(i,j,k,m1)
         ! S(m1,m2) = .5*(uij_uji)
         Sw(m1)=Sw(m1)+.5*(uij+uji)*vor(m2)
         !wS(m2)=wS(m2)+.5*(uij+uji)*vor(m1)
         xtmp=(.5*(uij+uji))**2
         S2=S2 + xtmp
         S4=S4 + xtmp**2
      enddo
   enddo
   S2sum=S2sum+S2
   S4sum=S4sum+S4
   ! compute Sww = wi*(Sij*wj)
   Sww=Sww+Sw(1)*vor(1)+Sw(2)*vor(2)+Sw(3)*vor(3)

   xtmp=vor(1)**2+vor(2)**2+vor(3)**2
   ensave=ensave+xtmp
   S2w2 = S2*xtmp

   ! if we use gradu(i,j,k,1)**3, do we preserve the sign?  
   ! lets not put f90 to that test!
   uij=gradu(i,j,k,1)**2
   ux2(1)=ux2(1)+uij
   ux3(1)=ux3(1)+uij*gradu(i,j,k,1)
   ux4(1)=ux4(1)+uij*uij

   uij=gradv(i,j,k,2)**2
   ux2(2)=ux2(2)+uij
   ux3(2)=ux3(2)+uij*gradv(i,j,k,2)
   ux4(2)=ux4(2)+uij*uij

   uij=gradw(i,j,k,3)**2
   ux2(3)=ux2(3)+uij
   ux3(3)=ux3(3)+uij*gradw(i,j,k,3)
   ux4(3)=ux4(3)+uij*uij

   vor2=vor2 + vor**2
   vor3=vor3 + vor*vor**2  ! will **3 preserve sign?
   vor4=vor4 + vor**4
enddo
enddo
enddo

S2sum=S2sum/g_nx/g_ny/g_nz
S2w2=S2w2/g_nx/g_ny/g_nz
S4sum=S4sum/g_nx/g_ny/g_nz
Sww=Sww/g_nx/g_ny/g_nz
ux2=ux2/g_nx/g_ny/g_nz
ux3=ux3/g_nx/g_ny/g_nz
ux4=ux4/g_nx/g_ny/g_nz
vor2=vor2/g_nx/g_ny/g_nz
vor3=vor3/g_nx/g_ny/g_nz
vor4=vor4/g_nx/g_ny/g_nz
u2=u2/g_nx/g_ny/g_nz
uxx2=uxx2/g_nx/g_ny/g_nz

ensave=ensave/g_nx/g_ny/g_nz




ASSERT("compute_expensive_scalars: ns too small ",ns>=28)
do n=1,3
scalars(n)=ux2(n)
scalars(n+3)=ux3(n)
scalars(n+6)=ux4(n)
enddo
scalars(10)=Sww
do n=1,3
scalars(10+n)=u2(n)
enddo
scalars(14)=S2sum

scalars(15:17)=vor2
scalars(18:20)=vor3
scalars(21:23)=vor4

scalars(24)=S4sum
scalars(25)=S2w2
scalars(26:28)=uxx2

#ifdef USE_MPI
   scalars2=scalars
   call MPI_allreduce(scalars2,scalars,ns,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif


end subroutine










subroutine compute_expensive_pscalars(Q,np,grads,grads2,gradu,work,scalars,ns)
!
! INPUT:  Q,np, gradu = (ux,vy,wz)
! 
use params
use fft_interface
use transpose
implicit none
integer :: ns,np
real*8 :: scalars(ns)
real*8 Q(nx,ny,nz,n_var)    
real*8 work(nx,ny,nz)
real*8 grads(nx,ny,nz,n_var)    
real*8 grads2(nx,ny,nz,n_var)    
real*8 gradu(nx,ny,nz,n_var)  

!local
real*8 :: scalars2(ns)
integer n1,n1d,n2,n2d,n3,n3d,ierr
integer i,j,k,n,m1,m2
real*8 :: ux2(3),ux3(3),ux4(3),u2,uij,su(3)
real*8 :: uxx2(3),uxx3(3),uxx4(3)

!
!  gradu = ux,vy,wz,
!
! grads = d/dx, d/dy and d/dz
! grads2 = d/dxx, d/dyy and d/dzz
do n=1,3
   call der(Q(1,1,1,np),grads(1,1,1,n),grads2(1,1,1,n),work,DX_AND_DXX,n)
enddo


! scalars
ux2=0
ux3=0
ux4=0
u2=0
su=0
uxx2=0
uxx3=0
uxx4=0


do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
    u2=u2+Q(i,j,k,np)**2

   ! if we use grads(i,j,k,1)**3, do we preserve the sign?  
   ! lets not put f90 to that test!
   do n=1,3
      uij=grads(i,j,k,n)**2
      ux2(n)=ux2(n)+uij
      ux3(n)=ux3(n)+uij*grads(i,j,k,n)
      ux4(n)=ux4(n)+uij*uij

      uij=grads2(i,j,k,n)**2
      uxx2(n)=uxx2(n)+uij
      uxx3(n)=uxx3(n)+uij*grads2(i,j,k,n)
      uxx4(n)=uxx4(n)+uij*uij

      su(n) = su(n) + gradu(i,j,k,n)*grads(i,j,k,n)*grads(i,j,k,n)
   enddo

enddo
enddo
enddo

u2=u2/g_nx/g_ny/g_nz
ux2=ux2/g_nx/g_ny/g_nz
ux3=ux3/g_nx/g_ny/g_nz
ux4=ux4/g_nx/g_ny/g_nz
uxx2=uxx2/g_nx/g_ny/g_nz
uxx3=uxx3/g_nx/g_ny/g_nz
uxx4=uxx4/g_nx/g_ny/g_nz
su=su/g_nx/g_ny/g_nz





ASSERT("compute_expensive_pscalars: ns too small ",ns>=23)

scalars(1)=schmidt(np)
scalars(2)=u2
i=2

do n=1,3
scalars(n+i)=ux2(n)
scalars(n+3+i)=ux3(n)
scalars(n+6+i)=ux4(n)
enddo
i=i+9

do n=1,3
scalars(n+i)=uxx2(n)
scalars(n+3+i)=uxx3(n)
scalars(n+6+i)=uxx4(n)
enddo
i=i+9

do n=1,3
scalars(n+i)=su(n)
enddo
i=i+3

#ifdef USE_MPI
   scalars2=scalars
   call MPI_allreduce(scalars2,scalars,i,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif


end subroutine





