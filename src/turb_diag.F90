#include "macros.h"
subroutine output_model(doit_model,doit_diag,time,Q,Qhat,q1,q2,q3,work1,work2)
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
logical :: doit_model,doit_diag

! local variables
integer,parameter :: nints_e=49,npints_e=51
real*8 :: ints_e(nints_e)
real*8 :: pints_e(npints_e,n_var)
real*8 :: x,zero_len
real*8 :: divx,divi
integer i,j,k,n,ierr,csig
integer :: n1,n1d,n2,n2d,n3,n3d
character(len=80) :: message
CPOINTER fid,fidj,fidS,fidcore




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
   call output_spec(time,time_initial)
   call output_helicity_spec(time,time_initial)  ! put all hel spec in same file

   !set this flag so that for next timestep, we will compute and save
   !spectral transfer functions:
   compute_transfer=.true.

   ! for incompressible equations, print divergence as diagnostic:
   if (equations==NS_UVW) then
      call compute_div(Q,q1,work1,work2,divx,divi)
      write(message,'(3(a,e12.5))') 'max(div)=',divx
      call print_message(message)	
   endif


endif
endif

! do PDF's and scalars if doit_model=.true., OR if this is a restart
! but we have computed new passive scalars.
if ((compute_passive_on_restart .and. time==time_initial) .or. &
    doit_model) then
   ! do the rest of this suburoutine
else
   return
endif




!
! the "expensive" scalars
!
   call print_message("Computing velocity stats...")
   call compute_expensive_scalars(Q,q1,q2,q3,work1,work2,ints_e,nints_e)
   if (minval(ints_e(1:3))>0) then
      write(message,'(a,3f14.8)') 'skewness ux,vw,wz: ',&
           (ints_e(n+3)/ints_e(n)**1.5,n=1,3)
      call print_message(message)
      
      write(message,'(a,f14.8)') 'wSw: ',&
           ints_e(10)/ ( (ints_e(1)**2 + ints_e(2)**2 + ints_e(3)**2)/3 )
      call print_message(message)
   endif
   if (npassive>0) then
      call print_message("Computing passive scalar stats...")
      ! copy data computed above so that q3 = (ux,vy,wz)
      q3(:,:,:,1)=q1(:,:,:,1)
      q3(:,:,:,2)=q2(:,:,:,2)
      !q3(:,:,:,3)=q3(:,:,:,3)
      do n=np1,np2
         call compute_expensive_pscalars(Q,n,q1,q2,q3,work1,&
              pints_e(1,n),npints_e)

      enddo
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
      x=npints_e; call cwrite8(fid,x,1)
      x=npassive; call cwrite8(fid,x,1)
      call cwrite8(fid,time,1)
      x=mu; call cwrite8(fid,x,1)
      do n=np1,np2	
         call cwrite8(fid,pints_e(1,n),npints_e)
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
   call compute_all_pdfs(Q,q1)

   if (my_pe==io_pe) then
      write(message,'(f10.4)') 10000.0000 + time
      message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".sf"
      call copen(message,"w",fid,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "output_model(): Error opening .sf file errno=",ierr
         call abort(message)
      endif

      write(message,'(f10.4)') 10000.0000 + time
      message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".cores"
      call copen(message,"w",fidcore,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "output_model(): Error opening .cores file errno=",ierr
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
   call output_pdf(time,fid,fidj,fidS,fidcore)
   if (my_pe==io_pe) call cclose(fid,ierr)
   if (my_pe==io_pe) call cclose(fidcore,ierr)
   if (compute_uvw_jpdfs .and. my_pe==io_pe) call cclose(fidj,ierr)
   if (compute_passive_pdfs .and. my_pe==io_pe) call cclose(fidS,ierr)
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
real*8 :: vor(3),Sw(3),wS(3),Sww,ux2(3),ux3(3),ux4(3),uij,uji
real*8 :: u1(3),u2(3),u3(3),u4(3),u1tmp(3)
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
      uxx2(n)=uxx2(n)+work2(i,j,k)**2
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
u1=0
vor2=0
vor3=0
vor4=0

do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   do n=1,3
      u1(n)=u1(n)+Q(i,j,k,n)
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
uxx2=uxx2/g_nx/g_ny/g_nz
u1=u1/g_nx/g_ny/g_nz
ensave=ensave/g_nx/g_ny/g_nz


#ifdef USE_MPI
   call MPI_allreduce(u1,u1tmp,3,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#else
   u1tmp=u1
#endif
u2=0
u3=0
u4=0
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   do n=1,3
      xtmp=(Q(i,j,k,n)-u1tmp(n))**2
      u2(n)=u2(n)+xtmp
      u3(n)=u3(n)+xtmp*(Q(i,j,k,n)-u1tmp(n))
      u4(n)=u4(n)+xtmp*xtmp
   enddo
enddo
enddo
enddo

u2=u2/g_nx/g_ny/g_nz
u3=u3/g_nx/g_ny/g_nz
u4=u4/g_nx/g_ny/g_nz


ASSERT("compute_expensive_scalars: ns too small ",ns>=49)
scalars=0
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
scalars(29:31)=u1
scalars(32:34)=u3
scalars(35:37)=u4

#ifdef USE_MPI
   scalars2=scalars
   call MPI_allreduce(scalars2,scalars,ns,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif



! zero crossing:  (mean is in u1)
! long. only.  
call transpose_to_x(Q(1,1,1,1),work,n1,n1d,n2,n2d,n3,n3d)
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,u1(1),scalars(38))
call transpose_to_y(Q(1,1,1,2),work,n1,n1d,n2,n2d,n3,n3d)
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,u1(2),scalars(39))
call transpose_to_z(Q(1,1,1,3),work,n1,n1d,n2,n2d,n3,n3d)
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,u1(3),scalars(40))


i=40
u1=0  ! it's a deriavitve, so we know mean=0
! zero crossing of dissipation: 
call transpose_to_x(gradu(1,1,1,1),work,n1,n1d,n2,n2d,n3,n3d)
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,u1,scalars(i+1))
call transpose_to_y(gradv(1,1,1,2),work,n1,n1d,n2,n2d,n3,n3d)
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,u1,scalars(i+2))
call transpose_to_z(gradw(1,1,1,3),work,n1,n1d,n2,n2d,n3,n3d)
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,u1,scalars(i+3))
i=i+3
! i=43



end subroutine










subroutine compute_expensive_pscalars(Q,np,grads,grads2,gradu,work,scalars,ns)
!
! INPUT:  Q,np, gradu = (ux,vy,wz)
! 
use params
use fft_interface
use transpose
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
integer ii,i,j,k,n,m1,m2
real*8 :: ux1,ux2(3),ux3(3),ux4(3),u2,x1,x2,su(3),u1,u3,u4,u1tmp
real*8 :: uxx2(3),uxx3(3),uxx4(3),xtmp,u2ux2(3)
real*8 :: lnux1,lnux2,lnux3,lnux4,xln

!
! gradu = ux,vy,wz,
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
u1=0
u2=0
u3=0
u4=0
su=0
uxx2=0
uxx3=0
uxx4=0
lnux1=0
lnux2=0  
lnux3=0
lnux4=0

do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
    u1=u1+Q(i,j,k,np)

   ! if we use grads(i,j,k,1)**3, do we preserve the sign?  
   ! lets not put f90 to that test!
   do n=1,3
      x1=grads(i,j,k,n)**2
      ux2(n)=ux2(n)+x1
      ux3(n)=ux3(n)+x1*grads(i,j,k,n)
      ux4(n)=ux4(n)+x1*x1

      ! average over all three directions, so devide by 3 here:
      if (x1>0) then
         xln=log(x1)
      else
         xln=-100
      endif
      lnux1=lnux1 + xln/3
      lnux2=lnux2 + xln*xln/3
      lnux3=lnux3 + xln*xln*xln/3
      lnux4=lnux4 + xln*xln*xln*xln/3

      x2=grads2(i,j,k,n)**2
      uxx2(n)=uxx2(n)+x2
      uxx3(n)=uxx3(n)+x2*grads2(i,j,k,n)
      uxx4(n)=uxx4(n)+x2*x2

      su(n) = su(n) + gradu(i,j,k,n)*x1
   enddo

enddo
enddo
enddo

! <(a-a0)**2> = < a**2 > - a0**2
u1=u1/g_nx/g_ny/g_nz
ux2=ux2/g_nx/g_ny/g_nz
ux3=ux3/g_nx/g_ny/g_nz
ux4=ux4/g_nx/g_ny/g_nz
uxx2=uxx2/g_nx/g_ny/g_nz
uxx3=uxx3/g_nx/g_ny/g_nz
uxx4=uxx4/g_nx/g_ny/g_nz
su=su/g_nx/g_ny/g_nz

lnux1=lnux1/g_nx/g_ny/g_nz
lnux2=lnux2/g_nx/g_ny/g_nz
lnux3=lnux3/g_nx/g_ny/g_nz
lnux4=lnux4/g_nx/g_ny/g_nz


#ifdef USE_MPI
   call MPI_allreduce(u1,u1tmp,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#else
   u1tmp=u1
#endif


u2ux2=0
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
    xtmp=(Q(i,j,k,np)-u1tmp)**2
    u2=u2+xtmp
    u3=u3+xtmp*(Q(i,j,k,np)-u1tmp)
    u4=u4+xtmp*xtmp
    do n=1,3
       u2ux2(n) = u2ux2(n) + xtmp*grads(i,j,k,n)**2
    enddo
enddo
enddo
enddo
u2=u2/g_nx/g_ny/g_nz; 
u3=u3/g_nx/g_ny/g_nz; 
u4=u4/g_nx/g_ny/g_nz; 
u2ux2=u2ux2/g_nx/g_ny/g_nz; 




! we will sum over all pe's below, so do this for non-sums:
scalars=0
scalars(1)=0
if (my_pe==io_pe) scalars(1)=schmidt(np)

scalars(2)=u2
i=2

do n=1,3
scalars(n+i)=ux2(n)         ! 3,4,5
scalars(n+3+i)=ux3(n)       ! 6,7,8
scalars(n+6+i)=ux4(n)       ! 9,10,11
enddo
i=i+9

! i = 11
do n=1,3
scalars(n+i)=uxx2(n)       ! 12,13,14
scalars(n+3+i)=uxx3(n)     ! 15,16,17
scalars(n+6+i)=uxx4(n)     ! 18,19,20
enddo
i=i+9

! i=20
do n=1,3
scalars(n+i)=su(n)         ! 21,22,23
enddo
i=i+3

i=i+1
scalars(i)=u3              ! 24
i=i+1
scalars(i)=u4              ! 25

scalars(26)=0               ! opps, forgot this one


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! zero crossing:  (mean is in pints_e(24,n))
! other values:
ux1=.8; ux2=.9
call transpose_to_x(Q(1,1,1,np),work,n1,n1d,n2,n2d,n3,n3d)
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,scalars(24),scalars(27))
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,ux1,scalars(28))
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,ux2,scalars(29))

call transpose_to_y(Q(1,1,1,np),work,n1,n1d,n2,n2d,n3,n3d)
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,scalars(24),scalars(30))
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,ux1,scalars(31))
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,ux2,scalars(32))

call transpose_to_z(Q(1,1,1,np),work,n1,n1d,n2,n2d,n3,n3d)
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,scalars(24),scalars(33))
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,ux1,scalars(34))
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,ux2,scalars(35))
i=35

ux1=0  ! it's a deriavitve, so we know mean=0
do n=1,3
! zero crossing of dissipation: 
call transpose_to_x(grads(1,1,1,n),work,n1,n1d,n2,n2d,n3,n3d)
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,ux1,scalars(i+1))
call transpose_to_y(grads(1,1,1,n),work,n1,n1d,n2,n2d,n3,n3d)
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,ux1,scalars(i+2))
call transpose_to_z(grads(1,1,1,n),work,n1,n1d,n2,n2d,n3,n3d)
call compute_zero_crossing(work,n1,n1d,n2,n2d,n3,n3d,ux1,scalars(i+3))
i=i+3
enddo
! compute_zero_crossing does the MPI_allreduce for us, but we
! have to do it again below.  To keep the second one from having any
! effect, set all but one of the values to zero:
if (io_pe /= my_pe) then
   scalars(27:44)=0
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! i = 44
do n=1,3
   i=i+1
   scalars(i)=u2ux2(n)        
enddo
! i=47

scalars(48)=lnux1
scalars(49)=lnux2
scalars(50)=lnux3
scalars(51)=lnux4
i=51



if (i/=ns) then
   call abort("compute_expensive_pscalars: Error: i/=ns")
endif


#ifdef USE_MPI
   scalars2=scalars
   call MPI_allreduce(scalars2,scalars,i,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
!su=u2*uxx2/(ux2*ux2)   scalars(2)*scalars(12,13,14)/scalars(3,4,5)
!print *,'G_theta=',su


end subroutine







subroutine compute_zero_crossing(s,n1,n1d,n2,n2d,n3,n3d,avedata,N0)
use mpi
use params
implicit none
! input
integer :: n1,n1d,n2,n2d,n3,n3d
real*8 :: s(n1d,n2d,n3d)
! output
real*8 :: len,N0,xtmp
!local
integer i,j,k,n,count,raycount,ierr,count1
real*8 :: avedata,ave2,avelen

! compute zero crossing length in first direction 
! (all on processor:)
raycount=0
N0=0
avelen=0
do k=1,n3
do j=1,n2
   raycount=raycount+1  ! number of rays on this cpu
   call zero_crossing(s(1,j,k),n1,avedata,len,count1)
   N0=N0+count1   ! number of zero crossings along j,k ray

   ! len not computed anymore.  just average 'count1'
!   if (count1>0) then
!      count=count+1
!      avelen=avelen+len
!   endif
enddo
enddo

! number of rays: raycount
! number of zero crossings sumed over all rays:  N0
!    (use floating point, at 2048^3 could overflow interger*4)

#ifdef USE_MPI
xtmp=N0
call MPI_allreduce(xtmp,N0,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
!ave2=avelen
!call MPI_allreduce(ave2,avelen,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
count1=raycount
call MPI_allreduce(count1,raycount,1,MPI_INTEGER,MPI_SUM,comm_3d,ierr)
#endif

N0 = N0/raycount
!avelen=avelen/count

end subroutine





subroutine zero_crossing(data,n,ave,len,count)
integer :: n,count
real*8 :: data(n),len,ave
real*8 :: y1,y2
real*8 :: zero_location(n)

count=0
! find first crossing:
do i=1,n
   i1=i+1
   if (i1>n) i1=1
   y1=data(i)-ave
   y2=data(i1)-ave
   if ( (y1*y2<=0) .and. y2/=0 ) then
      ! find the location of the zero between [ data(i),data(i+1) )
      count=count+1
!      zero_location(count)=i +y1/(y1-y2)
   endif
enddo

return


if (count>0) then
   len=0
   do i=1,count
      if (i==count) then
         len=len+(1+zero_location(1))-zero_location(i)
      else
         len=len+zero_location(i+1)-zero_location(i)
      endif
   enddo
   len=len/count
endif
end subroutine




