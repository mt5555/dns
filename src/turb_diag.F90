#include "macros.h"
subroutine output_model(time,Q,Qhat,q1,q2,q3,work1,work2)
use params
use structf
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(*)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: q2(nx,ny,nz,n_var)
real*8 :: q3(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time

! local variables
integer,parameter :: nints_e=14
real*8 :: ints_e(nints_e)
real*8 :: x
integer i,j,k,n,ierr
character(len=80) :: message
CPOINTER fid,fidj,fidS





!
! output structure functions and time averaged forcing
! 
if (compute_struct==1) then
   call compute_all_pdfs(Q,q1,q2,q3,work1,work2,ints_e,nints_e)
   
   
   write(message,'(a,3f14.8)') 'skewness ux,vw,wz: ',&
        (ints_e(n+3)/ints_e(n)**1.5,n=1,3)
   call print_message(message)
   
   write(message,'(a,3f14.8)') 'wSw: ',&
        (ints_e(10)/ints_e(1)**2)
   call print_message(message)
   
   
   if (structf_init==1) then
   if (my_pe==io_pe) then
      write(message,'(f10.4)') 10000.0000 + time
      message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".sf"
      call copen(message,"w",fid,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "output_model(): Error opening .sf file errno=",ierr
         call abort(message)
      endif

      write(message,'(f10.4)') 10000.0000 + time
      message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".jpdf"
      call copen(message,"w",fidj,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "output_model(): Error opening .jpdf file errno=",ierr
         call abort(message)
      endif

      write(message,'(f10.4)') 10000.0000 + time
      message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".s2v2"
      call copen(message,"w",fidS,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "output_model(): Error opening .s2v2 file errno=",ierr
         call abort(message)
      endif

   endif
   call output_pdf(time,fid,fidj,fidS)
   if (my_pe==io_pe) call cclose(fid,ierr)
   if (my_pe==io_pe) call cclose(fidj,ierr)
   if (my_pe==io_pe) call cclose(fidS,ierr)
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
      if (ierr/=0) then
         write(message,'(a,i5)') "diag_output(): Error closing .scalars-turb file errno=",ierr
         call abort(message)
      endif
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

