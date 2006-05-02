#include "macros.h"



subroutine multfile_io(time,Q,iflag)
!
!  iflag=0     write .header file
!        1     write .u file  
!        2     write .v file  
!        3     write .w file  
!
use params
implicit none
real*8 :: Q(nx,ny,nz)
real*8 :: time
integer :: iflag

! local variables
integer i,j,k,n
real*8 xnx,xny,xnz,xnv
real*4 :: buf(nx)
character(len=80) message,ftime
character(len=240) fname
character(len=20) tmp
character(len=4) :: extension="uvwX"
CPOINTER :: fid
integer ierr

n=max(mpidims(1),mpidims(2),mpidims(3))
if (n<10) then
   n=5
else if (n<100) then
   n=4
else if (n<1000) then
   n=3
else if (n<10000) then
   n=2
else 
   call abort("opps, we assumed no more than 10000 cpus along one direction!")
endif
write(tmp,'(i5)') 10000+my_x
message="-" // tmp(n:5)
write(tmp,'(i5)') 10000+my_y
message=message(1:len_trim(message)) // "-" // tmp(n:5)
write(tmp,'(i5)') 10000+my_z
message=message(1:len_trim(message)) // "-" // tmp(n:5) // "-"

write(tmp,'(f10.4)') 10000.0000 + time
ftime=message(1:len_trim(message)) // tmp(2:10)


if (iflag==0) then
fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) //  &
        ftime(1:len_trim(ftime)) // ".header"
call copen(fname,"w",fid,ierr)
if (ierr/=0) then
   write(message,'(a,i5)') "multfile_io(): Error opening file errno=",ierr
   call abort(message)
endif

call cwrite8(fid,time,1)
xnv=n_var
xnx=nslabx
xny=nslaby
xnz=nslabz
call cwrite8(fid,xnx,1)
call cwrite8(fid,xny,1)
call cwrite8(fid,xnz,1)
call cwrite8(fid,xnv,1)
call cclose(fid,ierr)
endif

if (iflag>0) then
fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) &
   // ftime(1:len_trim(ftime)) // "." // extension(iflag:iflag)
if (output_size==4) then
   fname = fname(1:len_trim(fname)) // "4"
endif

call copen(fname,"w",fid,ierr)
if (ierr/=0) then
   write(message,'(a,i5)') "multfile_io(): Error opening file errno=",ierr
   call abort(message)
endif

if (output_size==8) then
   do k=nz1,nz2
   do j=ny1,ny2
      call cwrite8(fid,Q(nx1,j,k),nx2-nx1+1)
   enddo
   enddo
else
   do k=nz1,nz2
   do j=ny1,ny2
      buf=Q(:,j,k)
      call cwrite4(fid,buf(nx1),nx2-nx1+1)
   enddo
   enddo
endif

call cclose(fid,ierr)
endif







end subroutine













subroutine output_scalars(time,ints_save,maxs_save,nv,nscalars)
use params
implicit none
real*8 :: time
integer nv,nscalars
real*8 :: ints_save(nv,nscalars)
real*8 :: maxs_save(nv,nscalars)


! local variables
real*8 :: x
integer i,j,k,n,ierr
character(len=80) :: message
character,save :: access="0"
CPOINTER fid

! append to output files, unless this is first call.
if (access=="0") then
   access="w"
else
   access="a"
endif


if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time_initial
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".scalars"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      call print_message(message)
      write(message,'(a,i5)') "diag_output(): Error opening .scalars file errno=",ierr
      call abort(message)
   endif
   x=nv; call cwrite8(fid,x,1)
   x=nscalars; call cwrite8(fid,x,1)
   call cwrite8(fid,mu,1)
   call cwrite8(fid,alpha_value,1)
   call cwrite8(fid,ints_save,nv*nscalars);
   call cwrite8(fid,maxs_save,nv*nscalars);

   x=0; call cwrite8(fid,x,1)  ! for historical file format reasons
   call cwrite8(fid,time,1)

   call cclose(fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "diag_output(): Error closing .scalars file errno=",ierr
      call abort(message)
   endif
endif

end subroutine






subroutine singlefile_io(time,p,fname,work,work2,io_read,fpe)
!
! I/O routines where all data goes through a single PE and is
! written to a single file
!
! io_read=0    write data to file fname
! io_read=1    read data from file fname
!
! fpe       processor to do the file I/O
!
use params
use mpi
use transpose
implicit none
integer :: io_read  ! =1 for read, 0 for write
integer :: fpe
real*8 :: time
real*8 :: p(nx,ny,nz)
real*8 :: work2(nx,ny,nz),work(nx,ny,nz)
character(len=*) :: fname

call singlefile_io3(time,p,fname,work,work2,io_read,fpe,.false.,1)

end subroutine


subroutine singlefile_io2(time,p,fname,work,work2,io_read,fpe,output_spec)
!
! I/O routines where all data goes through a single PE and is
! written to a single file
!
! io_read=0    write data to file fname
! io_read=1    read data from file fname
!
! fpe       processor to do the file I/O
!
use params
use mpi
use transpose
implicit none
integer :: io_read  ! =1 for read, 0 for write
integer :: fpe
real*8 :: time
real*8 :: p(nx,ny,nz)
real*8 :: work2(nx,ny,nz),work(nx,ny,nz)
character(len=*) :: fname
logical :: output_spec

call singlefile_io3(time,p,fname,work,work2,io_read,fpe,output_spec,1)

end subroutine




subroutine singlefile_io3(time,p,fname,work,work2,io_read,fpe,&
output_spec,header_type)
!
! I/O routines where all data goes through a single PE and is
! written to a single file
!
! io_read=0    write data to file fname
! io_read=1    read data from file fname
!
! fpe       processor to do the file I/O
!
! output_spec=.true.     i/o on 2/3 dealiased spectral coefficients
! output_spec=.false.    i/o on full 3d array p 
!                     (grid point data or non-dealiased spec coeff.)
!
!
! header_type:    1   my default header, periodic extension
!                 2   no headers, no periodic extension      
!                 3   ensight headers, no periodic extension
!                 4   4 byte header  (fortran record?), no periodic extension
!                 5   same as 1, but also include Lz scale parameter
!   
! for header_type 2-4:  also disable output of the periodic extension data
!
use params
use mpi
use transpose
implicit none
integer :: io_read  ! =1 for read, 0 for write
integer :: fpe
real*8 :: time
real*8 :: p(nx,ny,nz)
real*8 :: work2(nx,ny,nz),work(nx,ny,nz)
character(len=*) :: fname
logical :: output_spec
integer im_max,km_max,jm_max

! local variables
integer :: i,j,k,n, header_type
real*8 :: xnx,xny,xnz
character(len=80) :: message
integer :: n_var_start,ierr
logical :: do_mpi_io_save
integer :: o_nx_save,o_ny_save,o_nz_save
CPOINTER :: fid


! disable periodic extension for header_type<>1 (turn it back on
! when we are done)
o_nx_save=o_nx; 
o_ny_save=o_ny; 
o_nz_save=o_nz;
if (header_type==2 .or. header_type==3 .or. header_type==4) then
   o_nx=g_nx; 
   o_ny=g_ny;
   o_nz=g_nz;
endif

if (do_mpi_io .and. .not. output_spec ) then
   call singlefile_mpi_io(time,p,fname,work,work2,io_read,header_type)
   o_nx=o_nx_save; o_ny=o_ny_save; o_nz=o_nz_save;
   return
endif


! disable MPI_IO for the calls below:
do_mpi_io_save=do_mpi_io
do_mpi_io=.false.

xnx=o_nx
xny=o_ny
xnz=o_nz


if (output_spec) then
   ! determine number of spectral coefficients to output:
   if (spec_max<0) then
      ! user did not specify:
      xnx=g_nx
      xny=g_ny
      xnz=g_nz
      if (dealias==1) then
         xnx=2+2*(g_nx/3)
         xny=2+2*(g_ny/3)
         xnz=2+2*(g_nz/3)
      endif
   else
      ! user specified
      xnx=min(spec_max,g_nx)
      xny=min(spec_max,g_ny)
      xnz=min(spec_max,g_nz)
   endif
endif


if (my_pe==fpe) then
   if (io_read==1) then
      call copen(fname,"r",fid,ierr)
   else
      call copen(fname,"w",fid,ierr)
   endif


   if (ierr/=0) then
      write(message,'(a,i5)') "singlefile_io(): Error opening file. Error no=",ierr
      call print_message(message)
      call print_message(fname)
      call abort("")
   endif
   if (header_type==1) call header1_io(io_read,output_spec,fid,time,xnx,xny,xnz)
   if (header_type==2) call header2_io(io_read,fid,time,xnx,xny,xnz)
   if (header_type==3) call header3_io(io_read,fid,time,xnx,xny,xnz)
   if (header_type==4) call header4_io(io_read,fid,time,xnx,xny,xnz)
   if (header_type==5) call header5_io(io_read,output_spec,fid,time,xnx,xny,xnz)
endif

#ifdef USE_MPI
if (io_read==1) then
   call mpi_bcast(time,1,MPI_REAL8,io_pe,comm_3d ,ierr)
   call mpi_bcast(xnx,1,MPI_REAL8,io_pe,comm_3d ,ierr)
   call mpi_bcast(xny,1,MPI_REAL8,io_pe,comm_3d ,ierr)
   call mpi_bcast(xnz,1,MPI_REAL8,io_pe,comm_3d ,ierr)
endif
#endif

! max wave number to read/write for spectral I/O
im_max=xnx/2-1
jm_max=xny/2-1
km_max=xnz/2-1

if (io_read==1) then
   if (output_spec) then
      call input1_spec(p,work,work2,fid,fpe,im_max,jm_max,km_max)
   else
      call input1(p,work,work2,fid,fpe,.false.,-1)
   endif
else
   if (output_spec) then
      call output1_spec(p,work,work2,fid,fpe,im_max,jm_max,km_max)
   else
      call output1(p,work,work2,fid,fpe,-1)
   endif

endif

if (my_pe==fpe) call cclose(fid,ierr)

o_nx=o_nx_save; o_ny=o_ny_save; o_nz=o_nz_save;
do_mpi_io=do_mpi_io_save


end subroutine







subroutine singlefile_mpi_io(time,p,fname,work,work2,io_read,header_type)
!
! I/O routines using mpi_io to write a single file
!
! io_read=0    write data to file fname
! io_read=1    read data from file fname
!
! fpe       processor to do the file I/O
!
!
!
use params
use mpi
use transpose
implicit none
integer :: io_read  ! =1 for read, 0 for write
integer :: fpe,header_type
real*8 :: time
real*8 :: p(nx,ny,nz)
real*8 :: work2(nx,ny,nz),work(nx,ny,nz)
character(len=*) :: fname
integer im_max,km_max,jm_max


! local variables
integer i,j,k,n
real*8 xnx,xny,xnz
character(len=80) message
integer :: n_var_start,ierr,offset=0
integer*8 infoin



CPOINTER fid

#ifndef USE_MPI_IO
call abort("singlefile_mpi_io: error: code not compiled with MPI-IO")   
#else



!
! do a check to see if all cpus can read /scratch2:
! if not, recompute io_nodes:
! call mpi_init_io(0)
!


fpe=io_nodes(0)

xnx=o_nx
xny=o_ny
xnz=o_nz


if (io_nodes(my_z)==my_pe) then
   call print_message("MPI-IO open...")
   if (io_read==1) then
      call mpi_file_open(comm_io,fname, &
           MPI_MODE_RDONLY ,&
           MPI_INFO_NULL, fid,ierr)
   else
      ! if file exists, then new stripe and stride settings are ignored,
      ! so delete file first
      call mpi_file_delete(fname,MPI_INFO_NULL,ierr)
      call mpi_info_create(infoin,ierr)
      call mpi_info_set(infoin, "striping_factor", mpi_stripe,ierr) 	
      call mpi_info_set(infoin, "striping_unit", mpi_stride ,ierr)
      call mpi_file_open(comm_io,fname, &
           MPI_MODE_WRONLY + MPI_MODE_CREATE ,&
           infoin, fid,ierr)
   endif
   if (ierr/=0) then
      write(message,'(a,i5)') "singlefile_mpi_io(): Error opening file. Error no=",ierr
      call print_message(message)
      call print_message(fname)
      call abort("")
   endif
   call print_message("MPI-IO open finished")
endif



if (my_pe==fpe) then
   if (header_type==1) then
      call header1_io(io_read,.false.,fid,time,xnx,xny,xnz)
   endif
   if (header_type==2) then
      call header2_io(io_read,fid,time,xnx,xny,xnz)  
   endif
   if (header_type==3) then
      call header3_io(io_read,fid,time,xnx,xny,xnz)
   endif
   if (header_type==5) then
      call header5_io(io_read,.false.,fid,time,xnx,xny,xnz)
   endif
endif

if (header_type==1) then
   offset=4 + o_nx + o_ny + o_nz
endif
if (header_type==2) then
   offset=0
endif
if (header_type==3) then
   if (io_read==1) then
      offset=240/output_size  ! offset multiplied by output_size later 
   else
      offset=240/input_size  ! offset multiplied by output_size later 
   endif
endif
if (header_type==4) then
   if (io_read==1) then
      offset=4/output_size  ! offset multiplied by output_size later 
   else
      offset=4/input_size  ! offset multiplied by output_size later 
   endif
   if (offset==0) then
      ! we should convert offset to total number of bytes!
      ! need to change this code and input1(), output1()
      call abort("Error: MPI-IO can not handle header_type==4 with real*8 data") 
   endif
endif
if (header_type==5) then
   offset=5 + o_nx + o_ny + o_nz
endif

#ifdef USE_MPI
call mpi_bcast(time,1,MPI_REAL8,fpe,comm_3d ,ierr)
#endif

if (io_read==1) then
   call input1(p,work,work2,fid,fpe,.false.,offset)
else
   call output1(p,work,work2,fid,fpe,offset)
endif

if (io_nodes(my_z)==my_pe) then
   call mpi_file_close(fid,ierr)
endif
#endif
end subroutine


























subroutine input_uvw(time,Q,Qhat,work1,work2,header_type)
!
! low wave number, quasi isotropic initial condition
!
use params
use mpi
use tracers
use fft_interface
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)

!local
character(len=80) message
character(len=280) fname
character(len=280) base
integer :: n,header_type,im,jm,km,i,j,k,ierr
real*8 :: time
real*8 :: time_in,mx(3),mx2,ke,ens,uy,uz,vx,vz,wx,wy,u2,xfac


Q=0


if (time<0) then
   base=rundir(1:len_trim(rundir)) // "restart"
else
   write(message,'(f10.4)') 10000.0000 + time
   base=rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) 
endif

if (r_compressed) then
   call uvw_compressed_io(base,time,1,Q,Qhat,work1,work2)
   return
endif


if (equations==NS_UVW) then

if (udm_input) then
   call udm_read_uvw(time_in,base,Q,work1,work2)
else
   if (r_spec) then
      if (header_type /= 1) call abort("Error: spectral I/O requires header_type==1");
      fname = base(1:len_trim(base)) // ".us"
      call print_message("Input: ")
      call print_message(trim(fname))
      call singlefile_io2(time_in,Q(1,1,1,1),fname,work1,work2,1,io_pe,r_spec)
      fname = base(1:len_trim(base)) // ".vs"
      call print_message(fname)
      call singlefile_io2(time_in,Q(1,1,1,2),fname,work1,work2,1,io_pe,r_spec)
      if (ndim==3) then
         fname = base(1:len_trim(base)) // ".ws"
         call print_message(trim(fname))
         call singlefile_io2(time_in,Q(1,1,1,ndim),fname,work1,work2,1,io_pe,r_spec)
      endif
      ! compute and print KE, ens
      ke=0; ens=0
      do k=nz1,nz2
         km=kmcord(k)
         do j=ny1,ny2
            jm=jmcord(j)
            do i=nx1,nx2
               im=imcord(i)
               ! u_x term
               vx = - im*pi2*Q(i+imsign(i),j,k,2)
               wx=0
               if (ndim==3) wx = - im*pi2*Q(i+imsign(i),j,k,ndim)
               uy = - jm*pi2*Q(i,j+jmsign(j),k,1)
               wy=0
               if (ndim==3) wy = - jm*pi2*Q(i,j+jmsign(j),k,ndim)
               uz =  - km*pi2*Q(i,j,k+kmsign(k),1)
               vz =  - km*pi2*Q(i,j,k+kmsign(k),2)
               
               u2=0
               do n=1,ndim
                  u2=u2+Q(i,j,k,n)*Q(i,j,k,n)
               enddo
            
               xfac = 2*2*2
               if (km==0) xfac=xfac/2
               if (jm==0) xfac=xfac/2
               if (im==0) xfac=xfac/2
               
               ke = ke + .5*xfac*u2
               ens = ens + xfac* &           
                    ((wy-vz)**2 + (uz-wx)**2 + (vx-uy)**2)   

            enddo
         enddo
      enddo
#ifdef USE_MPI
      mx2=ke
      call mpi_allreduce(mx2,ke,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
      mx2=ens
      call mpi_allreduce(mx2,ens,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
      do n=1,ndim
         call ifft3d(Q(1,1,1,n),work1)
         mx(n)=maxval(abs(Q(nx1:nx2,ny1:ny2,nz1:nz2,n)))
#ifdef USE_MPI
         mx2=mx(n)
         call mpi_allreduce(mx2,mx(n),1,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif
      enddo
      write(message,'(a,3f18.14)') 'spec_io: maxW = ',mx
      call print_message(message)
      write(message,'(a,3f18.14)') 'spec_io: ke =   ',ke
      call print_message(message)
      write(message,'(a,3f18.14)') 'spec_io: ens =  ',ens
      call print_message(message)


   else
      fname = base(1:len_trim(base)) // ".u"
      call print_message("Input: ")
      call print_message(trim(fname))
      call singlefile_io3(time_in,Q(1,1,1,1),fname,work1,work2,1,io_pe,r_spec,header_type)
      fname = base(1:len_trim(base)) // ".v"
      call print_message(trim(fname))
      call singlefile_io3(time_in,Q(1,1,1,2),fname,work1,work2,1,io_pe,r_spec,header_type)
      if (ndim==3) then
         fname = base(1:len_trim(base)) // ".w"
         call print_message(trim(fname))
         call singlefile_io3(time_in,Q(1,1,1,ndim),fname,work1,work2,1,io_pe,r_spec,header_type)
      endif
   endif
endif

else if (equations==SHALLOW) then
   fname = base(1:len_trim(base)) // ".u"
   if (r_spec) fname=fname(1:len_trim(fname)) // "s"
   call print_message("Input: ")
   call print_message(trim(fname))
   call singlefile_io3(time_in,Q(1,1,1,1),fname,work1,work2,1,io_pe,r_spec,header_type)

   fname = base(1:len_trim(base)) // ".v"
   if (r_spec) fname=fname(1:len_trim(fname)) // "s"
   call print_message(trim(fname))
   call singlefile_io3(time_in,Q(1,1,1,2),fname,work1,work2,1,io_pe,r_spec,header_type)
   fname = base(1:len_trim(base)) // ".h"
   if (r_spec) fname=fname(1:len_trim(fname)) // "s"
   call print_message(trim(fname))
   call singlefile_io3(time_in,Q(1,1,1,n_var),fname,work1,work2,1,io_pe,r_spec,header_type)

   if (r_spec) then
      ! transform to grid space
      do n=1,n_var
         call ifft3d(Q(1,1,1,n),work1)
      enddo
   endif


else if (equations==NS_PSIVOR) then
   if (header_type /= 1) call abort("Error: NS_PSIVOR I/O requires header_type==1");
   call tracers_restart(io_pe)
   fname = base(1:len_trim(base)) // ".vor"
   call print_message("Input: ")
   call print_message(trim(fname))
   call singlefile_io(time_in,Qhat(1,1,1,1),fname,work1,work2,1,io_pe)
   !fname = base(1:len_trim(base)) // ".psi"
   !call singlefile_io(time_in,Qhat,fname,work1,work2,1,io_pe)

else if (equations==CNS) then
   fname = base(1:len_trim(base)) // ".u"
   call print_message("Input: ")
   call print_message(trim(fname))
   call singlefile_io3(time_in,Q(1,1,1,1),fname,work1,work2,1,io_pe,r_spec,header_type)
   fname = base(1:len_trim(base)) // ".v"
   call print_message(trim(fname))
   call singlefile_io3(time_in,Q(1,1,1,2),fname,work1,work2,1,io_pe,r_spec,header_type)
   if (ndim==3) then
      fname = base(1:len_trim(base)) // ".w"
      call print_message(trim(fname))
      call singlefile_io3(time_in,Q(1,1,1,ndim),fname,work1,work2,1,io_pe,r_spec,header_type)
   endif
   if (n_var-ndim==2) then
      fname = base(1:len_trim(base)) // ".p"
      n=ndim+1
      call singlefile_io3(time_in,Q(1,1,1,n),fname,work1,work2,1,io_pe,r_spec,header_type)
      fname = base(1:len_trim(base)) // ".rho"
      n=ndim+2
      call singlefile_io3(time_in,Q(1,1,1,n),fname,work1,work2,1,io_pe,r_spec,header_type)
   endif  
endif

!these header types can read the time from the input file:
! for other headers, dont change 'time'
if (header_type==1 .or. header_type==5) then
   time=time_in
endif
end subroutine




subroutine output_uvw(basename,time,Q,q1,work1,work2,header_type)
use params
use tracers
use mpi
use fft_interface
use transpose
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
character(len=*) :: basename
integer :: header_type

!local
character(len=80) message
character(len=280) fname
character(len=80) base
integer :: n
real*8 :: time

write(message,'(f10.4)') 10000.0000 + time

if (w_compressed) then
   fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) 
   call uvw_compressed_io(fname,time,0,Q,q1,work1,work2)
   return
endif




if (equations==NS_UVW .and. w_spec) then
   
   ! NS, primitive variables
   fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".us"
   call singlefile_io3(time,Q(1,1,1,1),fname,work1,work2,0,io_pe,.true.,header_type)
   
   fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".vs"
   call singlefile_io3(time,Q(1,1,1,2),fname,work1,work2,0,io_pe,.true.,header_type)
   if (ndim==3) then
      fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".ws"
      call singlefile_io3(time,Q(1,1,1,ndim),fname,work1,work2,0,io_pe,.true.,header_type)
   endif
   
else if (equations==NS_UVW) then
   if (udm_output) then
      fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".h5"
      call udm_write_uvw(fname,time,Q,q1)
   else
      ! NS, primitive variables
      fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".u"
      call singlefile_io3(time,Q(1,1,1,1),fname,work1,work2,0,io_pe,.false.,header_type)
      fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".v"
      call singlefile_io3(time,Q(1,1,1,2),fname,work1,work2,0,io_pe,.false.,header_type)
      if (ndim==3) then
         fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".w"
         call singlefile_io3(time,Q(1,1,1,ndim),fname,work1,work2,0,io_pe,.false.,header_type)
      endif
      if (ndim==2) then
         call vorticity2d(q1,Q,work1,work2)
         fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".vor"
         call singlefile_io3(time,q1,fname,work1,work2,0,io_pe,.false.,header_type)
      endif
   endif
   
else if (equations==SHALLOW) then
   ! shallow water 2D
   fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".u"
   if (w_spec) fname=fname(1:len_trim(fname)) // "s"
   call singlefile_io3(time,Q(1,1,1,1),fname,work1,work2,0,io_pe,w_spec,header_type)

   fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".v"
   if (w_spec) fname=fname(1:len_trim(fname)) // "s"
   call singlefile_io3(time,Q(1,1,1,2),fname,work1,work2,0,io_pe,w_spec,header_type)

   fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".h"
   if (w_spec) fname=fname(1:len_trim(fname)) // "s"
   call singlefile_io3(time,Q(1,1,1,n_var),fname,work1,work2,0,io_pe,w_spec,header_type)

   if (ndim==2 .and. .not. w_spec) then
      call vorticity2d(q1,Q,work1,work2)
      fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".vor"
      call singlefile_io(time,q1,fname,work1,work2,0,io_pe)
   endif

else if (equations==NS_PSIVOR) then
   call tracers_save(io_pe,time)
   ! 2D NS psi-vor formulation
   fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".vor"
   call singlefile_io(time,Q(1,1,1,1),fname,work1,work2,0,io_pe)
   fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".psi"
   call singlefile_io(time,Q(1,1,1,2),fname,work1,work2,0,io_pe)
else if (equations==CNS) then
   ! NS, primitive variables
   fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".u"
   call singlefile_io3(time,Q(1,1,1,1),fname,work1,work2,0,io_pe,.false.,header_type)
   fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".v"
   call singlefile_io3(time,Q(1,1,1,2),fname,work1,work2,0,io_pe,.false.,header_type)
   if (ndim==3) then
      fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".w"
      call singlefile_io3(time,Q(1,1,1,ndim),fname,work1,work2,0,io_pe,.false.,header_type)
   endif
   if (ndim==2) then
      call vorticity2d(q1,Q,work1,work2)
      fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".vor"
      call singlefile_io3(time,q1,fname,work1,work2,0,io_pe,.false.,header_type)
   endif
   if (n_var-ndim==2) then
      fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".p"
      n=ndim+1
      call singlefile_io3(time,Q(1,1,1,n),fname,work1,work2,0,io_pe,.false.,header_type)
      fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".rho"
      n=ndim+2
      call singlefile_io3(time,Q(1,1,1,n),fname,work1,work2,0,io_pe,.false.,header_type)
      
   endif
endif

if (ndim==3 .and. output_vorticity/=0) then
   call print_message("outputting vorticity...")
   call vorticity(q1,Q,work1,work2)
   fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".vor1"
   call singlefile_io3(time,q1(1,1,1,1),fname,work1,work2,0,io_pe,.false.,header_type)
   fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".vor2"
   call singlefile_io3(time,q1(1,1,1,2),fname,work1,work2,0,io_pe,.false.,header_type)
   if (ndim==3) then
      fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) // message(2:10) // ".vor3"
      call singlefile_io3(time,q1(1,1,1,ndim),fname,work1,work2,0,io_pe,.false.,header_type)
   endif
endif


end subroutine






subroutine uvw_compressed_io(basename,time,io_read,Q,Qhat,work1,work2)
!
!  compressed I/O routine
!  io_read=0   write data
!  io_read=1   read data
!
use params
use tracers
use mpi
use fft_interface
use transpose
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(g_nz2,nslabx,ny_2dz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
integer :: io_read
character(len=*) :: basename

!local
character(len=80) message
character(len=280) fname
character(len=80) base
integer :: n,im,jm,km,ierr,i,j,k
integer n1,n1d,n2,n2d,n3,n3d
logical :: do_mpi_io_save
real*8 :: time,mx,mx2,ke,vx,wx,uy,wy,uz,vz,ens,xfac,u2
CPOINTER :: fid


write(message,'(f10.4)') 10000.0000 + time
if (equations/=NS_UVW) then
   call abort("compressed I/O requires NS_UVW equations")
endif
if (ndim/=3) then
   call abort("compressed I/O requires ndim==3")
endif
if (ncpu_x*ncpu_y>1) then
   call abort("compressed I/O requies ncpu_x=ncpu_y=1")
endif
if (nx2-nx1+1 /= g_nx) then
   call abort("compressed I/O x dimension error")
endif
if (ny2-ny1+1 /= g_ny) then
   call abort("compressed I/O y dimension error")
endif

if (io_read) then
   call print_message("reading compressed I/O:")
else
   call print_message("writing compressed I/O:")
endif

fname = basename(1:len_trim(basename)) // ".uc"
call print_message(trim(fname))
call singlefile_io3(time,Q(1,1,1,1),fname,work1,work2,io_read,io_pe,.false.,2)
if (io_read) then
   mx=maxval(abs(Q(nx1:nx2,ny1:ny2,nz1:nz2,1)))
#ifdef USE_MPI
   mx2=mx
   call mpi_allreduce(mx2,mx,1,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif
write(message,'(a,3f18.14)') 'compressed_io: maxU = ',mx
call print_message(message)
endif


   
fname = basename(1:len_trim(basename)) // ".vc"
call print_message(trim(fname))
call singlefile_io3(time,Q(1,1,1,2),fname,work1,work2,io_read,io_pe,.false.,2)
if (io_read) then
   mx=maxval(abs(Q(nx1:nx2,ny1:ny2,nz1:nz2,2)))
#ifdef USE_MPI
   mx2=mx
   call mpi_allreduce(mx2,mx,1,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif
write(message,'(a,3f18.14)') 'compressed_io: maxV = ',mx
call print_message(message)
endif



fname = basename(1:len_trim(basename)) // ".wc"
call print_message(trim(fname))

! temp. turn off MPI-IO, if it was enabled 
do_mpi_io_save=do_mpi_io
do_mpi_io=.false.


if (io_read==0) then
   ! output u,v and z-averaged w
   ! this subroutine requires ncpux=ncpuy=1, so now have everyone compute
   ! their z-ave and then reduce to io_pe
   work1=0
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      work1(i,j,1)=work1(i,j,1)+Q(i,j,k,ndim)/g_nz
   enddo
   enddo
   enddo
#ifdef USE_MPI
   call mpi_reduce(work1,work2,nx*ny,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
#else
   work2=work1
#endif


   if (my_pe==io_pe) then
      call copen(fname,"w",fid,ierr)
      do j=ny1,ny2
         call mwrite8(fid,work2(nx1,j,1),nx2-nx1+1)
      enddo
      call cclose(fid,ierr)
   endif
else
   ! read in u,v and z-averaged w
   ! then recover w from div(u)=0
   work2=0
   if (my_pe==io_pe) then
      call copen(fname,"r",fid,ierr)
      do j=ny1,ny2
         call mread8(fid,work2(nx1,j,1),nx2-nx1+1)
      enddo
      call cclose(fid,ierr)
   endif
#ifdef USE_MPI
   call mpi_bcast(work2,nx*ny,MPI_REAL8,io_pe,comm_3d,ierr)
#endif
   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      Q(i,j,k,ndim)=work2(i,j,1)
   enddo
   enddo
   enddo

   ! now use the divergence relation to recover the rest of Q(:,:,:,3)
   ! in wave number space:  (l,m,n)   
   !   i l*u + m*v + n*w = 0    so for n<>0    w = -(l*u + m*v)/n
   ! and the n=0 component has already been computed and stored in
   ! Q(:,:,:,3) above.  
   call print_message("recovering u_3 from div(u)=0 relation...")
   do n=1,3
      work1=Q(:,:,:,n)
      call z_fft3d_trashinput(work1,Qhat(1,1,1,n),work2)
   enddo

   do j=1,ny_2dz
      jm=z_jmcord(j)
      do i=1,nslabx
         im=z_imcord(i)
         do k=1,g_nz
            km=z_kmcord(k)
            
            ! the divergence
            ! 0 =  - im*rhs(k,i+z_imsign(i),j,1) &
            !     - jm*rhs(k,i,j+z_jmsign(j),2) &
            !     - km*rhs(k+z_kmsign(k),i,j,3)/Lz

            if (km/=0) then
               Qhat(k+z_kmsign(k),i,j,3) = Qhat(k+z_kmsign(k),i,j,3) + &
(-Lz/km) * ( im*Qhat(k,i+z_imsign(i),j,1) + jm*Qhat(k,i,j+z_jmsign(j),2) )
            endif

         enddo
      enddo
   enddo


   ke=0
   ens=0
   do j=1,ny_2dz
      jm=z_jmcord(j)
      do i=1,nslabx
         im=z_imcord(i)
         do k=1,g_nz
            km=z_kmcord(k)
            xfac = 2*2*2
            if (km==0) xfac=xfac/2
            if (jm==0) xfac=xfac/2
            if (im==0) xfac=xfac/2
            
            u2=Qhat(k,i,j,1)*Qhat(k,i,j,1) + &
                 Qhat(k,i,j,2)*Qhat(k,i,j,2) + &
                 Qhat(k,i,j,3)*Qhat(k,i,j,3)
            
            ke = ke + .5*xfac*u2

            ! u_x term
            vx = - pi2*im*Qhat(k,i+z_imsign(i),j,2)
            wx = - pi2*im*Qhat(k,i+z_imsign(i),j,3)
            uy = - pi2*jm*Qhat(k,i,j+z_jmsign(j),1)
            wy = - pi2*jm*Qhat(k,i,j+z_jmsign(j),3)
            uz =  - pi2*km*Qhat(k+z_kmsign(k),i,j,1)/Lz
            vz =  - pi2*km*Qhat(k+z_kmsign(k),i,j,2)/Lz
            ! vorcity: ( (wy - vz), (uz - wx), (vx - uy) )
            ens = ens + xfac* &           
                 ((wy-vz)**2 + (uz-wx)**2 + (vx-uy)**2)   

         enddo
      enddo
   enddo



#ifdef USE_MPI
   mx2=ke
   call mpi_allreduce(mx2,ke,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   mx2=ens
   call mpi_allreduce(mx2,ens,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
   do n=1,3
      call z_ifft3d(Qhat(1,1,1,n),Q(1,1,1,n),work1)
   enddo

   mx=maxval(abs(Q(nx1:nx2,ny1:ny2,nz1:nz2,3)))
#ifdef USE_MPI
   mx2=mx
   call mpi_allreduce(mx2,mx,1,MPI_REAL8,MPI_MAX,comm_3d,ierr)
#endif
   write(message,'(a,3f18.14)') 'compressed_io: maxW = ',mx
   call print_message(message)
   write(message,'(a,3f18.14)') 'compressed_io: ke =   ',ke
   call print_message(message)
   write(message,'(a,3f18.14)') 'compressed_io: ens =  ',ens
   call print_message(message)

   call print_message("done with uncompress")
endif


! restore MPI-IO flag
do_mpi_io=do_mpi_io_save
   

end subroutine







subroutine output_passive(basename,time,Q,q1,work1,work2)
use params
use tracers
use mpi
use fft_interface
use transpose
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
character(len=*) :: basename

!local
character(len=80) message
character(len=80) fname
character(len=80) base,ext,ext2
integer :: n
real*8 :: time,dummy

if (npassive==0) return
write(message,'(f10.4)') 10000.0000 + time

do n=np1,np2
   write(ext,'(f8.3)') 1000 + schmidt(n) ! 000.000
   write(ext2,'(i3)') 100+passive_type(n)
   fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) &
        // message(2:10) // '.t' // ext2(2:3) // '.s' // ext(2:len_trim(ext))
   call singlefile_io(time,Q(1,1,1,n),fname,work1,work2,0,io_pe)
enddo

! Ray/tmix passive scalars: output plane dissipation
! output plane dissipaiton:  (mu/schmidt) * c_x**2 + c_y**2
! for passive scalars of type 0 or 1 (Gaussian or KE correlated)
if (passive_type(np1)<=1) then
   do n=np1,np2
      call print_message("outputting plane dissipation")
      call der(Q(1,1,1,n),work1,q1(1,1,1,1),q1(1,1,1,2),DX_ONLY,1)  ! c_x
      call der(Q(1,1,1,n),work2,q1(1,1,1,1),q1(1,1,1,2),DX_ONLY,2)  ! c_y
      q1(:,:,:,1)=mu*(work1**2 + work2**2)/schmidt(n)

      write(ext,'(f8.3)') 1000 + schmidt(n) ! 000.000
      write(ext2,'(i3)') 100+passive_type(n)
      fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) &
           // message(2:10) // '.t' // ext2(2:3) // '.s' // &
          ext(2:len_trim(ext))  // ".gradxy2"
      call singlefile_io(time,q1,fname,work1,work2,0,io_pe)
   enddo
endif


! SK passive scalars: output u+grad(s)
if (passive_type(np1)==2 .and. npassive==1) then
   call print_message("outputting u+grad(s)")
   do n=np1,np2
      call der(Q(1,1,1,n),q1(1,1,1,1),work1,work2,DX_ONLY,1)  ! c_x
      call der(Q(1,1,1,n),q1(1,1,1,2),work1,work2,DX_ONLY,2)  ! c_y
      call der(Q(1,1,1,n),q1(1,1,1,ndim),work1,work2,DX_ONLY,3)  ! c_z

      q1(:,:,:,1:ndim)=Q(:,:,:,1:ndim)+ q1(:,:,:,1:ndim)

      write(ext,'(f8.3)') 1000 + schmidt(n) ! 000.000
      write(ext2,'(i3)') 100+passive_type(n)
      fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) &
           // message(2:10) // '.gamma1'
      call singlefile_io(time,q1(1,1,1,1),fname,work1,work2,0,io_pe)
      fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) &
           // message(2:10) // '.gamma2'
      call singlefile_io(time,q1(1,1,1,2),fname,work1,work2,0,io_pe)
      if (ndim>=3) then
      fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) &
           // message(2:10) // '.gamma3'
      call singlefile_io(time,q1(1,1,1,ndim),fname,work1,work2,0,io_pe)
      endif
   enddo
endif

if (passive_type(np1)==4 .and. npassive==1) then
   ! compute PV in q1(:,:,:,1)
   call vorticity(q1,Q,work1,work2)  ! q1 = vorticity

   ! now compute    ( curl(u) + f-khat )  dot grad(theta) 
   ! q1(:,:,:,1) = vor1*theta_x + vor2*theta_y + (vor3+f)*theta_z

   ! overwrite q1(:,:,:,1) with x-component of vorticity times d(theta)/dx
   call der(Q(1,1,1,np1),work1,dummy,work2,DX_ONLY,1)  ! work1 = theta_x
   q1(:,:,:,1)=q1(:,:,:,1)*work1                       

   call der(Q(1,1,1,np1),work1,dummy,work2,DX_ONLY,2)  ! work1 = theta_y
   q1(:,:,:,1)=q1(:,:,:,1) + q1(:,:,:,2)*work1   

   call der(Q(1,1,1,np1),work1,dummy,work2,DX_ONLY,3)  ! work1 = theta_z
   q1(:,:,:,1)=q1(:,:,:,1) + (q1(:,:,:,3)+fcor)*work1   

   ! output into a file suffixed with ".pv"   
   fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) &
        // message(2:10) // '.pv' 
   call singlefile_io(time,q1,fname,work1,work2,0,io_pe)
endif

end subroutine







subroutine input_passive(basename,time,Q,work1,work2)
use params
use tracers
use mpi
use fft_interface
use transpose
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
character(len=*) :: basename

!local
character(len=80) message
character(len=80) fname
character(len=80) base,ext,ext2
integer :: n
real*8 :: time,mn,mx

if (npassive==0) return

!     <-time--> <t>  <-sch->
! name0000.0000.t00.s000.000
!
do n=np1,np2
   write(message,'(f10.4)') 10000.0000 + time
   write(ext,'(f8.3)') 1000 + schmidt(n) ! 000.000
   write(ext2,'(i3)') 100+passive_type(n)
   fname = rundir(1:len_trim(rundir)) // basename(1:len_trim(basename)) &
        // message(2:10) // '.t' // ext2(2:3) // '.s' // ext(2:8)
   call print_message(trim(fname))	
   call singlefile_io(time,Q(1,1,1,n),fname,work1,work2,1,io_pe)


   call global_min(Q(1,1,1,n),mn)
   call global_max(Q(1,1,1,n),mx)

   write(message,'(a,i3,a,2f17.5)') 'passive scalar n=',n,' min/max: ',mn,mx
   call print_message(message)	
enddo

end subroutine














subroutine header1_io(io_read,output_spec,fid,time,xnx,xny,xnz)
use params
use mpi
use transpose
implicit none
CPOINTER fid
integer :: io_read  ! =1 for read, 0 for write
real*8 :: time,xnx,xny,xnz
logical :: output_spec

! local
integer :: ierr,	i
character(len=80) message
real*8 :: temp(g_nz+1)


if (io_read==1) then
   call mread8e(fid,time,1,ierr)
   if (ierr/=1) then
      write(message,'(a,i5)') "header1_io(): Error reading file"
      call print_message(message)
      call abort("")
   endif
   call mread8(fid,xnx,1)
   call mread8(fid,xny,1)
   call mread8(fid,xnz,1)
   if (output_spec) then
      print *,'Spectrum input data'
      write(*,'(a,3f5.0)') 'number of real coefficients: ',xnx,xny,xnz
      if (xnx>g_nx .or. xny>g_ny .or. xnz>g_nz) then
         ! we can only upsample low-res data.
         ! to run with high-res data, output a trucated form.  
         call print_message("Error: spectral input requires downsampling to lower resolution")
         call print_message("Input routines can only upsample.") 
         call print_message("Output routines can only downsample.") 
         call print_message("Run code at higher resolution, calling Output to downsample")
         call abort("error in header1_io")
      endif
   else
      print *,'grid input data'
      write(*,'(a,3f7.0)') 'number of grid points: ',xnx,xny,xnz
      if (int(xnx)/=o_nx) call abort("Error: data file nx <> nx set in params.h");
      if (int(xny)/=o_ny) call abort("Error: data file ny <> ny set in params.h");
      if (int(xnz)/=o_nz) call abort("Error: data file nz <> nz set in params.h");
      call mread8(fid,g_xcord(1),o_nx)
      call mread8(fid,g_ycord(1),o_ny)
      call mread8(fid,g_zcord(1),o_nz)
   endif
else
   call mwrite8(fid,time,1)
   call mwrite8(fid,xnx,1)
   call mwrite8(fid,xny,1)
   call mwrite8(fid,xnz,1)
   if ( output_spec) then
      print *,'Spectrum output data'
      write(*,'(a,3f5.0)') 'number of real coefficients: ',xnx,xny,xnz
      if (xnx>g_nx .or. xny>g_ny .or. xnz>g_nz) then
         ! we can only downsample to lower res data on output.
         call print_message("Error: spectral output requires zero padding") 
         call print_message("Output routines can only downsample.") 
         call print_message("Input routines can input this data directly")
         call abort("error in header1_io")
      endif
   else
      call mwrite8(fid,g_xcord(1),o_nx)
      call mwrite8(fid,g_ycord(1),o_ny)
      call mwrite8(fid,temp(1),o_nz)
   endif
   
endif
end subroutine header1_io



subroutine header5_io(io_read,output_spec,fid,time,xnx,xny,xnz)
use params
use mpi
use transpose
implicit none
CPOINTER fid
integer :: io_read  ! =1 for read, 0 for write
real*8 :: time,xnx,xny,xnz
logical :: output_spec

! local
integer :: ierr,	i
character(len=80) message
real*8 :: temp(g_nz+1)


if (io_read==1) then
   call mread8e(fid,time,1,ierr)
   if (ierr/=1) then
      write(message,'(a,i5)') "header5_io(): Error reading file"
      call print_message(message)
      call abort("")
   endif
   call mread8(fid,Lz,1)
   ! run sanity check on Lz
   if (Lz < .001  .or. Lz > 20 ) then
      print *,'Error fileio.F90: input file aspect ratio out of range?'
      print *,'Lz = ',Lz
      call abort("abort...")
   endif

   call mread8(fid,xnx,1)
   call mread8(fid,xny,1)
   call mread8(fid,xnz,1)
   if (output_spec) then
      print *,'Spectrum input data'
      write(*,'(a,3f5.0)') 'number of real coefficients: ',xnx,xny,xnz
      if (xnx>g_nx .or. xny>g_ny .or. xnz>g_nz) then
         ! we can only upsample low-res data.
         ! to run with high-res data, output a trucated form.  
         call print_message("Error: spectral input requires downsampling to lower resolution")
         call print_message("Input routines can only upsample.") 
         call print_message("Output routines can only downsample.") 
         call print_message("Run code at higher resolution, calling Output to downsample")
         call abort("error in header5_io")
      endif
   else
      print *,'grid input data'
      write(*,'(a,3f7.0)') 'number of grid points: ',xnx,xny,xnz
      if (int(xnx)/=o_nx) call abort("Error: data file nx <> nx set in params.h");
      if (int(xny)/=o_ny) call abort("Error: data file ny <> ny set in params.h");
      if (int(xnz)/=o_nz) call abort("Error: data file nz <> nz set in params.h");
      call mread8(fid,g_xcord(1),o_nx)
      call mread8(fid,g_ycord(1),o_ny)
      call mread8(fid,g_zcord(1),o_nz)
   endif
else
   call mwrite8(fid,time,1)
   call mwrite8(fid,Lz,1)
   call mwrite8(fid,xnx,1)
   call mwrite8(fid,xny,1)
   call mwrite8(fid,xnz,1)
   if ( output_spec) then
      print *,'Spectrum output data'
      write(*,'(a,3f5.0)') 'number of real coefficients: ',xnx,xny,xnz
      if (xnx>g_nx .or. xny>g_ny .or. xnz>g_nz) then
         ! we can only downsample to lower res data on output.
         call print_message("Error: spectral output requires zero padding") 
         call print_message("Output routines can only downsample.") 
         call print_message("Input routines can input this data directly")
         call abort("error in header5_io")
      endif
   else
      call mwrite8(fid,g_xcord(1),o_nx)
      call mwrite8(fid,g_ycord(1),o_ny)
      call mwrite8(fid,temp(1),o_nz)
   endif
   
endif
end subroutine header5_io





subroutine header3_io(io_read,fid,time,xnx,xny,xnz)
!
! header is 3 lines of text, with each line exactly 80 characters
!
use params
use mpi
use transpose
implicit none
CPOINTER fid
integer :: io_read  ! =1 for read, 0 for write
real*8 :: time,xnx,xny,xnz
logical :: output_spec

! local
integer :: ierr,i
character(len=80) :: message

xnx=-1
xny=-1
xnz=-1

if (io_read==1) then
   call mread1e(fid,message,80,ierr)
   if (ierr/=1) then
      write(message,'(a,i5)') "header3_io(): Error reading file"
      call print_message(message)
      call abort("")
   endif
   call mread1e(fid,message,80,ierr)
   call mread1e(fid,message,80,ierr)
else
   message="turbulence" // char(0)
   call mwrite1(fid,message,80)
   message="part" // char(0)
   call mwrite1(fid,message,80)
   message="block" // char(0)
   call mwrite1(fid,message,80)
endif
end subroutine header3_io



subroutine header4_io(io_read,fid,time,xnx,xny,xnz)
!
! header is 4 bytes, fortran record counter
!
use params
use mpi
use transpose
implicit none
CPOINTER fid
integer :: io_read  ! =1 for read, 0 for write
real*8 :: time,xnx,xny,xnz
logical :: output_spec

! local
integer :: ierr,i
character(len=80) :: message

xnx=-1
xny=-1
xnz=-1

if (io_read==1) then
   ! read 4 bytes
   call mread1e(fid,message,4,ierr)
   if (ierr/=4) then
      write(message,'(a,i5)') "header4_io(): Error reading file"
      call print_message(message)
      call abort("")
   endif
else
   ! dont know how to write a fortran header
   call abort("ERROR: header_type==4 output not supported")
endif
end subroutine header4_io







subroutine header2_io(io_read,fid,time,xnx,xny,xnz)
use params
use mpi
use transpose
implicit none
CPOINTER fid
integer :: io_read  ! =1 for read, 0 for write
real*8 :: time,xnx,xny,xnz
logical :: output_spec

! local
integer :: ierr
character(len=80) message

if (io_read==1) time=-1  ! we cannot read time, so set to bad value
                         ! calling program must obtain time some other way
xnx=-1
xny=-1
xnz=-1

end subroutine header2_io





