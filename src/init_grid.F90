#include "macros.h"
subroutine init_model
use params
use fft_interface
use transpose
implicit none
character(len=80) ::  message
integer :: k



call params_init()
call transpose_init()
call fft_interface_init()
call init_input_file()
call init_grid()          


!
! scale alpha now that we know delx
!
if (alpha_value>=1e10) then
   infinite_alpha=1
else if (alpha_value>=1.0) then
   alpha_value=alpha_value*min(delx,dely,delz)
endif
if (infinite_alpha==0) then
   write(message,'(a,f14.8,f10.4)') "NS-Alpha:  alpha, alpha/h :",&
        alpha_value,alpha_value/min(delx,dely,delz)
else
   write(message,'(a,f14.8,f10.4)') "NS-Alpha:  infinite"
endif
call print_message(message)

write(message,'(a,i6,a,i6,a,i6)') "Global grid: ",g_nx," x",g_ny," x",g_nz
call print_message(message)
write(message,'(a,i6,a,i6,a,i6)') "Local grid (with padding): ",nx," x",ny," x",nz
call print_message(message)




end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! initialization for fft grid
! domain is from 0..1
! first point is at grid point 0, last point is at gridpoint
! 1-delta_x
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_grid
use params
use fft_interface
use mpi
implicit none

!local variables
real*8 :: one=1
integer i,j,k,l,ierr
character(len=80) message


if (numerical_method==FOURIER) then
!
! this is needed to gaurentee that the sine and cosine modes
! are on the same processor.  
!

if (mod(nslabx,2)/=0) then
   call abort("nslabx must be even")
endif
if (mod(nslaby,2)/=0) then
   call abort("nslaby must be even")
endif
if (mod(nslabz,2)/=0 .and. g_nz>1) then
   call abort("nslabz must be even if g_nz>1")
endif
endif


! periodic FFT case:  for output, we include the point at x=1 (same as x=0)
o_nx=g_nx
if (g_bdy_x1==PERIODIC) o_nx=g_nx+1

o_ny=g_ny
if (g_bdy_y1==PERIODIC) o_ny=g_ny+1

o_nz=g_nz
if (g_bdy_z1==PERIODIC) o_nz=g_nz+1
if (g_nz==1) o_nz=1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local grid, bounds over interior (non-boundary) points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
intx1=nx1
intx2=nx2
inty1=ny1
inty2=ny2
intz1=nz1
intz2=nz2

bx1=nx1
bx2=nx2
by1=ny1
by2=ny2
bz1=nz1
bz2=nz2

if (offset_bdy==1) then
   if (my_x==ncpu_x-1) bx2=nx2+1; o_nx=g_nx+1
   if (my_y==ncpu_y-1) by2=ny2+1; o_ny=g_ny+1
   if (my_z==ncpu_z-1 .and. g_nz>1) then
      bz2=nz2+1; o_nz=g_nz+1
   endif
endif


bdy_x1=g_bdy_x1
bdy_x2=g_bdy_x2
bdy_y1=g_bdy_y1
bdy_y2=g_bdy_y2
bdy_z1=g_bdy_z1
bdy_z2=g_bdy_z2

if (my_x/=0) bdy_x1=INTERNAL
if (my_x/=ncpu_x-1) bdy_x2=INTERNAL
if (my_y/=0) bdy_y1=INTERNAL
if (my_y/=ncpu_y-1) bdy_y2=INTERNAL
if (my_z/=0) bdy_z1=INTERNAL
if (my_z/=ncpu_z-1) bdy_z2=INTERNAL


if REALBOUNDARY(bdy_x1) then
   intx1=bx1+1
endif
if REALBOUNDARY(bdy_x2) then
   intx2=bx2-1
endif
if REALBOUNDARY(bdy_y1) then
   inty1=by1+1
endif
if REALBOUNDARY(bdy_y2) then
   inty2=by2-1
endif
if REALBOUNDARY(bdy_z1) then
   intz1=bz1+1
endif
if REALBOUNDARY(bdy_z2) then
   intz2=bz2-1
endif










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! global grid data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
delx = xscale/(o_nx-1)
dely = yscale/(o_ny-1)
if (o_nz==1) then
   delz=zscale
else
   delz = zscale/(o_nz-1)
endif

do i=1,o_nx
   g_xcord(i)=(i-1)*delx	
enddo

do j=1,o_ny
   g_ycord(j)=(j-1)*dely	
enddo
do k=1,o_nz
   g_zcord(k)=(k-1)*delz	
enddo

call fft_get_mcord(g_imcord,g_nx)
call fft_get_mcord(g_jmcord,g_ny)
call fft_get_mcord(g_kmcord,g_nz)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local grid data  3D decomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

imsine=0              
imsign=0   
imcord=0
do i=bx1,bx2
   l = i-nx1+1 + nslabx*my_x
   xcord(i)=g_xcord(l)
   imcord(i)=g_imcord(l)
   imsine(i)=l-1
   imsign(i)=sign(1,imcord(i))
   if (imcord(i)==0) imsign(i)=0
   if (imcord(i)==g_nx/2) imsign(i)=0
enddo


jmcord=0
jmsign=0
jmsine=0
do j=by1,by2
   l = j-ny1+1 + nslaby*my_y
   ycord(j)=g_ycord(l)
   jmsine(j)=l-1
   jmcord(j)=g_jmcord(l)
   jmsign(j)=sign(1,jmcord(j))
   if (jmcord(j)==0) jmsign(j)=0
   if (jmcord(j)==g_ny/2) jmsign(j)=0
enddo

kmcord=0
kmsign=0
kmsine=0
do k=bz1,bz2
   l = k-nz1+1 + nslabz*my_z
   zcord(k)=g_zcord(l)
   kmcord(k)=g_kmcord(l)
   kmsine(k)=l-1
   kmsign(k)=sign(1,kmcord(k))
   if (kmcord(k)==0) kmsign(k)=0
   if (kmcord(k)==g_nz/2) kmsign(k)=0
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! special data used for y-decomposition sine transform
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,nx_2dy
   k= -1 + nx1 + i + my_y*nx_2dy    ! nslabx/ncpu_y   
   if (k>nx2) then
      y_imsine(i)=0  ! this is possible if domain decompostion is not perfect
   else
      y_imsine(i)=imsine(k)
   endif
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local grid data, z-decomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=nx1,nx2
   z_imcord(i-nx1+1)=imcord(i)
   z_imsign(i-nx1+1)=imsign(i)
enddo

z_kmcord=g_kmcord
do k=1,g_nz
   z_kmsign(k)=sign(1,z_kmcord(k))
   if (z_kmcord(k)==0) z_kmsign(k)=0
   if (z_kmcord(k)==g_nz/2) z_kmsign(k)=0
enddo


do j=1,ny_2dz
   L = -1 + ny1 + j + my_z*ny_2dz      ! ncpy_z*ny_2dz = nslaby
   if (L>ny2) then
      z_jmcord(j)=0  ! possible, for not-perfect load balance cases
   else
      z_jmcord(j)=jmcord(L)
   endif
   z_jmsign(j)=sign(1,z_jmcord(j))
   if (z_jmcord(j)==0) z_jmsign(j)=0
   if (z_jmcord(j)==g_ny/2) z_jmsign(j)=0
enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local grid data, x-decomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do k=nz1,nz2
   x_kmcord(k-nz1+1)=kmcord(k)
   x_kmsign(k-nz1+1)=kmsign(k)
enddo

x_imcord=g_imcord
do i=1,g_nx
   x_imsign(i)=sign(1,x_imcord(i))
   if (x_imcord(i)==0) x_imsign(i)=0
   if (x_imcord(i)==g_nx/2) x_imsign(i)=0
enddo


do j=1,ny_2dx
   L = -1 + ny1 + j + my_x*ny_2dx      ! ncpy_z*ny_2dz = nslaby
   if (L>ny2) then
      x_jmcord(j)=0  ! possible, for not-perfect load balance cases
   else
      x_jmcord(j)=jmcord(L)
   endif
   x_jmsign(j)=sign(1,x_jmcord(j))
   if (x_jmcord(j)==0) x_jmsign(j)=0
   if (x_jmcord(j)==g_ny/2) x_jmsign(j)=0
enddo







end subroutine





subroutine init_input_file()
use params
use transpose
use fft_interface
use mpi
implicit none

!local variables
real*8 :: one=1,xfac,diff_2h,kmode
integer i,j,k,l,ierr
character(len=80) message
integer input_file_type
integer,external :: iargc
if (my_pe==io_pe) then
   !
   ! command line parameters
   ! ./dns -i -r -d rundir  runname  
   !
   !   -d <dirname>  put all output files in dirname/
   !
   !   -i <infile>   take input from infile instead of stdin
   !
   !   -t   enable LSF timelimit feature.  Only needed if you are 
   !        running under LSF and dont want the code to stop with 30min left
   !
   !   -r   use a restart file for the initial conditions
   !
   !   -b   enable bytewsap_input
   !
   !
   !   -ui  use UDM for input
   !   -uo  use UDM for output
   !
   !   -s   restart and output file are 2/3 dealiased spectral coefficieints,
   !        not the default grid point values
   !
   rundir="./"
   runname="temp"
   i=0
   do 
      i=i+1
      if (i>iargc()) exit
      call getarg(i,message)

      if (message(1:2)=="-r") then
         restart=1
      else if (message(1:2)=="-t") then
         enable_lsf_timelimit=1
      else if (message(1:2)=="-s") then
         rw_spec=.true.
      else if (message(1:3)=="-ui") then
         udm_input=.true.
      else if (message(1:3)=="-uo") then
         udm_output=.true.
      else if (message(1:4)=="-mio") then
         use_mpi_io=.true.
      else if (message(1:2)=="-d") then
         i=i+1
         if (i>iargc()) exit
         call getarg(i,rundir)
         j=len_trim(rundir)
         if (j>0) then
            if (rundir(j:j)/="/") then
               rundir(j+1:j+1)="/"
            endif
         endif
      else if (message(1:2)=="-i") then
         i=i+1
         if (i>iargc()) exit
         call getarg(i,inputfile)
         j=len_trim(inputfile)
      else if (message(1:1)/="-") then
         ! this must be the runname
         runname=message(1:len_trim(message))
      else
         print *,'Invalid options.'  
         print *,'./dns [-[t,r,b,s,ui,uo]]  [-i params.inp] [-d rundir] [runname] '
         print *,'-t  enable LSF time remaining check, if compiled in'
         print *,'-r  restart from rundir/restart.[uvw,h5] '
         print *,'-s  output spectral coefficients instead of grid data'
         print *,'-u1 restart file is 1st generation UDM'
         print *,'-ui restart file is UDM'
         print *,'-ui output to UDM'
         print *,'-mio use MPI-IO'
      endif
   enddo
   print *,'Run name:         ',runname(1:len_trim(runname))
   print *,'output directory: ',rundir(1:len_trim(rundir))

   if (len_trim(inputfile)>0) then
      open(5,file=inputfile)
      print *,'reading input from file: ',inputfile	
   endif	
   print *,'Enter input file type: '
   read(5,*) input_file_type
   print *,'input_file_type=',input_file_type	

   if (input_file_type==0) then
      call read_type0()
   else if (input_file_type==2) then
      call read_type2()
   else if (input_file_type==3) then
      call read_type3()
   else
      ! if you add a new variable and new input type, be sure to add
      ! it to the MPI broadcast!
      call abort("bad input file")
   endif

   if (rw_spec .and. dealias/=1) then
      call abort("Error: -s option only works with fft-dealias method")
   endif

endif




#ifdef USE_MPI
call MPI_bcast(runname,80,MPI_CHARACTER,io_pe,comm_3d ,ierr)
call MPI_bcast(rundir,80,MPI_CHARACTER,io_pe,comm_3d ,ierr)
call MPI_bcast(mu,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(rw_spec,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(udm_input,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(udm_output,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(equations,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(mu_hyper,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(dealias,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(numerical_method,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(time_final,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(cfl_adv,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(cfl_vis,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(delt_min,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(delt_max,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(restart_dt,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(diag_dt,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(screen_dt,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(output_dt,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(restart,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(enable_lsf_timelimit,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(init_cond,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(init_cond_subtype,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(forcing_type,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(g_bdy_x1,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(g_bdy_x2,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(g_bdy_y1,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(g_bdy_y2,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(g_bdy_z1,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(g_bdy_z2,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(compute_struct,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(alpha_value,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(smagorinsky,1,MPI_REAL8,io_pe,comm_3d ,ierr)


call MPI_bcast(ncustom,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
if (.not. allocated(custom)) allocate(custom(ncustom))
call MPI_bcast(custom,ncustom,MPI_REAL8,io_pe,comm_3d,ierr)
call MPI_Barrier(comm_3d,ierr)

#endif




! Diffusion of mode k:
! d/dt (.5 u**2)  = - mu 2pi 2pi k k u**2
!
! d/dt (KE) / KE = - 2 mu 2pi 2pi k**2,  
! USE k = the dealiased highest mode, k = g_nx/3
!
! if mu is chosen for a given value of diff_2h,
! 
!    mu = diff_2h / (2 2pi 2pi k**2)
! 


write(message,'(a,e10.4)') 'Diffusion coefficient mu=',mu
call print_message(message)


if (mu_hyper>1) then
write(message,'(a,e10.4)') 'Hyper diffusion coefficient mu=',mu
call print_message(message)
kmode=1
xfac = 2* (2*pi*2*pi*kmode**2)**4
diff_2h =  mu*xfac
write(message,'(a,f5.0,a,f8.2)') 'Diffusion d/dt(KE)/KE on mode k = ',kmode,': ',diff_2h
call print_message(message)

kmode = sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
if (ndim==2) kmode = sqrt( (g_nx**2 + g_ny**2 )/4.0)
xfac = 2*(2*pi*2*pi*kmode**2)**4
diff_2h =  mu*xfac
write(message,'(a,f5.0,a,f8.2)') 'Diffusion d/dt(KE)/KE on mode k = ',kmode,': ',diff_2h
call print_message(message)

else

if (mu>0) then
kmode=1
xfac = 2*2*pi*2*pi*kmode**2
diff_2h =  mu*xfac
write(message,'(a,f5.0,a,f8.2)') 'Diffusion d/dt(KE)/KE on mode k = ',kmode,': ',diff_2h
call print_message(message)

kmode = sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
if (ndim==2) kmode = sqrt( (g_nx**2 + g_ny**2 )/4.0)
xfac = 2*2*pi*2*pi*kmode**2
diff_2h =  mu*xfac
write(message,'(a,f5.0,a,f8.2)') 'Diffusion d/dt(KE)/KE on mode k = ',kmode,': ',diff_2h
call print_message(message)
endif

endif

end subroutine














subroutine read_type0
use params
implicit none

!local variables
character(len=80) message
character(len=20) sdata
real*8 rvalue
integer i


read(5,'(a)') message
do i=1,80
   if (message(i:i)==' ') then
      runname=message(1:i-1)
      exit
   endif
enddo


read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='value') then
   mu=rvalue
else 
   call abort("only viscosity type 'value' supported for input format=0")
endif

read(5,'(a12)') sdata
if (sdata=='fft') then
   numerical_method=FOURIER
   dealias=0
else if (sdata=='fft-dealias') then
   numerical_method=FOURIER
   dealias=1
else if (sdata=='fft-sphere') then
   numerical_method=FOURIER
   dealias=2
else if (sdata=='4th') then
   numerical_method=FOURTH_ORDER
   dealias=0
else
   call abort("only 'fft' derivative method supported")
endif


read(5,'(a12)') sdata
if (sdata=='periodic') then
else
   call abort("only 'perodic' b.c. supported")
endif

read(5,'(a12)') sdata
if (sdata=='periodic') then
else
   call abort("only 'perodic' b.c. supported")
endif

read(5,'(a12)') sdata
if (sdata=='periodic') then
else
   call abort("only 'perodic' b.c. supported")
endif

read(5,*) time_final
read(5,*) cfl_adv
read(5,*) cfl_vis
read(5,*) delt_min
read(5,*) delt_max
read(5,*) restart_dt
read(5,*) diag_dt
read(5,*) screen_dt
read(5,*) output_dt
read(5,*) ncustom
allocate(custom(ncustom))
do i=1,ncustom
   read(5,*) custom(i)
enddo

init_cond=1           ! KH analytic - default initial condition
init_cond_subtype=0   ! default parameters
forcing_type=0   ! no forcing
compute_struct=0
alpha_value=0


end subroutine








subroutine read_type2
use params
implicit none

!local variables
character(len=80) message
character(len=20) sdata
real*8 rvalue,xfac,kmode
integer i


read(5,'(a12)') sdata
print *,'initial condition: ',sdata
if (sdata=='KH-blob') then
   init_cond=0
else if (sdata=='KH-anal') then
   init_cond=1
else if (sdata=='iso12') then
   init_cond=2
else if (sdata=='sht') then
   init_cond=3
else if (sdata=='vxpair') then
   init_cond=4
else if (sdata=='iso12e') then
   init_cond=5
else 
   print *,'value = >>',sdata,'<<'
   call abort("invalid initial condtion specified on line 3 on input file")
endif

read(5,*) init_cond_subtype


read(5,'(a12)') sdata
print *,'forcing: ',sdata
if (sdata=='none') then
   forcing_type=0
else if (sdata=='iso12') then
   forcing_type=1
else if (sdata=='iso12w') then
   forcing_type=2
else 
   call abort("invalid forcing type specified on line 4 on input file")
endif


read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='value') then
   mu=rvalue
else if (sdata=='kediff') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2*2*pi*2*pi*kmode**2
   ! diff_2h =  mu*xfac
   mu = rvalue/xfac
else if (sdata=='hyper') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2* (2*pi*kmode)**8
   mu = rvalue/xfac
   mu_hyper = 4
else 
   call abort("non supported viscosity type")
endif

read(5,*) alpha_value
read(5,*) compute_struct


read(5,'(a12)') sdata
print *,'method: ',sdata
if (sdata=='fft') then
   numerical_method=FOURIER
   dealias=0
else if (sdata=='fft-dealias') then
   numerical_method=FOURIER
   dealias=1
else if (sdata=='fft-sphere') then
   numerical_method=FOURIER
   dealias=2
else if (sdata=='4th') then
   numerical_method=FOURTH_ORDER
   dealias=0
else
   print *,'value=',sdata
   call abort("only 'fft' derivative method supported")
endif


read(5,'(a12)') sdata
print *,'b.c. x-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_x1=PERIODIC
   g_bdy_x2=PERIODIC
else if (sdata=='in0/onesided') then
   g_bdy_x1=INFLOW0_ONESIDED
   g_bdy_x2=INFLOW0_ONESIDED
else if (sdata=='in0') then
   g_bdy_x1=INFLOW0
   g_bdy_x2=INFLOW0
else
   call abort("x boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. y-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_y1=PERIODIC
   g_bdy_y2=PERIODIC
else if (sdata=='reflect') then
   g_bdy_y1=REFLECT
   g_bdy_y2=REFLECT
else if (sdata=='custom0') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0_ONESIDED
else if (sdata=='custom1') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0
else
   call abort("y boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. z-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_z1=PERIODIC
   g_bdy_z2=PERIODIC
else
   call abort("only 'perodic' b.c. supported in z-direction")
endif

read(5,*) time_final
read(5,*) cfl_adv
read(5,*) cfl_vis
read(5,*) delt_min
read(5,*) delt_max
read(5,*) restart_dt
read(5,*) diag_dt
read(5,*) screen_dt
read(5,*) output_dt
read(5,*) ncustom
allocate(custom(ncustom))
do i=1,ncustom
   read(5,*) custom(i)
enddo


if (numerical_method==FOURIER) then
   ! make sure boundary conditons are periodic:
   if (g_bdy_x1/=PERIODIC .or. g_bdy_x2/=PERIODIC .or. &
       g_bdy_y1/=PERIODIC .or. g_bdy_y2/=PERIODIC .or. &
       g_bdy_z1/=PERIODIC .or. g_bdy_z2/=PERIODIC) then
      call abort("FOURIER method requires all boundary conditions be PERIODIC")
   endif
endif


end subroutine








subroutine read_type3
use params
implicit none

!local variables
character(len=80) message
character(len=20) sdata
real*8 rvalue,xfac,kmode
integer i

read(5,'(a12)') sdata
write(*,'(a,a)') 'equations: ',sdata
if (sdata=='ns_uvw') then
   equations=NS_UVW
else if (sdata=='ns_psivor') then
   equations=NS_PSIVOR
else if (sdata=='shallow') then
   equations=SHALLOW
else 
   print *,'value = >>',sdata,'<<'
   call abort("invalid equations specified")
endif

read(5,'(a12)') sdata
print *,'initial condition: ',sdata
if (sdata=='KH-blob') then
   init_cond=0
else if (sdata=='KH-anal') then
   init_cond=1
else if (sdata=='iso12') then
   init_cond=2
else if (sdata=='sht') then
   init_cond=3
else if (sdata=='vxpair') then
   init_cond=4
else if (sdata=='iso12e') then
   init_cond=5
else if (sdata=='zero') then
   init_cond=6
else if (sdata=='decay2048') then
   init_cond=7
else if (sdata=='decay2048_e') then
   init_cond=8
else if (sdata=='decay2048_s') then
   init_cond=9
else 
   print *,'value = >>',sdata,'<<'
   call abort("invalid initial condtion specified on line 3 on input file")
endif

read(5,*) init_cond_subtype


read(5,'(a12)') sdata
print *,'forcing: ',sdata
if (sdata=='none') then
   forcing_type=0
else if (sdata=='iso12') then
   forcing_type=1
else if (sdata=='iso12w') then
   forcing_type=2
else if (sdata=='iso') then
   forcing_type=3
else if (sdata=='iso23w') then
   forcing_type=4
else 
   call abort("invalid forcing type specified on line 4 on input file")
endif


read(5,'(a12)') sdata
read(5,*) rvalue
if (sdata=='value') then
   mu=rvalue
else if (sdata=='kediff') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2*2*pi*2*pi*kmode**2
   ! diff_2h =  mu*xfac
   mu = rvalue/xfac
else if (sdata=='hyper') then
   if (ndim==2) then
      kmode=sqrt( (g_nx**2 + g_ny**2 )/4.0)
   else
      kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   endif
   xfac = 2* (2*pi*kmode)**8
   mu = rvalue/xfac
   mu_hyper = 4
else 
   call abort("non supported viscosity type")
endif

read(5,*) alpha_value
read(5,*) smagorinsky
read(5,*) compute_struct


read(5,'(a12)') sdata
print *,'method: ',sdata
if (sdata=='fft') then
   numerical_method=FOURIER
   dealias=0
else if (sdata=='fft-dealias') then
   numerical_method=FOURIER
   dealias=1
else if (sdata=='fft-sphere') then
   numerical_method=FOURIER
   dealias=2
else if (sdata=='4th') then
   numerical_method=FOURTH_ORDER
   dealias=0
else
   print *,'value=',sdata
   call abort("only 'fft' derivative method supported")
endif


read(5,'(a12)') sdata
print *,'b.c. x-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_x1=PERIODIC
   g_bdy_x2=PERIODIC
else if (sdata=='in0/onesided') then
   g_bdy_x1=INFLOW0_ONESIDED
   g_bdy_x2=INFLOW0_ONESIDED
else if (sdata=='in0') then
   g_bdy_x1=INFLOW0
   g_bdy_x2=INFLOW0
else
   call abort("x boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. y-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_y1=PERIODIC
   g_bdy_y2=PERIODIC
else if (sdata=='reflect') then
   g_bdy_y1=REFLECT
   g_bdy_y2=REFLECT
else if (sdata=='custom0') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0_ONESIDED
else if (sdata=='custom1') then
   g_bdy_y1=REFLECT_ODD
   g_bdy_y2=INFLOW0
else
   call abort("y boundary condition not supported")
endif

read(5,'(a12)') sdata
print *,'b.c. z-direction: ',sdata
if (sdata=='periodic') then
   g_bdy_z1=PERIODIC
   g_bdy_z2=PERIODIC
else
   call abort("only 'perodic' b.c. supported in z-direction")
endif

read(5,*) time_final
read(5,*) cfl_adv
read(5,*) cfl_vis
read(5,*) delt_min
read(5,*) delt_max
read(5,*) restart_dt
read(5,*) diag_dt
read(5,*) screen_dt
read(5,*) output_dt
read(5,*) ncustom
allocate(custom(ncustom))
do i=1,ncustom
   read(5,*) custom(i)
enddo


if (numerical_method==FOURIER) then
   ! make sure boundary conditons are periodic:
   if (g_bdy_x1/=PERIODIC .or. g_bdy_x2/=PERIODIC .or. &
       g_bdy_y1/=PERIODIC .or. g_bdy_y2/=PERIODIC .or. &
       g_bdy_z1/=PERIODIC .or. g_bdy_z2/=PERIODIC) then
      call abort("FOURIER method requires all boundary conditions be PERIODIC")
   endif
endif


end subroutine








