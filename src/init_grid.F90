#include "macros.h"
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
use mpi
use fft_interface
implicit none

!local variables
real*8 :: one=1,xfac,diff_2h,kmode
integer i,j,k,l,ierr
character(len=80) message
integer input_file_type


call params_init()
call fft_interface_init()

if (my_pe==io_pe) then
   print *,'Enter input file type: '
   read(*,*) input_file_type
   if (input_file_type==0) then
      print *,'Calling read_type0()'
      call read_type0()
   else if (input_file_type==1) then
      print *,'Calling read_type1()'
      call read_type1()
   else
      call abort("bad input file")
   endif
endif

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

kmode=1
xfac = 2*2*pi*2*pi*kmode**2
diff_2h =  mu*xfac
write(message,'(a,f5.0,a,f8.2)') 'Diffusion d/dt(KE)/KE on mode k = ',kmode,': ',diff_2h
call print_message(message)

kmode = sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
xfac = 2*2*pi*2*pi*kmode**2
diff_2h =  mu*xfac
write(message,'(a,f5.0,a,f8.2)') 'Diffusion d/dt(KE)/KE on mode k = ',kmode,': ',diff_2h
call print_message(message)



#ifdef USE_MPI
call MPI_bcast(runname,80,MPI_CHARACTER,io_pe,comm_3d ,ierr)
call MPI_bcast(mu,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(dealias,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(time_final,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(cfl_adv,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(cfl_vis,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(delt_min,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(delt_max,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(restart_dt,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(diag_dt,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(screen_dt,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(output_dt,1,MPI_REAL8,io_pe,comm_3d ,ierr)
call MPI_bcast(ncustom,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(init_cond,1,MPI_INTEGER,io_pe,comm_3d ,ierr)
call MPI_bcast(forcing_type,1,MPI_INTEGER,io_pe,comm_3d ,ierr)


if (.not. allocated(custom)) allocate(custom(ncustom))
call MPI_bcast(custom,ncustom,MPI_REAL8,io_pe,comm_3d,ierr)
call MPI_Barrier(comm_3d,ierr)

#endif



! periodic FFT case:  for output, we include the point at x=1 (same as x=0)
o_nx=g_nx+1
o_ny=g_ny+1
o_nz=g_nz+1
if (g_nz==1) o_nz=1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! global grid data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
delx = one/g_nx
dely = one/g_ny
delz = one/g_nz

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

do i=nx1,nx2
   l = i-nx1+1 + nslabx*my_x
   xcord(i)=g_xcord(l)
   imcord(i)=g_imcord(l)
   imsign(i)=sign(1,imcord(i))
   if (imcord(i)==0) imsign(i)=0
   if (imcord(i)==g_nx/2) imsign(i)=0
enddo
do j=ny1,ny2
   l = j-ny1+1 + nslaby*my_y
   ycord(j)=g_ycord(l)
   jmcord(j)=g_jmcord(l)
   jmsign(j)=sign(1,jmcord(j))
   if (jmcord(j)==0) jmsign(j)=0
   if (jmcord(j)==g_ny/2) jmsign(j)=0
enddo
do k=nz1,nz2
   l = k-nz1+1 + nslabz*my_z
   zcord(k)=g_zcord(l)
   kmcord(k)=g_kmcord(l)
   kmsign(k)=sign(1,kmcord(k))
   if (kmcord(k)==0) kmsign(k)=0
   if (kmcord(k)==g_nz/2) kmsign(k)=0
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
   z_jmcord(j)=jmcord(L)
   z_jmsign(j)=sign(1,z_jmcord(j))
   if (z_jmcord(j)==0) z_jmsign(j)=0
   if (z_jmcord(j)==g_ny/2) z_jmsign(j)=0
enddo




write(message,'(a,i6,a,i6,a,i6)') "Global grid: ",g_nx," x",g_ny," x",g_nz
call print_message(message)
write(message,'(a,i6,a,i6,a,i6)') "Local grid (with padding): ",nx," x",ny," x",nz
call print_message(message)



end subroutine



subroutine read_type0
use params
implicit none

!local variables
character(len=80) message
character(len=20) sdata
real*8 rvalue
integer i


read(*,'(a)') message
do i=1,80
   if (message(i:i)==' ') then
      runname=message(1:i-1)
      exit
   endif
enddo


read(*,'(a12)') sdata
read(*,*) rvalue
if (sdata=='value') then
   mu=rvalue
else 
   call abort("only viscosity type 'value' supported for input format=0")
endif

read(*,'(a12)') sdata
if (sdata=='fft') then
   dealias=.false.
else if (sdata=='fft-dealias') then
   dealias=.true.
else
   call abort("only 'fft' derivative method supported")
endif


read(*,'(a12)') sdata
if (sdata=='periodic') then
else
   call abort("only 'perodic' b.c. supported")
endif

read(*,'(a12)') sdata
if (sdata=='periodic') then
else
   call abort("only 'perodic' b.c. supported")
endif

read(*,'(a12)') sdata
if (sdata=='periodic') then
else
   call abort("only 'perodic' b.c. supported")
endif

read(*,*) time_final
read(*,*) cfl_adv
read(*,*) cfl_vis
read(*,*) delt_min
read(*,*) delt_max
read(*,*) restart_dt
read(*,*) diag_dt
read(*,*) screen_dt
read(*,*) output_dt
read(*,*) ncustom
allocate(custom(ncustom))
do i=1,ncustom
   read(*,*) custom(i)
enddo

init_cond=1      ! KH analytic - default initial condition
forcing_type=0   ! no forcing

end subroutine








subroutine read_type1
use params
implicit none

!local variables
character(len=80) message
character(len=20) sdata
real*8 rvalue,xfac,kmode
integer i


read(*,'(a)') message
do i=1,80
   if (message(i:i)==' ') then
      runname=message(1:i-1)
      exit
   endif
enddo

read(*,'(a12)') sdata
if (sdata=='KH-blob') then
   init_cond=0
else if (sdata=='KH-anal') then
   init_cond=1
else if (sdata=='iso12') then
   init_cond=2
else 
   call abort("invalid initial condtion specified on line 3 on input file")
endif


read(*,'(a12)') sdata
if (sdata=='none') then
   forcing_type=0
else if (sdata=='iso12') then
   forcing_type=1
else 
   call abort("invalid forcing type specified on line 4 on input file")
endif

read(*,*) struct_nx
if (struct_nx==1) call abort("struct_nx = 1 in input file not allowed")
read(*,*) struct_ny
if (struct_ny==1) call abort("struct_ny = 1 in input file not allowed")
read(*,*) struct_nz


read(*,'(a12)') sdata
read(*,*) rvalue
if (sdata=='value') then
   mu=rvalue
else if (sdata=='kediff') then
   kmode=sqrt( (g_nx**2 + g_ny**2 + g_nz**2)/9.0)
   xfac = 2*2*pi*2*pi*kmode**2
   ! diff_2h =  mu*xfac
   mu = rvalue/xfac
else 
   call abort("only viscosity type 'value' supported")
endif

read(*,'(a12)') sdata
if (sdata=='fft') then
   dealias=.false.
else if (sdata=='fft-dealias') then
   dealias=.true.
else
   call abort("only 'fft' derivative method supported")
endif


read(*,'(a12)') sdata
if (sdata=='periodic') then
else
   call abort("only 'perodic' b.c. supported")
endif

read(*,'(a12)') sdata
if (sdata=='periodic') then
else
   call abort("only 'perodic' b.c. supported")
endif

read(*,'(a12)') sdata
if (sdata=='periodic') then
else
   call abort("only 'perodic' b.c. supported")
endif

read(*,*) time_final
read(*,*) cfl_adv
read(*,*) cfl_vis
read(*,*) delt_min
read(*,*) delt_max
read(*,*) restart_dt
read(*,*) diag_dt
read(*,*) screen_dt
read(*,*) output_dt
read(*,*) ncustom
allocate(custom(ncustom))
do i=1,ncustom
   read(*,*) custom(i)
enddo

end subroutine








