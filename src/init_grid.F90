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
use fft_interface
implicit none

!local variables
real*8 :: one=1
integer i,j,k,l
character*80 message
integer input_file_type


call params_init()
call fft_interface_init()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! global grid data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
delx = one/g_nx
dely = one/g_ny
delz = one/g_nz

do i=1,g_nx
   g_xcord(i)=(i-1)*delx	
enddo
do j=1,g_ny
   g_ycord(j)=(j-1)*dely	
enddo
do k=1,g_nz
   g_zcord(k)=(k-1)*delz	
enddo

call fft_get_mcord(g_imcord,g_nx)
call fft_get_mcord(g_jmcord,g_ny)
call fft_get_mcord(g_kmcord,g_nz)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local grid data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call init_mpi_comm3d()

do i=nx1,nx2
   l = i-nx1+1 + nslabx*my_x
   xcord(i)=g_xcord(l)
   imcord(i)=g_imcord(l)
enddo
do j=ny1,ny2
   l = j-ny1+1 + nslaby*my_y
   ycord(j)=g_ycord(l)
   jmcord(j)=g_jmcord(l)
enddo
do k=nz1,nz2
   l = k-nz1+1 + nslabz*my_z
   zcord(k)=g_zcord(l)
   kmcord(k)=g_kmcord(l)
enddo

write(message,'(a,i6,a,i6,a,i6)') "Global grid: ",g_nx," x",g_ny," x",g_nz
call print_message(message)
write(message,'(a,i6,a,i6,a,i6)') "Local grid (with padding): ",nx," x",ny," x",nz
call print_message(message)


if (my_pe==io_pe) then
   read(*,*) input_file_type
   if (input_file_type==0) then
      call read_type0()
   else
      call abort("bad input file")
   endif
   !  u_t = mu * (2 pi k)^2
   print *,'Diffusion damping time on smallest scale: ',&
     1/(mu*pi*pi*(g_nx**2 + g_ny**2 + g_nz**2))
endif

#ifdef MPI
call MPI_BROADCAST
mu
time_final
cfl_adv
cfl_vis
delt_min
delt_max
restart_dt
diag_dt
screen_dt
output_dt
ncustom
custom(ncustom)

#endif


end subroutine



subroutine read_type0
use params
implicit none

!local variables
character*80 message
character*20 sdata
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
