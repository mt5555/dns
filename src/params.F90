#include "macros.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  module containing all global parameters like mesh data
!
! To create a valid file:  
! Set nx1,nx2,ny1,ny2,nz1,nz2   and nx,ny,nz
! pick l,m,n, and then set:
!
! ncpu_x*l=(nz2-nz1+1)
! ncpu_y*m=(nx2-nx1+1)
! ncpu_z*n=(ny2-ny1+1)
!
!
! params_init() will set the global grid by:
! g_nx=ncpu_x*(nx2-nx1+1) = l*nslabz*nslabx
! g_ny=ncpu_y*(ny2-ny1+1) = m*nslabx*nslaby
! g_nz=ncpu_z*(nz2-nz1+1) = n*nslaby*nslabz
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module params

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8  :: mu=0   !viscosity
real*8  :: pi,pi2_squared



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! global dimensions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
integer :: g_nx,g_ny,g_nz    ! dimension of global grid (unique data)
integer :: g_nx2,g_ny2,g_nz2 ! dimension used by fft
                             ! (because FFT99 requires 2 extra points)
 


! mesh dimensions on a single processor
#include "user_params.h"
!integer,parameter :: n_var=3                  ! number of prognostic variables
!integer,parameter :: nx=18,ny=18,nz=18         ! dimension of grid & data
!integer,parameter :: nx1=1,nx2=16              ! upper and lower bounds of non-ghost data
!integer,parameter :: ny1=1,ny2=16             ! upper and lower bounds of non-ghost data
!integer,parameter :: nz1=1,nz2=16             ! upper and lower bounds of non-ghost data

! number of actual data points
integer,parameter :: nslabx=nx2-nx1+1
integer,parameter :: nslaby=ny2-ny1+1
integer,parameter :: nslabz=nz2-nz1+1

! mesh coordinates
real*8 :: xcord(nx),delx
real*8 :: ycord(ny),dely
real*8 :: zcord(nz),delz
real*8,allocatable :: g_xcord(:)  
real*8,allocatable :: g_ycord(:)  
real*8,allocatable :: g_zcord(:)  
integer :: imcord(nx),jmcord(ny),kmcord(nz)  ! fft modes local
integer,allocatable :: g_imcord(:)           ! fft modes global
integer,allocatable :: g_jmcord(:)  
integer,allocatable :: g_kmcord(:)  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! time stepping
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8  :: delt
integer :: itime_max=1000   ! number of time steps
integer :: itime_output=10  ! output every itime_output time steps


! set non-zero to cause time stepping loop to exit
! error_code = 1   advection CFL violated
! error_code = 2   sound speed CFL violated
! error_code = 3   diffusion CFL violated
integer :: error_code 





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! parallel decompositions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!integer :: ncpu_x=1,ncpu_y=1,ncpu_z=1
integer :: ncpu_z
integer :: nx_2d,ny_2d,nz_2d
integer :: io_pe
integer :: my_world_pe,my_pe,mpicoords(3),mpidims(3)
integer :: initial_live_procs
integer :: comm_3d                 ! the MPI cartesian communcator



!
! cpu decomposition, 3D:   ncpu_x * ncpu_y * ncpu_z
! ncpu_x = g_nx/(nx2-nx1+1)
! ncpu_y = g_ny/(ny2-ny1+1)
! ncpu_z = g_nz/(nz2-nz1+1)
!
! 2D cpu decomposition,  GRID                        CPU
! (+2 needed because of FFT99)
! x-direction:  g_nx2,(nslaby*nslabz)/ncpu_x     ncpu_y*(ncpu_x*ncpu_z)
! y-direction:  g_ny2,(nslabz*nslabx)/ncpu_y     ncpu_z*(ncpu_y*ncpu_x)
! z-direction:  g_nz2,(nslabx*nslaby)/ncpu_z     ncpu_x*(ncpu_z*ncpu_y)
!
!
! for simplicity, we require:
!    ncpu_x divides (nz2-nz1+1)   nz_2d=(nz2-nz1+1)/ncpu_x
!    ncpu_y divides (nx2-nx1+1)   nx_2d=(nx2-nx1+1)/ncpu_y
!    ncpu_z divides (ny2-ny1+1)   ny_2d=(ny2-ny1+1)/ncpu_z
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




contains

subroutine params_init
character*80 message
integer :: fail=0

real*8 :: one=1
pi=4*atan(one)
pi2_squared=4*pi*pi

g_nx=(nx2-nx1+1)*ncpu_x
g_ny=(ny2-ny1+1)*ncpu_y
g_nz=(nz2-nz1+1)*ncpu_z
g_nx2=g_nx+2
g_ny2=g_ny+2
g_nz2=g_nz+2
if (g_nz==1) g_nz2=1


! these values must divide with no remainder:
nz_2d=(nz2-nz1+1)/ncpu_x
nx_2d=(nx2-nx1+1)/ncpu_y
ny_2d=(ny2-ny1+1)/ncpu_z
if (0/=mod(nz2-nz1+1,ncpu_x)) then
   fail=1
   call print_message("ncpu_x does not divide nz");
endif
if (0/=mod(nx2-nx1+1,ncpu_y)) then
   fail=1
   call print_message("ncpu_y does not divide nx");
endif
if (0/=mod(ny2-ny1+1,ncpu_z)) then
   fail=1
   call print_message("ncpu_z does not divide nz");
endif



! memory contraint: 2D decomposition should fit in a 3D decomposition array
! check if: (g_nz2)*nslabx*ny_2d <= nx*ny*nz)
!
if ((g_nx2)*nz_2d*real(nslaby) > nx*nz*real(ny) )  then
   fail=1
   call print_message("insufficient storage.");
   write(message,'(a,3i6,a,f10.0)') "nx,ny,nz=",nx,ny,nz,"   nx*ny*nz=",nx*ny*real(nz)
   call print_message(message)	
   write(message,'(a,f10.0)') "storage needed for 2D x-decomposition: ", &
     (g_nx2)*nz_2d*real(nslaby)
   call print_message(message)	

endif
if ( (g_ny2)*nx_2d*real(nslabz) > nx*ny*real(nz)) then
   fail=1
   call print_message("insufficient storage.");
   write(message,'(a,3i6,a,f10.0)') "nx,ny,nz=",nx,ny,nz,"   nx*ny*nz=",nx*ny*real(nz)
   call print_message(message)	
   write(message,'(a,f10.0)') "storage needed for 2D y-decomposition: ", &
     (g_ny2)*nx_2d*real(nslabz)
   call print_message(message)	

endif

if ( (g_nz2)*ny_2d*real(nslabx) >  ny*nz*real(nx)) then
   fail=1
   call print_message("insufficient storage.");
   write(message,'(a,3i6,a,f10.0)') "nx,ny,nz=",nx,ny,nz,"   nx*ny*nz=",nx*ny*real(nz)
   call print_message(message)	
   write(message,'(a,f10.0)') "storage needed for 2D z-decomposition: ", &
     (g_nz2)*ny_2d*real(nslabx)
   call print_message(message)	
endif






if ( nx2>nx) then
   fail=1
   call print_message("nx is too small. nx must be >=nx2")	
endif
if ( ny2>ny) then
   fail=1
   call print_message("ny is too small. ny must be >=ny2")	
endif
if ( nz2>nz) then
   fail=1
   call print_message("nz is too small. nz must be >=nz2")	
endif

if (fail/=0) call abort("params.F90 dimension settings failure")

allocate(g_xcord(g_nx))
allocate(g_ycord(g_ny))
allocate(g_zcord(g_nz))
allocate(g_imcord(g_nx))
allocate(g_jmcord(g_ny))
allocate(g_kmcord(g_nz))


end subroutine
                      

end ! module mod_params








