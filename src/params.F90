#include "macros.h"
#include "transpose.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  module containing all global parameters like mesh data
!
! To create a valid file:  
! pick ncpu_x,ncpu_y,ncpu_z	
! Then l,m,n, and set:
!
! #ifdef TRANSPOSE_X_SPLIT_Z   (NO LONGER SUPPORTED)
!   
!    ncpu_x*l=nslabz
!    ncpu_y*m=nslabx
!    ncpu_z*n=nslaby
!
! #ifdef TRANSPOSE_X_SPLIT_Y        (usefull if nslabz=1)
!    ncpu_x*l=nslaby
!    ncpu_y*m=nslabx
!    ncpu_z*n=nslaby  <==>    n=l*ncpu_x/ncpu_z
!
! from nslab*, pick suitable values of nx1,nx2,ny1,ny2,nz1,nz2   and nx,ny,nz)
!
! params_init() will set the global grid by:
! g_nx=ncpu_x*nslabx
! g_ny=ncpu_y*nslaby
! g_nz=ncpu_z*nslabz
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module params

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8  :: mu=0           !viscosity
real*8  :: alpha_value=.1      !for the alpha model  
integer,parameter :: r8kind=kind(mu)
logical :: dealias       
character(len=80) :: runname
real*8  :: pi,pi2,pi2_squared



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initial condition and forcing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: init_cond      ! 0 = KH vortex blob
                          ! 1 = KH analytic
                          ! 2 = random, isotropic, .5 < wave no. < 2.5
                          !     E(1)=1, E(2)=2**(-5/3)
integer :: forcing_type   ! 0 = none
                          ! 1 = relax back to E(1)=1, E(2)=2**(-5/3)
                          !     can only be used by the z-decomp model!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! global dimensions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
integer :: g_nx,g_ny,g_nz    ! dimension of global grid (unique data)
integer :: o_nx,o_ny,o_nz    ! dimensions of plotting output data
                             ! For periodic FFT case, o_nx=g_nx+1 because we do not
                             ! store the point at x=1.  for 4th order case, we 
                             ! may store this point.  
                             
integer :: g_nx2,g_ny2,g_nz2 ! dimension used by fft
                             ! (because FFT99 requires 2 extra points)
 


! mesh dimensions on a single processor
#include "params.h"
!integer,parameter :: n_var=3                  ! number of prognostic variables
!integer,parameter :: nxd=18,nyd=18,nzd=18         ! dimension of grid & data
!integer,parameter :: nx1=1,nx2=16              ! upper and lower bounds of non-ghost data
!integer,parameter :: ny1=1,ny2=16             ! upper and lower bounds of non-ghost data
!integer,parameter :: nz1=1,nz2=16             ! upper and lower bounds of non-ghost data

! NOTE: if these are parameters, then all automatic arrays
! in subroutines will be allocated at run time.  This doubles
! the amount of memory needed:
!integer,parameter :: nx=nxd,ny=nyd,nz=nzd      ! dimension of grid & data

! NOTE: if these are NOT parameters, then all automatic arrays
! are placed on the stack, and thus dont take up memory between calls.
! this halves the amount of memory needed.  
! any performance penaulty?
integer           :: nx=nxd,ny=nyd,nz=nzd      ! dimension of grid & data


! number of actual data points
integer,parameter :: nslabx=nx2-nx1+1
integer,parameter :: nslaby=ny2-ny1+1
integer,parameter :: nslabz=nz2-nz1+1

! mesh coordinates
real*8 :: xcord(nxd),delx
real*8 :: ycord(nyd),dely
real*8 :: zcord(nzd),delz
real*8,allocatable :: g_xcord(:)  
real*8,allocatable :: g_ycord(:)  
real*8,allocatable :: g_zcord(:)  

! fft modes, global
integer,allocatable :: g_imcord(:)
integer,allocatable :: g_jmcord(:)  
integer,allocatable :: g_kmcord(:)  

! fft modes, local 3D decompostion
integer :: imcord(nxd),jmcord(nyd),kmcord(nzd)  ! fft modes local
integer :: imsign(nxd),jmsign(nyd),kmsign(nzd)  ! fft modes local

! fft modes, local z-decompostion
integer,allocatable :: z_imcord(:),Z_jmcord(:),z_kmcord(:)  ! fft modes local
integer,allocatable :: z_imsign(:),z_jmsign(:),z_kmsign(:)  ! fft modes local


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! time stepping
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8  :: delt
real*8  :: cfl_adv = 1.2
real*8  :: cfl_vis = 1.2
real*8  :: delt_min = 0
real*8  :: delt_max = 1
real*8  :: time_final = 1 
integer :: struct_nx=10          ! compute structure functions 
integer :: struct_ny=10          ! every 10 time steps
integer :: struct_nz=10          !


real*8 :: output_dt = 0    ! netcdf output for plotting
integer :: ncustom =0
real*8, allocatable :: custom(:)
real*8 :: diag_dt = 0       ! diagnostics
real*8 :: screen_dt = 0     ! screen output
real*8 :: restart_dt = 0    ! restart 

! set to > 0 to cause time stepping loop to exit.  if more than one cpu sets
! error_code > 0, the largest value will be reported.  
!
! error_code = 1   u > 1000
! error_code = 2   
! error_code = 3   
!
integer :: error_code =0



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! scalar quantities of current state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer,parameter :: nints=8
real*8 :: ints(nints)=0,maxs(nints)=0
real*8 :: ints_timeU,ints_timeDU

!
! Quantities involving only Q are computed in the timestep
! routine, at the current time = ints_timeU:  
!
! ints(1) = ke
! maxs(1,2,3) = max U,V,W
! maxs(4) = max (U/delx + V/dely + W/delz)  used for CFL
!
! Quantities involving derivatives are computed when the RHS is computed
! and thus the data is at the prevous timestep = ints_timeDU
! 
! ints(2) = ke dissapation from forcing
! ints(3) = ke dissapation from diffusion
! ints(4) = z-component of vorticity
! ints(5) = helicity
! ints(6) = delke_tot = Actual d(KE)/dt computed from KE(t) and KE(t+1)
! ints(7) = enstrophy (vorticity**2)
! ints(8) = ke dissapation from alpha model term 
!
! maxs(5) = max vorticity
!
! for convience, we store the time data in maxs(6:7)
! maxs(6) = ints_timeU
! maxs(7) = ints_timeDU
!

integer,parameter :: ntimers=12
real*8 :: tims(ntimers)=0
!  tims(1)    time for initialization
!  tims(2)    total runtime after initialization
!  tims(3)    time spent in time_control()
!  tims(4)    not used
!  tims(5)    time spent in RHS
!  tims(6)    transpose_to_z
!  tims(7)    transpose_from_z
!  tims(8)    transpose_to_x
!  tims(9)    transpose_from_x
!  tims(10)    transpose_to_y
!  tims(11)    transpose_from_y
!  tims(12)    compute_pdf
!
!







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! parallel decompositions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!integer :: ncpu_x=1,ncpu_y=1,ncpu_z=1
integer :: nx_2dy,ny_2dz,nz_2dx,ny_2dx
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
!
! #ifdef TRANSPOSE_X_SPLIT_Z 
!    ncpu_x divides (nz2-nz1+1)   nz_2dx=(nz2-nz1+1)/ncpu_x
!    ncpu_y divides (nx2-nx1+1)   nx_2dy=(nx2-nx1+1)/ncpu_y
!    ncpu_z divides (ny2-ny1+1)   ny_2dz=(ny2-ny1+1)/ncpu_z
!
! #ifdef TRANSPOSE_X_SPLIT_Y        (usefull if ncpu_z=1)
!    ncpu_x divides (ny2-ny1+1)   ny_2dx=(ny2-ny1+1)/ncpu_x
!    ncpu_y divides (nx2-nx1+1)   nx_2dy=(nx2-nx1+1)/ncpu_y
!    ncpu_z divides (ny2-ny1+1)   ny_2dz=(ny2-ny1+1)/ncpu_z
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




contains

subroutine params_init
character(len=80) message
integer :: fail=0

g_nx=nslabx*ncpu_x
g_ny=nslaby*ncpu_y
g_nz=nslabz*ncpu_z
g_nx2=g_nx+2
g_ny2=g_ny+2
g_nz2=g_nz+2
if (g_nz==1) g_nz2=1


! these values must divide with no remainder:
nz_2dx=nslabz/ncpu_x     ! TRANSPOSE_X_SPLIT_Z 
ny_2dx=nslaby/ncpu_x     ! TRANSPOSE_X_SPLIT_Y  (usefull if nslabz=1)
nx_2dy=nslabx/ncpu_y     ! TRANSPOSE_Y_SPLIT_X  (always used)
ny_2dz=nslaby/ncpu_z     ! TRANSPOSE_Z_SPLIT_Y  (always used)

#if (!defined TRANSPOSE_X_SPLIT_Z && !defined TRANSPOSE_X_SPLIT_Y)
   call abort("define TRANSPOSE_X_SPLIT_Y or TRANSPOSE_X_SPLIT_Z in transpose.h") 
#endif

#ifdef TRANSPOSE_X_SPLIT_Z
ny_2dx=-1  ! set to an invalid value
if (0/=mod(nslabz,ncpu_x)) then
   fail=1
   call print_message("ncpu_x does not divide nz");
endif
#endif

#ifdef TRANSPOSE_X_SPLIT_Y
nz_2dx=-1  ! set to an invalid value
if (0/=mod(nslaby,ncpu_x)) then
   fail=1
   call print_message("ncpu_x does not divide ny");
endif
#endif

if (0/=mod(nslabx,ncpu_y)) then
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
#ifdef TRANSPOSE_X_SPLIT_Z
if ((g_nx2)*nz_2dx*real(nslaby,r8kind) > nx*nz*real(ny,r8kind) )  then
   fail=1
   call print_message("insufficient storage.");
   write(message,'(a,3i6,a,f10.0)') "nx,ny,nz=",nx,ny,nz,"   nx*ny*nz=",nx*ny*real(nz,r8kind)
   call print_message(message)	
   write(message,'(a,f10.0)') "storage needed for 2D x-decomposition: ", &
     (g_nx2)*nz_2dx*real(nslaby,r8kind)
   call print_message(message)	
   call print_message("You might also try #define TRANSPOSE_X_SPLIT_Y in transpose.h")

endif
#endif


#ifdef TRANSPOSE_X_SPLIT_Y
if ((g_nx2)*ny_2dx*real(nslabz,r8kind) > nx*nz*real(ny,r8kind) )  then
   fail=1
   call print_message("insufficient storage.");
   write(message,'(a,3i6,a,f10.0)') "nx,ny,nz=",nx,ny,nz,"   nx*ny*nz=",nx*ny*real(nz,r8kind)
   call print_message(message)	
   write(message,'(a,f10.0)') "storage needed for 2D x-decomposition: ", &
     (g_nx2)*ny_2dx*real(nslabz,r8kind)
   call print_message(message)	
   call print_message("You might also try #define TRANSPOSE_X_SPLIT_Z in transpose.h")

endif
#endif



if ( (g_ny2)*nx_2dy*real(nslabz,r8kind) > nx*ny*real(nz,r8kind)) then
   fail=1
   call print_message("insufficient storage.");
   write(message,'(a,3i6,a,f10.0)') "nx,ny,nz=",nx,ny,nz,"   nx*ny*nz=",nx*ny*real(nz,r8kind)
   call print_message(message)	
   write(message,'(a,f10.0)') "storage needed for 2D y-decomposition: ", &
     (g_ny2)*nx_2dy*real(nslabz,r8kind)
   call print_message(message)	

endif

if ( (g_nz2)*ny_2dz*real(nslabx,r8kind) >  ny*nz*real(nx,r8kind)) then
   fail=1
   call print_message("insufficient storage.");
   write(message,'(a,3i6,a,f10.0)') "nx,ny,nz=",nx,ny,nz,"   nx*ny*nz=",nx*ny*real(nz,r8kind)
   call print_message(message)	
   write(message,'(a,f10.0)') "storage needed for 2D z-decomposition: ", &
     (g_nz2)*ny_2dz*real(nslabx,r8kind)
   call print_message(message)	
endif


!
! this is needed to gaurentee that the sine and cosine modes
! are on the same processor.  
!
if (mod(nslabx,2)/=0) then
   fail=1
   call print_message("nslabx must be even")
endif
if (mod(nslaby,2)/=0) then
   fail=1
   call print_message("nslaby must be even")
endif
if (mod(nslabz,2)/=0 .and. g_nz>1) then
   fail=1
   call print_message("nslabz must be even if g_nz>1")
endif


!
!  for the full spectral method, we will be working mostly in
!  the z-transform fourier space, 
!     pt(g_nz2,nslabx,ny_2dz)
!  so we also need that ny_2dz is even
if (mod(ny_2dz,2)/=0) then
   fail=1
   call print_message("ny_2dz is not even.  cant run the ns3dspectral model")
   call print_message("but other models should work fine")
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

allocate(g_xcord(g_nx+1))
allocate(g_ycord(g_ny+1))
allocate(g_zcord(g_nz+1))
allocate(g_imcord(g_nx))
allocate(g_jmcord(g_ny))
allocate(g_kmcord(g_nz))

!  Z-decomposition modes:   p(g_nz2,nslabx,ny_2dz)
allocate(z_imcord(nx))
allocate(z_jmcord(ny_2dz))
allocate(z_kmcord(g_nz))
allocate(z_imsign(nx))
allocate(z_jmsign(ny_2dz))
allocate(z_kmsign(g_nz))



pi=1
pi=4*atan(pi)
pi2=2*pi
pi2_squared=4*pi*pi


end subroutine
                      

end ! module mod_params








