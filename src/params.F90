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
real*8  :: mu=0               !viscosity
real*8  :: mu_hyper=0         !set to 1 to enable hyper viscosity
real*8  :: alpha_value=0      !for the alpha model  
real*8  :: H0=0               ! for shallow water model
integer,parameter :: r8kind=kind(mu)
integer :: enable_lsf_timelimit = 0
integer :: equations=NS_UVW        ! NS_UVW    = NS / NS-alpha
                                   ! SHALLOW   = Shallow water / shallow water-alpha
                                   ! NS_PSIVOR = NS Streamfunction-Vorticity
integer :: numerical_method=FOURIER ! FOURIER
                                    ! FORTH_ORDER
                                    ! others?
logical :: dealias=.false.       

! parameter used by psi-vor model:
! before using xscale,yscale,zscale, we need to update all FFT
! derivatives and replace pi by pi/scale
real*8 :: xscale=1
real*8 :: yscale=1
real*8 :: zscale=1
real*8 :: biotsavart_cutoff
real*8 :: ubar=0


! local boundary conditions
! usually these will be INTERNAL, meaning the boundary is just
! an internal processor boundary
integer :: bdy_x1=PERIODIC
integer :: bdy_x2=PERIODIC
integer :: bdy_y1=PERIODIC
integer :: bdy_y2=PERIODIC
integer :: bdy_z1=PERIODIC
integer :: bdy_z2=PERIODIC

! global boundary conditions.  boundary conditions at x=0, 1, y=0,1
integer :: g_bdy_x1=PERIODIC
integer :: g_bdy_x2=PERIODIC
integer :: g_bdy_y1=PERIODIC
integer :: g_bdy_y2=PERIODIC
integer :: g_bdy_z1=PERIODIC
integer :: g_bdy_z2=PERIODIC


character(len=80) :: runname
character(len=80) :: rundir
real*8  :: pi,pi2,pi2_squared

real*8  :: grav=0  
real*8  :: fcor=0  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initial condition and forcing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: restart=0
integer :: init_cond      ! 0 = KH vortex blob
                          ! 1 = KH analytic
                          ! 2 = random, isotropic, .5 < wave no. < 2.5
                          !     E(1)=1, E(2)=2**(-5/3)
integer :: init_cond_subtype   !  0 = default
                               !  other vaules: set other parameters
                               !  in the initial condition
                               !
                               ! for KH anaylitic: 
                               !     0 = thin shear layer
                               !     1 = thick shear layer (E and Liu)
                               !
                               ! for iso12 with forcing:
                               !     0  tau=5
                               !     1  tau=20
                               !
                     
integer :: forcing_type   ! 0 = none
                          ! 1 = relax back to E(1)=1, E(2)=2**(-5/3)
                          !     can only be used by the z-decomp model!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! mesh dimensions on a single processor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "params.h"
!integer,parameter :: n_var=3                  ! number of prognostic variables
!integer,parameter :: nxd=18,nyd=18,nzd=18         ! dimension of grid & data
!integer,parameter :: nx1=1,nx2=16              ! upper and lower bounds of non-ghost data
!integer,parameter :: ny1=1,ny2=16             ! upper and lower bounds of non-ghost data
!integer,parameter :: nz1=1,nz2=16             ! upper and lower bounds of non-ghost data

! dimensions of interior, non-boundary points.  
! for problems with b.c. (other than periodic or reflections)
! use these bounds if you want to loop over only the non-boundary points.
integer :: intx1,intx2
integer :: inty1,inty2
integer :: intz1,intz2



! compiler has three choices, somewhat controllabe with options:
!  static  (if all arrays are static, memory use is doubled)
!  stack   Ideal, but some systems have limited stack sizes
!  heap    malloc()'d and free()'d between subroutine calls.  
!          Performance penalty? And causes problems with Elan
!          on ASCI Q when using more than 500M per cpu.  
! 
integer,parameter :: nx=nxd,ny=nyd,nz=nzd      ! dimension of grid & data


integer :: ndim       ! 2 or 3 dimensions.  ndim = (g_nz==1 ? 2 : 3)

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

! sine modes, local 3D decomp:
integer :: imsine(nxd)
integer :: jmsine(nyd)
integer :: kmsine(nzd)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! global dimensions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
integer,parameter :: g_nx=nslabx*ncpu_x
integer,parameter :: g_ny=nslaby*ncpu_y
integer,parameter :: g_nz=nslabz*ncpu_z


integer :: o_nx,o_ny,o_nz    ! dimensions of plotting output data
                             ! For periodic FFT case, o_nx=g_nx+1 because we do not
                             ! store the point at x=1.  for 4th order case, we 
                             ! may store this point.  
                             
integer :: g_nx2,g_ny2,g_nz2 ! dimension used by fft
                             ! (because FFT99 requires 2 extra points)
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! time stepping
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8  :: delt
real*8  :: cfl_adv = 1.2
real*8  :: cfl_vis = 1.2
real*8  :: delt_min = 0
real*8  :: delt_max = 1
real*8  :: time_final = 1 
real*8  :: time_initial = 0


real*8 :: output_dt = 0    ! netcdf output for plotting
integer :: ncustom =0
real*8, allocatable :: custom(:)
real*8 :: diag_dt = 0       ! diagnostics
real*8 :: screen_dt = 0     ! screen output
real*8 :: restart_dt = 0    ! restart 

#if 0
! set to > 0 to cause time stepping loop to exit.  if more than one cpu sets
! error_code > 0, the largest value will be reported.  
!
! error_code = 1   u > 1000
! error_code = 2   
! error_code = 3   
!
integer :: error_code =0
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! scalar quantities of current state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer,parameter :: nints=10
real*8 :: ints(nints),maxs(nints)

!
! For NON-ALPHA-MODEL
! KE = ints(6)
! KE dissapation:  ints(10) + ints(3) + ints(8)
!      ints(3) = < u,f >     
!      ints(8) = 0
!      ints(10) = < u, del u >  (or <u, -del**4 u)    
!
! for ALPHA-MODEL
!
! KE dissapation:  mu*ints(10) + ints(3) + ints(8)
!      ints(10) = < u_x,u_x >
!      ints(3) = < u,f' >           Helmholz(f')=f
!      ints(8) = < u,div(tau)' >
!
! E-alpha = ints(6) + .5*alpha**2 ints(2)
! E-alpha dissapation =   ints(10) + ints(9) - mu*alpha**2ints(1)
!      ints(9) = < u,f >
!      ints(1)= < u_xx,u_xx>
!  
! Note: we are going to change the forcing so that it always appears
! on the RHS as just f.  Then:
!    KE dissapation term:              <u,f>
!    E_alpha term:          <u,H f > = <u,f> + alpha**2 <uxx,f>
!   so make ints(3) = <u,f>
!           ints(9) = <uxx,f>
!
!
! These quantities are computed in the timestep
! routine, at the current time T:
!
! maxs(1,2,3) = max U,V,W
! maxs(4) = max (U/delx + V/dely + W/delz)  used for CFL
!
! These quantities computed as the RHS is computed
! and thus the data is at the prevous timestep  T-1:
! 
! ints(1) = < u_xx,u_xx >   (or < u_xx, -del**4 u>)
! ints(2) = < u_x,u_x >   
! ints(3) = < u,F >  Where F=forcing term which appears on RHS.
!                    non-alpha-mode: F=f.  Alpha model F=f'
! ints(4) = integral of z-component of vorticity
! ints(5) = helicity  (or: mu*<w,w_xx> enstrophy diffusion)
! ints(6) = .5 < u,u >
! ints(7) = enstrophy (vorticity**2)
! ints(8) = < u,div(tau)' >   (alpha model only)
! ints(9)  = < u,f >  (alpha model only)
! ints(10) = mu*< u, del u >   diffusion (or mu*<u,-del**4 u>)
! maxs(5) = max vorticity
!
!
! for convience, we store the time data in maxs(6:7)
! maxs(6) = time (at end of time step)  T
! maxs(7) = time at begining of time step T-1
! maxs(8) = time (in min) remaining from LSF job (or -1).
! maxs(9) = 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for the shallow water model, we modify the above slightly:
! ints(1) = < h u_xx,-del**4 u >   
! ints(2) = < h u_x,u_x >   
! ints(5) = .5 < h u , u >                 KE
! ints(6) = .5 < h u , u >  + g H^2        KE+PE
! ints(8) = < uH,div(tau)' >   (alpha model only. expensive to compute)  
!                              
! ints(10) = mu*< H u, -del**4 u >   hyper diffusion term
!
! E = ints(6)
! dE/dt = ints(10)   + ints(8)
!
! E_alpha = ints(6) + .5*alpha**2 ints(2)
! E_alpha dissapation: ints(10) - mu*alpha**2 * ints(1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer,parameter :: ntimers=14
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
!  tims(13)    ghost_update
!  tims(14)    biot savart
!
!







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! parallel decompositions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!integer :: ncpu_x=1,ncpu_y=1,ncpu_z=1
integer,parameter :: nz_2dx=nslabz/ncpu_x     ! TRANSPOSE_X_SPLIT_Z 
integer,parameter :: ny_2dx=nslaby/ncpu_x     ! TRANSPOSE_X_SPLIT_Y  (usefull if nslabz=1)
integer,parameter :: nx_2dy=nslabx/ncpu_y     ! TRANSPOSE_Y_SPLIT_X  (always used)
integer,parameter :: ny_2dz=nslaby/ncpu_z     ! TRANSPOSE_Z_SPLIT_Y  (always used)

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

g_nx2=g_nx+2
g_ny2=g_ny+2
g_nz2=g_nz+2
ndim=3

if (g_nz==1) then
   g_nz2=1
   ndim=2
endif



! these values must divide with no remainder:
!nz_2dx=nslabz/ncpu_x     ! TRANSPOSE_X_SPLIT_Z 
!ny_2dx=nslaby/ncpu_x     ! TRANSPOSE_X_SPLIT_Y  (usefull if nslabz=1)
!nx_2dy=nslabx/ncpu_y     ! TRANSPOSE_Y_SPLIT_X  (always used)
!ny_2dz=nslaby/ncpu_z     ! TRANSPOSE_Z_SPLIT_Y  (always used)

#if (!defined TRANSPOSE_X_SPLIT_Z && !defined TRANSPOSE_X_SPLIT_Y)
   call abort("define TRANSPOSE_X_SPLIT_Y or TRANSPOSE_X_SPLIT_Z in transpose.h") 
#endif

#ifdef TRANSPOSE_X_SPLIT_Z
if (0/=mod(nslabz,ncpu_x)) then
   fail=1
   call print_message("ncpu_x does not divide nz");
endif
#endif

#ifdef TRANSPOSE_X_SPLIT_Y
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








