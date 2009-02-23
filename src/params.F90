!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Copyright 2007.  Los Alamos National Security, LLC. This material was
!produced under U.S. Government contract DE-AC52-06NA25396 for Los
!Alamos National Laboratory (LANL), which is operated by Los Alamos
!National Security, LLC for the U.S. Department of Energy. The
!U.S. Government has rights to use, reproduce, and distribute this
!software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
!LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
!FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
!derivative works, such modified software should be clearly marked, so
!as not to confuse it with the version available from LANL.
!
!Additionally, this program is free software; you can redistribute it
!and/or modify it under the terms of the GNU General Public License as
!published by the Free Software Foundation; either version 2 of the
!License, or (at your option) any later version. Accordingly, this
!program is distributed in the hope that it will be useful, but WITHOUT
!ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
!for more details.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "macros.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  module containing all global parameters like mesh data
!
! To create a valid file:  
! pick ncpu_x,ncpu_y,ncpu_z	
! Then l,m,n, and set:
!
!
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
integer :: mu_hyper=1         !viscosity = (del**2)**mu_hyper
integer :: hyper_implicit=0   ! 1 = time-split, implicit hyper viscosity
real*8  :: mu_hyper_value=0   ! hyper viscosity value
real*8  :: k_Gt=0             ! JZ hyper viscosity scaling parameter
real*8  :: max_hyper          ! some diagnostics
integer :: mu_hypo =0         !large wave number dissipation = (del**-2)**mu_hyper
real*8  :: mu_hypo_value=0    !viscosity value 
real*8  :: alpha_value=0      !for the alpha model  
real*8  :: alpha_B=0          !for the alpha Bardina model  
integer :: infinite_alpha=0   !flag for infinite alpha case
real*8  :: smagorinsky=0      !for smagorinsky term
real*8  :: H0=0               ! used by shallow water model
integer,parameter :: r8kind=kind(mu)
integer :: equations=NS_UVW        ! NS_UVW    = NS / NS-alpha
                                   ! SHALLOW   = Shallow water / shallow water-alpha
                                   ! NS_PSIVOR = NS Streamfunction-Vorticity
integer :: numerical_method=FOURIER ! FOURIER
                                    ! FOURTH_ORDER
                                    ! others?

integer :: dealias=0                ! 0 = none
                                    ! 1 = 2/3 rule (exact)
                                    ! 2 = sqrt(2)*N/3 rule (approximate)

logical :: use_phaseshift=.false.   ! .true. = compute nonlinear term with phase shifting

!
! I/O options
!
logical :: r_spec=.false.           ! set to .true. to read spectral coefficieints 
logical :: w_spec=.false.           ! set to .true. to write spectral coefficieints 
                                    ! instead of grid point values 

integer :: spec_max=-1              ! if w_spec=.true. (spectral output)
                                    ! this is the max number of spectral modes to 
                                    ! write to output file
                                    ! -1 = use dealiased maximum
                                    ! 
                                    ! It is also used by the -cout trunc
                                    ! option in convert.F90.  Here it is a
                                    ! parameter to the truncation filter
                                    !  
                                 
logical :: r_compressed=.false.       ! NS_UVW Kerr compression
logical :: w_compressed=.false.       ! 
                                 


logical :: udm_input=.false.        ! use UDM for input/restart files
logical :: udm_output=.false.       ! use UDM for output files
logical :: do_mpi_io=.false.     
logical :: byteswap_input=.false.   ! byteswap all binary input files
integer :: output_size=8            ! output real*8
integer :: input_size=8             ! input real*8
logical :: compute_passive_on_restart = .false.

logical :: data_x_pencils=.false.   ! Q currently stored in x-pencil decomp
                                    ! not our regular (nx,ny,nz) decomp
logical :: use_vorticity3=.true.    ! compute vorticity during u iFFT
!
! options used by various utility programs:
!
                                    ! used by the 'convert' program:
integer :: convert_opt=-1           ! 0   output uvw  
                                    !     usefull to convert to/from UDM
                                    ! 1   output vor1,2,3 
                                    ! 2   output vor magnitude 
integer :: user_specified_isodir=-1   ! set number of directions and 
integer :: user_specified_isodel=-1   ! number of seperations in isoave


real*8  :: g_u2xave=0               ! <ux,ux> updated after each time step  

! xscale and yscale are used by psi-vor model
! they are used to modify the computational domain from the standard
! [0..1] x [0..1] to [0..xscale] x [0..yscale]
real*8 :: xscale=1
real*8 :: yscale=1
real*8 :: zscale=1


! Lz is used by the DNS model.  The code is solving
! the equations in the physical domain [0..1]x[0..1]x[0..Lz], but
! the computational domain remains [0..1]x[0..1]x[0..1]
! The equations are mapped to the computational domain
! (all d/dz derivatives become d/dz (1/Lz) See rotation.tex
real*8 :: Lz=1


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
character(len=80) :: inputfile=''
real*8  :: pi,pi2,pi2_squared

real*8  :: grav=0  
real*8  :: fcor=0                ! rotation in z-axis 
real*8  :: bous=0                ! bousenesque paramter
integer :: diag_struct=0         ! compute structure funtions on the fly
integer :: diag_pdfs=0           ! 0  disabled
                                 ! 1  compute pdfs on the fly
                                 ! -1  compute pdfs if simulation time >= 1.0 


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
real*8  :: init_cond_param1 = 0  ! Initial condition parameters
real*8  :: init_cond_param2 = 0  ! 0 = default  
                                 ! value interpreted by initial condition
                                 ! subroutine. see source. 
                                 !
                     
integer :: forcing_type   ! 0 = none
                          ! 1 = relax back to E(1)=1, E(2)=2**(-5/3)
                          !     can only be used by the z-decomp model!
integer :: forcing_peak_waveno = 0
real*8  :: ffval = 1    !  normalize some forcings so <f,f> = ffval
real*8  :: fparam1 = 0  !  another forcing parameter. meaning depends
                        ! on forcing_type, may be ignored 


integer :: header_user=1    ! I/O header type from input file


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! mesh dimensions on a single processor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "params.h"
!integer,parameter :: n_var=3                  ! number of prognostic variables
!integer,parameter :: nxd=18,nyd=18,nzd=18         ! dimension of grid & data
!integer,parameter :: nx1=1,nx2=16              ! upper and lower bounds of non-ghost data
!integer,parameter :: ny1=1,ny2=16             ! upper and lower bounds of non-ghost data
!integer,parameter :: nz1=1,nz2=16             ! upper and lower bounds of non-ghost data

! dimensions of interior, non-boundary points:
!   do k=intz1,intz2  
!   do j=inty1,inty2
!   do i=intx1,intx2
!
! to loop over all interior points and boundary points:
!   do k=bz1,bz2  
!   do j=by1,by2
!   do i=bx1,bx2
!
! Code that deals with boundary conditions (other than periodic or reflections)
! should use these bounds and NOT nx1,nx2...
!
! In most cases,
!   INTERIOR boundaries:   intx1=bx1=nx1    intx2=bx2=nx2
!   PERIODIC boundaries:   intx1=bx1=nx1    intx2=bx2=nx2
!   REFLECT boundaries:    intx1=bx1=nx1    intx2=bx2=nx2
!   Real boundaries at x2: intx1=bx1=nx1    intx2+1=bx2=nx2
!   Real boundaries at x1: intx1-1=bx1=nx1  intx2=bx2=nx2
!
! But for the psi-vor model, non-periodic, it needs a sine transform.
! Normally, on a 400x400 grid would need a len=399 transform.  
! So for this model we sometimes set, along the bdy_x2 boundary:
!  intx2=nx2   bx2=nx2+1
! which places the boundary in the ghost cell region.  
! 
! This is okay, as long as the boundary code is not relying 
! on two rows of ghost cells beyond the boundary.  Also, the
! transpose() routines will not pick up the ghost cell data.
!
! For the sine transform, this is no problem since that end point
! data is required to be zero.
! It is a problem for output, since that boundary data in the
! ghost cell region will not be included.
! Parallel transpose operations always work just on the nx1,nx2 data.
! Normally this is not a problem since code with boundary conditions
! will not be using FFT's, except for the sine transform.  In that case,
! if bx2<>nx2, the boundary data will not make it into the transposed
! array.  But this is fine for the sine transform since it requires
! that data to be zero.  The output routines also rely on the 
! transpose operatorors, and they will not output the data
!
!
!

integer :: intx1,intx2
integer :: inty1,inty2
integer :: intz1,intz2
integer :: bx1,bx2
integer :: by1,by2
integer :: bz1,bz2
integer :: offset_bdy=0   ! flag indicating boundary data should be placed   
                          ! in ghost cell of bdy_x2 (or bdy_y2)



! compiler has three choices, somewhat controllabe with options:
!  static  (if all arrays are static, memory use is doubled)
!  stack   Ideal, but some systems have limited stack sizes
!  heap    malloc()'d and free()'d between subroutine calls.  
!          Performance penalty? And causes problems with Elan
!          on ASCI Q when using more than 500M per cpu.  
! 
integer,parameter :: nx=nxd,ny=nyd,nz=nzd      ! dimension of grid & data


integer :: ndim       ! 2 or 3 dimensions.  ndim = (g_nz==1 ? 2 : 3)
integer :: npassive = 0   ! number of passive scalars
integer :: np1=1,np2=0    ! passive scalars Q(:,:,:,np1:np2)
real*8  :: schmidt(n_var)       ! schmidt number for passive scalars
integer :: passive_type(n_var)  ! type of passive scalar:
                                !     0 = Gaussian i.c., regular equation
                                !     1 = KE-correlated i.c., regular equation
                                !     2 = 0 i.c., Kurien's gage equation

integer :: input_npassive = 0                   ! npassive in input file
real*8,allocatable :: input_schmidt(:)          ! data from input file
integer,allocatable :: input_passive_type(:)    ! data from input file

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

! fft modes, local 3D decompostion, sin/cos basis
integer :: imcord(nxd),jmcord(nyd),kmcord(nzd)  ! fft modes local
integer :: imsign(nxd),jmsign(nyd),kmsign(nzd)  ! fft modes local
! fft modes, local 3D decomposition, exponential basis
! these numbers are always used in pairs, so if nz2 is odd, max
! dimension needs to be nz2+1.  Since nz2<=nzd, just add 1 to be save:
integer :: imcord_exp(nxd+1),jmcord_exp(nyd+1),kmcord_exp(nzd+1) 

! fft modes, local z-decompostion
integer,allocatable :: z_imcord(:),Z_jmcord(:),z_kmcord(:)  ! fft modes local
integer,allocatable :: z_imsign(:),z_jmsign(:),z_kmsign(:)  ! fft modes local
! fft modes, local x-decompostion
integer,allocatable :: x_imcord(:),X_jmcord(:),x_kmcord(:)  ! fft modes local
integer,allocatable :: x_imsign(:),x_jmsign(:),x_kmsign(:)  ! fft modes local

! fft modes, p3dfft
integer, allocatable :: p3_1mcord(:),p3_2mcord(:),p3_3mcord(:)

! sine modes, local 3D decomp:
integer,allocatable :: y_imsine(:)  ! ny_2dz
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
integer           :: g_nmax
integer           :: g_nmin

! spherical dealias kmax
integer           :: dealias_sphere_kmax2 
integer           :: dealias_sphere_kmax2_1 
real*8            :: dealias_sphere_kmax 

integer           :: dealias_23sphere_kmax2 
integer           :: dealias_23sphere_kmax2_1 
real*8            :: dealias_23sphere_kmax 

integer           :: dealias_23_kmax2 
integer           :: dealias_23_kmax2_1 
real*8            :: dealias_23_kmax 
! actuall largest wave number in the code (including the 2pi factor
! because our domain is of size 1)
real*8            :: xkmax

integer :: o_nx,o_ny,o_nz    ! dimensions of plotting output data
                             ! For periodic FFT case, o_nx=g_nx+1 because we do not
                             ! store the point at x=1, but we output it.
                             ! for 4th order case, we may store this point.  
                             
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
logical :: E_kmax_cfl = .false.    ! cfl_adv = delt*sqrt(ke)*kmax

real*8 :: output_dt = 0           ! output prognostic variables
integer :: output_vorticity = 0   ! also output vorticity
integer :: ncustom =0
real*8, allocatable :: custom(:)
real*8 :: diag_dt = 0       ! diagnostics
real*8 :: model_dt = 0       ! unused
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
integer,parameter :: nints=16
real*8 :: ints(nints),maxs(nints)

!
! For alpha>0 and alpha==0.      H(f')=f   H(f)=Helmholz(f)
!
! KE = ints(6)
! KE dissapation:  ints(10) + ints(3) + ints(8)
!      ints(10) = -mu*< u_x,u_x >
!      ints(3) = < u,f  >           
!      ints(8) = < u,div(tau)' >
!
! E-alpha = ints(6) + .5*alpha**2 ints(2)
! E-alpha dissapation =   ints(10) - mu*alpha**2 * ints(1) + 
!                          ints(3) - alpha**2*ints(9) 
!      ints(1)= < u_xx,u_xx>
!      ints(9) = < uxx,f >
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
! ints(3) = < u,f > 
! ints(4) = integral of z-component of vorticity
! ints(5) = helicity  (or: mu*<w,w_xx> enstrophy diffusion)
! ints(6) = .5 < u,u >
! ints(7) = enstrophy (vorticity**2)
! ints(8) = < u,div(tau)' >   (alpha model only)
! ints(9)  = < uxx,f > 
! ints(10) = mu*< u, del u >   diffusion (or mu*<u,-del**4 u>)
! ints(11) = helicity dissipation
! ints(12) = enstrophy dissipation k^2 |w|^2
! ints(13) = enstrophy dissipation k^4 |w|^2
! ints(14) = enstrophy dissipation k^6 |w|^2
! ints(15) = potential energy (when running Boussinesque)
! ints(16) = potential energy dissipation rate
! maxs(5) = max vorticity
!
!
! for convience, we store the time data in maxs(6:7)
! maxs(6) = time (at end of time step)  T
! maxs(7) = time at begining of time step T-1
! maxs(8) = not used
! maxs(9) = error conditions.  0 = no error
!                              1 = kill signal detected by handler
! maxs(10) = max of 1st passive scalar (if npassive>0)
! maxs(11) = max of negative of 1st passive scalar (used to copute min)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for the shallow water model, we modify the above slightly:
! ints(1) = < u_xx,-del**4 u >   
! ints(2) = < h u_x,u_x >   
! ints(5) = .5 < h u , u >                 KE
! ints(6) = .5 < h u , u >  + g H^2        KE+PE
! ints(8) = < h u,div(tau)' >              NOT debugged yet
!                              
! ints(10) = mu*< u, -del**4 u >   hyper diffusion term + smag. diffusion
!
! E = ints(6) = KE+PE
! dE/dt = ints(10)   + ints(8)
!
! E_alpha = ints(6) + .5*alpha**2 ints(2)
! E_alpha dissapation: ints(10) - mu*alpha**2 * ints(1)
!
!
!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer,parameter :: ntimers=24
real*8 :: tims(ntimers)=0
integer :: ncalls(ntimers)=0
!  tims(1)    time for initialization
!  tims(2)    total runtime after initialization
!  tims(3)    time spent in time_control() after initialization
!  tims(4)    time spend in time_control() during initializaiton
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
!  tims(15)    compute_psi
!  tims(16)    tracer_advance
!  tims(17)    ellipitical contours
!  tims(18)    FFT
!  tims(19)    iFFT
!  tims(20)    transpose_to_z_from_refyxz
!  tims(21)    transpose_from_z_to_refyxz
!  tims(22)    transpose_to_x_from_refyxz
!  tims(23)    transpose_from_x_to_refyxz
!  tims(24)    timing for single call to zx_fft3d_trashinput()
!
!







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! parallel decompositions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!integer :: ncpu_x=1,ncpu_y=1,ncpu_z=1
integer :: ny_2dx=nslaby/ncpu_x       ! y-dim in x-pencil-decomposition
integer :: nx_2dy=nslabx/ncpu_y       ! x-dim in y-pencil-decomposition
integer :: ny_2dz=nslaby/ncpu_z       ! y-dim in z-pencil-decomposition 
integer :: nx_2dz=nslabx              ! x-dim in z-pencil-decomposition 

integer :: p3_n1, p3_n2, p3_n3        ! P3DFFT loop bounds

integer :: io_pe
integer :: my_pe,mpicoords(3)         ! processor IDs
integer :: my_world_pe,mpidims(3)     ! used by init_mpi ONLY
logical :: p3dfft_cartmap = .false.   ! try to mimic p3dfft cartesian map
  
integer :: initial_live_procs      ! total number of cpus
integer :: ncpus                   ! number of cpus in cartesian communicator
integer :: comm_3d                 ! the MPI cartesian communcator
integer :: comm_1                  ! communicator with just io_pe



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
!    ncpu_x divides (ny2-ny1+1)   ny_2dx=(ny2-ny1+1)/ncpu_x
!    ncpu_y divides (nx2-nx1+1)   nx_2dy=(nx2-nx1+1)/ncpu_y
!    ncpu_z divides (ny2-ny1+1)   ny_2dz=(ny2-ny1+1)/ncpu_z
! OR: 
!    ncpu_z divides (nx2-nx1+1)   nx_2dz=(nx2-nx1+1)/ncpu_z
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




contains




subroutine params_init

character(len=80) message
integer :: fail=0
integer :: dims(2)
integer :: i,j,k,iw,jw,kw

g_nx2=g_nx+2
g_ny2=g_ny+2
g_nz2=g_nz+2
ndim=3

g_nmax=max(g_nx,g_ny,g_nz)
g_nmin=min(g_nx,g_ny,g_nz)

if (g_nz==1) then
   g_nz2=1
   ndim=2
   g_nmax=max(g_nx,g_ny)
   g_nmin=min(g_nx,g_ny)
endif

! spherical dealiasing
! pick the sphere that fits inside l^2 + m^2 + n^2 < 2N^2/9
! but lets truncate so that we only have complete shells:
!     k-.5 <= sqrt(l^2 + m^2 + n^2))  < k+.5 
!   k^2-k+.25  <=  l^2 + m^2 + n^2   < k^2 + k + .25
! (since k and LHS are integers, this is the same as:)
!   k^2-k   <  l^2 + m^2 + n^2   <= k^2 + k 
!
! which implies remove waves:
!             l^2 + m^2 + n^2   > k^2 + k
!
! find integer (k+.5)^2   <   2N^2/9
!                      k  <  sqrt(2) N / 3 -.5 = (sqrt(2) N -1.5) / 3
!
!
! Example: N=24  phase shift dealiasing allows up to k^2 < 128  k<11.3
! mode=(11,3,0) k=11.4    (part of the k=11 shell)  REMOVED
! mode=(11,2,0) k=11.18   (part of the k=11 shell)  KEPT
!
! Thus in this case, the k=11 shell is incomplete.  Thus we dealias
! a little extra and remove *all* modes outside the k=10 shell
! kmax = 10
! kmax2 = 110    (instead of 128)
!
#define TRUNC_SHELL
#ifdef TRUNC_SHELL
! these formulas truncate at a complete k spectrum shell boundary k+.5
dealias_sphere_kmax = floor(  (sqrt(2d0)*g_nmin - 1.5d0)/3 )
dealias_sphere_kmax2   = (dealias_sphere_kmax+.5)**2  ! rounds to k^2 + k
dealias_sphere_kmax2_1 = (dealias_sphere_kmax-.5)**2  ! rounds to k^2 - k
#else
! this formula keeps all modes that are fully dealiased
! so there will be a prartial shell in the spectrum
! we remove  k^2  >= 2N^2/9
! if this is NOT an integer, that is the same as:
!           k^2  >= 1 + floor(2N^2/9)
! if this IS an integer, that is the same as  
!           k^2  >= nint(2N^2/9)
if (mod(g_nmin,3)==0) then
   dealias_sphere_kmax2 = 2*(g_nmin/3)**2
else
   dealias_sphere_kmax2 = 1 + floor(  (2*g_nmin*g_nmin)/9d0 )
endif
dealias_sphere_kmax = sqrt(real(dealias_sphere_kmax2))
dealias_sphere_kmax2_1 = floor( (dealias_sphere_kmax-1)**2 )
#endif


! this dealiasing is used mostly for debugging:
! pick the sphere that fits inside the 2/3 cube: kmax = N/3.  
! but lets truncate so that we only have complete shells:
!     k-.5 <= sqrt(l^2 + m^2 + n^2))  < k+.5 
!   k^2-k+.25  <=  l^2 + m^2 + n^2   < k^2 + k + .25
! (since k and LHS are integers, this is the same as:)
!   k^2-k   <  l^2 + m^2 + n^2   <= k^2 + k 
!
! which implies remove waves:
!             l^2 + m^2 + n^2   > k^2 + k
!
! find integer (k+.5)^2   <=  (N/3)^2
!                      k  <= (N-1.5)/3 
! so take k = floor( (N-1.5)/3 )   
dealias_23sphere_kmax = floor ( (g_nmin-1.5d0)/3 )
dealias_23sphere_kmax2   = (dealias_23sphere_kmax+.5)**2  ! rounds to k^2 + k
dealias_23sphere_kmax2_1 = (dealias_23sphere_kmax-.5)**2  ! rounds to k^2 - k


write(message,'(a,3i6)') "x-pencil decomp array size: ",g_nx2,ny_2dx,nslabz
call print_message(message)
write(message,'(a,3i6)') "y-pencil decomp array size: ",nx_2dy,g_ny2,nslabz
call print_message(message)
write(message,'(a,3i6)') "z-pencil decomp array size: ",nx_2dz,ny_2dz,g_nz2  
call print_message(message)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  z-pencil decomposition tweaks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z-pencil decomposition.  can be anything of the form
!   nx_2dz=nslabx/z1         mod(nslabx,nx_2dz)==0
!   ny_2dz=nslaby/z2         mod(nslaby,ny_2dz)==0
! with z1*z2 = ncpu_z
!
! default was to chop up y direction
! usually nslaby will be larger than nslabx, except in P3DFFT case.  
! so if nslabx > nslaby, then chop up the x direction:
if (nslabx > nslaby ) then
   ny_2dz=nslaby              ! y-dim in z-pencil-decomposition 
   nx_2dz=nslabx/ncpu_z       ! x-dim in z-pencil-decomposition 
   write(message,'(a,3i6)') "Modifying z-pencil decomp.  array size: ",nx_2dz,ny_2dz,g_nz2  
   call print_message(message)
endif

!
!
! if z-pencil decomposition has y dimesnion of only 1, lets
! double it.   (only transpose_to/from_z() routines have been 
! modified to handle this)
!
! This allows us to run a 1x1xN decompostion for a N^3 grid
! 
! check: 
!      ncpu_z = (nslabx/nx_2dz) * (nslaby/ny_2dz)
! which means that:  nslabx, nslaby must be even
!           
! Does this code work in the non-perfect load balanced case?
! Probably not - so lets check for that too.
!
if (ny_2dz == 1 .and. 2*ny_2dz <= nslaby .and. nx_2dz /= 1  ) then
   ! double the size in y, half the size in x:
   ny_2dz=2*ny_2dz
   nx_2dz=nx_2dz/2

   write(message,'(a,3i6)') "Modifying z-pencil decomp.  array size: ",nx_2dz,ny_2dz,g_nz2  
   call print_message(message)

   if (0/=mod(nslabx,nx_2dz)) then
      fail=1
      call print_message("Error: nslabx/nx_2dz problem");
   endif
   if (0/=mod(nslaby,ny_2dz)) then
      fail=1
      call print_message("Error: nslaby/ny_2dz problem");
   endif
   if (ncpu_z /= (nslabx/nx_2dz)*(nslaby/ny_2dz)) then
      fail=1
      call print_message("Strange error in tranpose_to_z setup - check!")
   endif

   ! also make sure we are perfectly load balanced:
   ! we allow this, but it is not tested with this tweaked transpose_to/from_z
   if (0/=mod(nslaby,ncpu_x)) then
      fail=1
      call print_message("Error: ncpu_x does not divide nslaby");
   endif
   if (0/=mod(nslabx,ncpu_y)) then
      fail=1
      call print_message("Error: ncpu_y does not divide nslabx");
   endif
   if (0/=mod(nslaby,ncpu_z)) then
      fail=1
      call print_message("Error:  ncpu_z does not divide nslaby");
   endif
endif






! for perfect load balancing in x,y or z pencil decomposition space,
! these values must divide with no remainder.
!
! violating these conditions should now work.  tested by testvx.sh
! but only for 2D problems?
!
if (ncpu_x*ny_2dx<nslaby) then
   call print_message("*WARNING*:  transpose_to_x not perfectly load balanced")
   ny_2dx=ny_2dx+1
   ny_2dx=ny_2dx + mod(ny_2dx,2)
endif
if (ncpu_y*nx_2dy<nslabx) then
   nx_2dy=nx_2dy+1
   nx_2dy=nx_2dy + mod(nx_2dy,2)
   call print_message("*WARNING*:  transpose_to_y not perfectly load balanced")
endif
if ( mod(nslabx,nx_2dz)/=0 ) then
   nx_2dz=nx_2dz+1
   nx_2dz=nx_2dz + mod(nx_2dz,2)
   call print_message("*WARNING*:  transpose_to_z not perfectly load balanced")
endif
if ( mod(nslaby,ny_2dz)/=0 ) then
   ny_2dz=ny_2dz+1
   ny_2dz=ny_2dz + mod(ny_2dz,2)
   call print_message("*WARNING*:  transpose_to_z not perfectly load balanced")
endif









!
! memory contraint: pencil decomposition should fit in a 3D decomposition array
! check if: (g_nz2)*nslabx*ny_2d <= nx*ny*nz)
!
if ( (g_nz2)*ny_2dz*real(nx_2dz,r8kind) >  ny*nz*real(nx,r8kind)) then
   fail=1
   call print_message("insufficient storage.");
   write(message,'(a,3i6,a,f10.0)') "nx,ny,nz=",nx,ny,nz,"   nx*ny*nz=",nx*ny*real(nz,r8kind)
   call print_message(message)	
   write(message,'(a,f10.0)') "storage needed for 2D z-decomposition: ", &
     (g_nz2)*ny_2dz*real(nslabx,r8kind)
   call print_message(message)	
endif
if ( (g_ny2)*nx_2dy*real(nslabz,r8kind) > nx*ny*real(nz,r8kind)) then
   fail=1
   call print_message("insufficient storage.");
   write(message,'(a,3i6,a,f10.0)') "nx,ny,nz=",nx,ny,nz,"   nx*ny*nz=",nx*ny*real(nz,r8kind)
   call print_message(message)	
   write(message,'(a,f10.0)') "storage needed for 2D y-decomposition: ", &
     (g_ny2)*nx_2dy*real(nslabz,r8kind)
   call print_message(message)	

endif
if ((g_nx2)*ny_2dx*real(nslabz,r8kind) > nx*nz*real(ny,r8kind) )  then
   fail=1
   call print_message("insufficient storage.");
   write(message,'(a,3i6,a,f10.0)') "nx,ny,nz=",nx,ny,nz,"   nx*ny*nz=",nx*ny*real(nz,r8kind)
   call print_message(message)	
   write(message,'(a,f10.0)') "storage needed for 2D x-decomposition: ", &
     (g_nx2)*ny_2dx*real(nslabz,r8kind)
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




if (fail/=0) call abortdns("params.F90 dimension settings failure")

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

!  X-decomposition modes:   p(g_nz2,nslabx,ny_2dz)
allocate(x_imcord(g_nx))
allocate(x_jmcord(ny_2dx))
allocate(x_kmcord(nz))
allocate(x_imsign(g_nx))
allocate(x_jmsign(ny_2dx))
allocate(x_kmsign(nz))

allocate(y_imsine(nx_2dy))



pi=1
pi=4*atan(pi)
pi2=2*pi
pi2_squared=4*pi*pi



end subroutine







logical function dealias_remove(im,jm,km)
implicit none
integer :: im,jm,km

ASSERT("dealias_remove: im >= 0 failed",im>=0)
ASSERT("dealias_remove: jm >= 0 failed",jm>=0)
ASSERT("dealias_remove: km >= 0 failed",km>=0)
if (dealias==1) then
   dealias_remove = ( (3*km>=g_nz)  .or.  (3*jm>=g_ny)  .or. (3*im>=g_nx) )
else if (dealias==2) then
#ifdef TRUNC_SHELL
   dealias_remove = ( ( im**2 + jm**2 + km**2 ) > dealias_sphere_kmax2 )
#else
   dealias_remove = ( ( im**2 + jm**2 + km**2 ) >= dealias_sphere_kmax2 )
#endif
else if (dealias==3) then
   dealias_remove = ( ( im**2 + jm**2 + km**2 ) > dealias_23sphere_kmax2 )
else
   dealias_remove = .false.
endif
   
end function

                      

end ! module mod_params








