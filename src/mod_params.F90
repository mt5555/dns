!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  module containing all global parameters
!  like mesh data
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_params

! mesh dimensions on a single processor
integer,parameter :: n_var=3                  ! number of prognostic variables
integer,parameter :: nx=18,ny=18,nz=18         ! dimension of grid & data
integer,parameter :: nx1=1,nx2=16              ! upper and lower bounds 
integer,parameter :: ny1=1,ny2=16              ! upper and lower bounds 
integer,parameter :: nz1=1,nz2=16              ! upper and lower bounds 

! For example:  4x4x4 grid of processors, each processor has
! 16x16x16 points + padding.  
!
! 2d Decomposition will look like:
! 64x16x4   (with pad=2 for FFT91:  66x16x4)
!
!  reshape3to2(p,pt,n1,n2,n3,index)
!
!  input:  p(nx,ny,nz)
!  output:  n1,n2,n3
!  output: pt(n1,n2,n3)   actual data: pt(n1-2,n2,n3) because FFT needs 
!                         padding of 2 on first dimension
!
!  we need to check that all possible 2D decompositions (n1,n2,n3)
!  n1*n2*n3 <= nx*ny*nz
!  Note: there are 3 possible 2D decompositions.  
!

! mesh coordinates
real*8 :: xcord(nx)
real*8 :: ycord(ny)
real*8 :: zcord(nz)

real*8  :: delt
integer :: itime_max,itime_output
parameter itime_max=1000     ! number of time steps
parameter itime_output=10    ! output every itime_output time steps

! set non-zero to cause time stepping loop to exit
! error_code = 1   advection CFL violated
! error_code = 2   sound speed CFL violated
! error_code = 3   diffusion CFL violated
integer :: error_code 
                      

end module mod_params








