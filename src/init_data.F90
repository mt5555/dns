subroutine init_data(Q)
use params
use fft_interface
implicit none
real*8 :: Q(nx,ny,nz,n_var)

! uniform flow to the right
Q=0
Q(:,:,:,1)=1



call bc_preloop


end subroutine