subroutine init_model
use params
use fft_interface
implicit none

real*8 :: one=1
pi=4*atan(one)

call params_init()
call fft_interface_init()

end subroutine