subroutine init_model
use params
use fft_interface
implicit none


call params_init()
call fft_interface_init()

end subroutine