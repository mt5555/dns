subroutine init_mpi
use params
use fft_interface
implicit none
#ifdef MPI
integer ierr1,ierr2,ierr3
#include "mpif.h"
#endif
character*80 message

ncpu_x=1
ncpu_y=1
ncpu_z=1

myproc=0
myproc_x=0
myproc_y=0
myproc_z=0

ioproc=0


#ifdef MPI
call mpi_init(ierr1)
call mpi_comm_rank(MPI_COMM_WORLD,myproc,ierr2)
call mpi_comm_size(MPI_COMM_WORLD,initial_live_procs,ierr3)

if(debug_mpi.ne.0) then 

write(message,*) "me= ",me_global,"total procs= ",initial_live_procs,"Error codes= ",ierr1,ierr2,ierr3
call print_message(message)

end if 
#endif


end subroutine
