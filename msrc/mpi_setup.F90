program mpi_setup
use mod-mpi
use mod-debug

call mpi_init(ierr1)
call mpi_comm_rank(MPI_COMM_WORLD,me_global,ierr2)
call mpi_comm_size(MPI_COMM_WORLD,initial_live_procs,ierr3)

if(debug_mpi.ne.0) then 

write(*,*) "me= ",me_global,"total procs= ",initial_live_procs,"Error codes= ",ierr1,ierr2,ierr3

end if 

end program mpi_setup 
