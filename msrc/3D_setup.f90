program 3D_setup
use mod_mpi
use mod_debug

integer :: me_global2
integer :: live_procs
integer :: old_comm

call mpi_comm_rank(MPI_COMM_WORLD,me_global2,ierr2)
call mpi_comm_size(MPI_COMM_WORLD,live_procs,ierr3)

old_com=MPI_COMM_WORD
call mpi_cart_create(old_comm,ndims,dims,periodic,reorganize,  &
                     comm_3d,ierr4)

call mpi_cart_rank(MPI_COMM_WORLD,me_cart,ierr2)

if(debug_mpi.ne.0) then 
ierr=ierr4+iedd3+ierr2
if(ierr>)) write(*,*) "PDNS ERROR in 3D_setup",ierr4,iedd3,ierr2

if(me_cart.ne.me_global2) write(*,*)"PDNS note on reorder in 3D_setup",me_global,me_global2,me_cart 

end if

end program 3D_setup

