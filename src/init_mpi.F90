#include "macros.h"


subroutine init_mpi
use params
use fft_interface
use mpi

implicit none
logical isperiodic(3),reorder
integer ierr1,ierr2,ierr3,rank
character*80 message

my_pe=0
mpicoords=0
mpidims=1
ioproc=0
initial_live_procs=1


#ifdef MPI

mpidims(1)=1
mpidims(2)=1
mpidims(3)=2



call mpi_init(ierr1)
if (ierr1<>0) call abort("mpi_init failure")
call mpi_comm_rank(MPI_COMM_WORLD,my_world_pe,ierr2)
if (ierr2<>0) call abort("mpi_comm_rank failure")
call mpi_comm_size(MPI_COMM_WORLD,initial_live_procs,ierr3)
if (ierr3<>0) call abort("mpi_comm_size failure")


isperiodic(1)=.false.
isperiodic(2)=.false.
isperiodic(3)=.false.
reorder=.false.

call mpi_cart_create(MPI_COMM_WORLD,3,mpidims,isperiodic,reorder,comm_3d,ierr1)
if (ierr1<>0) call abort("mpi_cart_create failure")

call mpi_cart_get(comm_3d,3,mpidims,isperiodic,mpicoords,ierr1)
if (ierr1<>0) call abort("mpi_cart_get failure")
print *,'me=',my_pe," mpi coords: ",mpicoords(1),mpicoords(2),mpicoords(3)

! get processor number with coords = mpicoords
call mpi_cart_rank(comm_3d,mpicoords,my_pe,ierr2)

#endif

print *, "me= ",my_pe,my_world_pe," total procs= ",initial_live_procs

end subroutine
