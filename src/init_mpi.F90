#include "macros.h"


subroutine init_mpi
use params
use mpi

implicit none
integer ierr1,ierr2,ierr3,rank
character*80 message

my_pe=0
mpicoords=0
mpidims=1
io_pe=0
initial_live_procs=1

ncpu_x=1
ncpu_y=1
ncpu_z=1
my_x=0
my_y=0
my_z=0



#ifdef USE_MPI
call mpi_init(ierr1)
if (ierr1/=0) call abort("mpi_init failure")
call mpi_comm_rank(MPI_COMM_WORLD,my_world_pe,ierr2)
if (ierr2/=0) call abort("mpi_comm_rank failure")
call mpi_comm_size(MPI_COMM_WORLD,initial_live_procs,ierr3)
if (ierr3/=0) call abort("mpi_comm_size failure")
#endif

if (ncpu_x * ncpu_y * ncpu_z /= initial_live_procs) then
   call print_message("Error: incorrect number of cpus");

   write(message,'(a,i5,a,i3,a,i3,a,i3)') "Parallel decomposition requested: ncpus= ", &
      ncpu_x*ncpu_y*ncpu_z," = ",ncpu_x," x",ncpu_y," x",ncpu_z
   call print_message(message)
   write(message,'(a,i5)') "Actual ncpus = ",initial_live_procs
   call print_message(message)
   call abort("Terminating.")
endif
end subroutine





subroutine init_mpi_comm3d
use params
use fft_interface
use mpi

implicit none
logical isperiodic(3),reorder
integer ierr1,ierr2,ierr3,rank
character*80 message


#ifdef USE_MPI
isperiodic(1)=.false.
isperiodic(2)=.false.
isperiodic(3)=.false.
reorder=.false.

call mpi_cart_create(MPI_COMM_WORLD,3,mpidims,isperiodic,reorder,comm_3d,ierr1)
if (ierr1/=0) call abort("mpi_cart_create failure")

call mpi_cart_get(comm_3d,3,mpidims,isperiodic,mpicoords,ierr1)
if (ierr1/=0) call abort("mpi_cart_get failure")
!print *,'me=',my_pe," mpi coords: ",mpicoords(1),mpicoords(2),mpicoords(3)

! get processor number with coords = mpicoords
call mpi_cart_rank(comm_3d,mpicoords,my_pe,ierr2)

write(message,'(a,i5,a,i3,a,i3,a,i3)') "Parallel decomposition: ncpus= ", &
   initial_live_procs," = ",ncpu_x," x",ncpu_y," x",ncpu_z
call print_message(message)
#else
call print_message("Running single threaded")
#endif
end subroutine




subroutine close_mpi
use params
use mpi
integer ierr

#ifdef USE_MPI
call MPI_Finalize(ierr)
#endif

end subroutine
