Module mod-mpi
implicit none
include "mpif.h"
integer :: me_global, initial_live_procs
integer :: init_comm_3d
integer :: ierr



!      |z
!      |up /y
!      |  /west
!      | /
!      |/_________ x nort(h) 


integer :: n_uppp,n_down,n_west,n_east,n_nort,n_sout

real*8 :: MPI_WTIME
external MPI_WTIME





end module mod-mpi
