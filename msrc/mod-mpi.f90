Module mod_mpi
implicit none
include "mpif.h"
integer :: me_global, initial_live_procs
integer :: comm_3d, me_cart
integer, parameter :: ndims=3
integer,  dimension(ndims) :: dims
logical, parameter, dimension(ndims) :: periodic=.true.
logical, parameter :: reorganize=.true.
integer :: ierr




!      |z
!      |up /y
!      |  /west
!      | /
!      |/_________ x nort(h) 


integer :: n_uppp,n_down,n_west,n_east,n_nort,n_sout

DOUBLE PRECISION :: MPI_WTIME
external MPI_WTIME





end module mod_mpi
