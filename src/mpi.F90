#ifdef USE_MPI



!
! SUN
!
#ifdef SUNOS
#include "/usr/local/src/mpich/include/mpif.f90"
#endif



!
! IRIX
!
#ifdef IRIX64
module mpi
#include "mpif.h"
end module
#endif





#else
!
! dummy module, not using mpi, but code can still have "use mpi" statment
!
module mpi
integer dummy_mpi_variable
end module

#endif







