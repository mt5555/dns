#ifdef USE_MPI


module mpi
#include "mpif.h"

#ifndef MPI_HAS_REAL8
   integer MPI_REAL8 
   parameter (MPI_REAL8=MPI_DOUBLE_PRECISION)
   integer MPI_REAL4 
   parameter (MPI_REAL4=MPI_REAL)
#endif

end module



#else
!
! dummy module, not using mpi, but code can still have "use mpi" statment
!
module mpi
integer dummy_mpi_variable
end module

#endif







