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


!
! Linux	
!
#ifdef LINUX
module mpi
#include "mpif.h"

#ifndef MPI_HAS_REAL8
integer MPI_REAL8 
parameter (MPI_REAL8=MPI_DOUBLE_PRECISION)
integer MPI_REAL4 
parameter (MPI_REAL4=MPI_REAL)
#endif

end module
#endif


!
! OSF1
!
#ifdef OSF1
module mpi
#include "mpif.h"
#if (defined USE_MPI_IO && defined USE_LAMPI)
#include "mpiof.h"
#endif
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







