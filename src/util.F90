#include "macros.h"

subroutine abort(message)
use mpi
use params
implicit none
integer ierr
character*(*) message
character*15 :: pre="ASSERT FAILURE "
write(*,'(a)') pre // message

#ifdef IBM
call flush_(6)
#else
call flush(6)
#endif

#ifdef USE_MPI
   call MPI_abort(comm_3d,1,ierr)
#endif
stop
end subroutine


subroutine print_message(message)
use params
implicit none
character*(*) message
if (my_pe==io_pe) then
   write(*,'(a)') message
endif
end subroutine




