#include "macros.h"

subroutine abort(message)
use mpi
implicit none
character*(*) message
character*15 :: pre="ASSERT FAILURE "
call print_message(pre // message)
#ifdef MPI
   call MPI_abort()
#endif

end subroutine


subroutine print_message(message)
use params
implicit none
character*(*) message
if (my_pe==io_pe) then
   write(*,'(a)') message
endif
end subroutine




