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


subroutine wallclock(tmx)
use params
implicit none
integer count,count_rate,count_max
real*8 tmx
#ifdef MPI
real*8 MPI_Wtime
tmx = MPI_Wtime()
#else
call system_clock(count,count_rate,count_max)
tmx=real(count,r8kind)/real(count_rate,r8kind)
#endif

end subroutine

