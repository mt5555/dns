#include "macros.h"

subroutine abort(message)
implicit none
character*(*) message
character*15 :: pre="ASSERT FAILURE "
call print_message(pre // message)
stop
end subroutine


subroutine print_message(message)
implicit none
character*(*) message
write(*,'(a)') message
end subroutine




