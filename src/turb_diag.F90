#include "macros.h"
subroutine output_model(time,Q,q1,q2,q3,work1,work2)
use params
use structf
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: q1(nx,ny,nz,n_var)
real*8 :: q2(nx,ny,nz,n_var)
real*8 :: q3(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 :: work2(nx,ny,nz)
real*8 :: time

! local variables
integer,parameter :: nints_e=14
real*8 :: ints_e(nints_e)
real*8 :: x
integer i,j,k,n,ierr
character(len=80) :: message
character :: access
CPOINTER fid

access="a"
if (time==0) access="w"




!
! output structure functions
!
if (compute_struct==1) then
   call compute_all_pdfs(Q,q1,q2,q3,work1,work2,ints_e,nints_e)
   
   
   write(message,'(a,3f14.8)') 'skewness ux,vw,wz: ',&
        (ints_e(n+3)/ints_e(n)**1.5,n=1,3)
   call print_message(message)
   
   write(message,'(a,3f14.8)') 'wSw: ',&
        (ints_e(10)/ints_e(1)**2)
   call print_message(message)
   
   
   if (structf_init==1) then
   if (my_pe==io_pe) then
      write(message,'(f10.4)') 10000.0000 + time_initial
      message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".sf"
      call copen(message,access,fid,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "outputSF(): Error opening file errno=",ierr
         call abort(message)
      endif
   endif
   call outputSF(time,fid)
   if (my_pe==io_pe) call cclose(fid,ierr)
   endif


if (my_pe==io_pe) then
   write(message,'(f10.4)') 10000.0000 + time_initial
   message = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // message(2:10) // ".scalars-turb"
   call copen(message,access,fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "diag_output(): Error opening .scalars-turb file errno=",ierr
      call abort(message)
   endif
   x=nints_e; call cwrite8(fid,x,1)
   call cwrite8(fid,time,1)
   call cwrite8(fid,ints_e,nints_e)


   call cclose(fid,ierr)
   if (ierr/=0) then
      write(message,'(a,i5)') "diag_output(): Error closing .scalars-turb file errno=",ierr
      call abort(message)
   endif
endif

endif

end subroutine

