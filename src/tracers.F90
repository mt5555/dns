#include "macros.h"
module tracers
use params
implicit none
!
!  Tracers module
!
!  all cpus have a copy of all tracers.  But each tracer "belongs"
!  to one pe, and the value on all other pe's will be 0.  
!
!  1. call init_tracers with the number of tracers desired
!  2. initialize tracers::xtracer, ytracer, ztracer 
!  3. call advance_tracers(rhs) to advance them one time step
!  
!
!
integer :: numt=0
real*8,allocatable :: tracer(:,:)
real*8,private,allocatable :: tracer_old(:,:)
real*8,private,allocatable :: tracer_tmp(:,:)
real*8,private,allocatable :: tracer_work(:,:)
real*8,private,allocatable :: tracer_rhs(:,:)



contains

subroutine allocate_tracers(in_numt)
implicit none
integer :: in_numt


numt=in_numt

allocate(tracer(numt,n_var))  
allocate(tracer_work(numt,n_var))  
allocate(tracer_old(numt,n_var))  
allocate(tracer_tmp(numt,n_var))  
allocate(tracer_rhs(numt,n_var))  

end subroutine


subroutine tracers_restart(fpe)
!
!  read==1   read tracers from restart file
!  read==0   write tracers to restart file
!
implicit none
integer :: fpe
real*8 :: time=0
call tracers_io(1,fpe,'restart.tracer')
end subroutine


subroutine tracers_save(fpe,time)
!
!  read==1   read tracers from restart file
!  read==0   write tracers to restart file
!
use params
implicit none
integer :: fpe
real*8 :: time
character(len=280) :: fname,message

if (numt==0) return

write(message,'(f10.4)') 10000.0000 + time
fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // &
     message(2:10) // ".tracer"

call tracers_io(0,fpe,fname)
end subroutine





subroutine tracers_io(read,fpe,fname)
!
!  read==1   read tracers from restart file
!  read==0   write tracers to restart file
!
use params
use mpi
implicit none

integer :: read,fpe
character(len=*) :: fname

! local
character(len=280) :: message
CPOINTER :: fid
integer :: ierr
real*8 :: x

if (my_pe==fpe) then

   if (read==1) then
      call copen(fname,"r",fid,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "tracer_io(): Error opening file. Error no=",ierr
         call print_message(message)
         call print_message(fname)
         call abort("")
      endif
      call cread8e(fid,x,1,ierr)
      if (ierr/=1) then
         write(message,'(a,i5)') "tracer_io(): Error reading file"
         call print_message(message)
         call print_message(fname)
         call abort("")
      endif
      numt=x
      
      call cread8(fid,x,1)
      if (n_var /= nint(x)) then
         call abort("Error: n_var in code not the same as n_var in tracers restart file")
      endif
      call allocate_tracers(numt)
      call cread8(fid,tracer,numt*n_var)
      call cclose(fid,ierr)
   else
      call copen(fname,"w",fid,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "tracer_io(): Error opening file. Error no=",ierr
         call print_message(message)
         call print_message(fname)
         call abort("")
      endif
      x=numt
      call cwrite8(fid,x,1)
      x=n_var
      call cwrite8(fid,x,1)
      call cwrite8(fid,tracer,numt*n_var)
      call cclose(fid,ierr)
   endif
endif

end subroutine










subroutine tracer_advance(psi,rk4stage)
!
! stage 1:
!    tracer_old=tracer
!    tracer_tmp=tracer
!    RHS = interpolate psi to position of tracer_tmp
!    tracer = tracer + delt*RHS/6
!    tracer_tmp = tracer_old + delt*RHS/2
!   
! stage 2:
!    RHS = interpolate pso to poisition of tracer_tmp
!    tracer = tracer + delt*RHS/3
!    tracer_tmp = tracer_old + delt*RHS/2
!
! stage 3:
!    RHS = interpolate psi to position of tracer_tmp
!    tracer = tracer + delt*RHS/3
!    tracer_tmp = tracer_old + delt*RHS
!    
! stage 4:
!    RHS = interpolate psi to position of tracer_tmp
!    tracer = tracer + delt*RHS/6
!    tracer_tmp = tracer_old + 0*RHS (doesn't matter)
!    
!
!
use params
use mpi
implicit none
real*8 :: psi(nx,ny)
integer :: rk4stage

!local
integer :: i,j,ii,jj
real*8 :: trhs(n_var),c1,c2
integer :: ierr

if (numt==0) return

if (rk4stage==1) then
   tracer_tmp=tracer
   tracer_old=tracer
   c1=delt/6
   c2=delt/2
else if (rk4stage==2) then
   c1=delt/3
   c2=delt/2
else if (rk4stage==3) then
   c1=delt/3
   c2=delt
else if (rk4stage==4) then
   c1=delt/6
   c2=0
endif



! interpolate psi to position in tracer_tmp
do i=1,numt
   ! find cpu which owns the grid point tracer(i,:)
   if (xcord(bx1)<=tracer_tmp(i,1) .and. tracer_tmp(i,1)<xcord(bx2) .and. &
       ycord(by1)<=tracer_tmp(i,2) .and. tracer_tmp(i,2)<ycord(by2) ) then

      ! find ii,jj so that point is in box:
      ! ii,ii+1,ii+2,ii+3   and jj,jj+1,jj+2,jj+3
      ! interpolate trhs
      do j=1,3
         tracer(i,j)=tracer(i,j)+c1*trhs(j)
         tracer_tmp(i,j)=tracer_old(i,j)+c2*trhs(j)
      enddo
   else
      ! point does not belong to my_pe, set xt=yt=xtrhs=ytrhs=0
      do j=1,3
         tracer(i,j)=0
         tracer_tmp(i,j)=0
      enddo
   endif
enddo



#ifdef USE_MPI
   tracer_work=tracer
   call MPI_allreduce(tracer_work,tracer,n_var*numt,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   tracer_work=tracer_tmp
   call MPI_allreduce(tracer_work,tracer_tmp,n_var*numt,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif

end subroutine






end module
