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










subroutine tracer_advance(psi,ugrid,rk4stage)
!
!  input:  psi  stream function
!          rk4stage  = 1,2,3,4 to denote which rk4 stage
!
!  work array: ugrid 
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
use ghost
implicit none
real*8 :: psi(nx,ny)
real*8 :: ugrid(nx,ny,2)
integer :: rk4stage

!local
integer :: i,j,ii,jj,igrid,jgrid
real*8 :: trhs(n_var),c1,c2
real*8 :: Qint(4,n_var)
real*8 :: xc,yc
integer :: jc
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


do j=inty1,inty2
do i=intx1,intx2
      
      ugrid(i,j,1)=( 2*(psi(i,j+1)-psi(i,j-1))/3 -  &
           (psi(i,j+2)-psi(i,j-2))/12          )/dely
      
      ugrid(i,j,2)=-( 2*(psi(i+1,j)-psi(i-1,j))/3 -  &
           (psi(i+2,j)-psi(i-2,j))/12          )/delx
      
enddo
enddo

call ghost_update_x(ugrid,2)
call ghost_update_y(ugrid,2)



! interpolate psi to position in tracer_tmp
do i=1,numt
   ! find cpu which owns the grid point tracer(i,:)
   if (xcord(intx1)<=tracer_tmp(i,1) .and. tracer_tmp(i,1)<xcord(intx2) .and. &
       ycord(inty1)<=tracer_tmp(i,2) .and. tracer_tmp(i,2)<ycord(inty2) ) then

      ! find igrid,jgrid so that point is in box:
      ! igrid-1,igrid,igrid+1,igrid+2   and jgrid-1,jgrid,jgrid+1,jgrid+2
      igrid = intx1 + floor( (tracer_tmp(i,1)-xcord(intx1))/delx )
      jgrid = inty1 + floor( (tracer_tmp(i,2)-ycord(inty1))/dely )
      ASSERT("advance_tracers(): igrid interp error",igrid<intx2)
      ASSERT("advance_tracers(): jgrid interp error",jgrid<inty2)


      ! interpolate trhs
      do jj=1,4
         ! interpolate xcord(igrid-1:igrid+2) to xcord=tracer(i,1)
         ! data  ugrid(igrid-1:igrid+2, jgrid-2+jj,:) 
         xc = 1 + (tracer_tmp(i,1)-xcord(igrid))/delx
         jc = jgrid-2+jj
         do j=1,ndim
            call interp4(ugrid(igrid-1,jc,j),ugrid(igrid,jc,j),&
                ugrid(igrid+1,jc,j),ugrid(igrid+2,jc,j),&
                xc,Qint(jj,j))
         enddo
      enddo
      ! interpolate ycord(jgrid-1:jgrid+2) to ycord=tracer(i,2)
      ! data:  Qint(1:4,j)
      yc = 1 + (tracer_tmp(i,2)-ycord(jgrid))/dely
      do j=1,ndim
         call interp4(Qint(1,j),Qint(2,j),Qint(3,j),Qint(4,j),yc,trhs(j))
      enddo


      ! advance
      do j=1,ndim
         tracer(i,j)=tracer(i,j)+c1*trhs(j)
         tracer_tmp(i,j)=tracer_old(i,j)+c2*trhs(j)
      enddo
   else
      ! point does not belong to my_pe, set position to 0
      do j=1,ndim
         tracer(i,j)=0
         tracer_tmp(i,j)=0
      enddo
   endif
enddo



#ifdef USE_MPI
   tracer_work=tracer
   call MPI_allreduce(tracer_work,tracer,n_var*numt,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   if (rk4sgate/=4) then
      ! not necessary on last stage - we no longer need tracer_tmp
      tracer_work=tracer_tmp
      call MPI_allreduce(tracer_work,tracer_tmp,n_var*numt,MPI_REAL8,MPI_SUM,comm_3d,ierr)
   endif
#endif



if (rk4stage==4) then
   ! insert points into tracer() if necessary:

endif





end subroutine






      SUBROUTINE interp4(y0,y1,y2,y3,newx,ynew)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     Interpolates "y" to xloc position
!
!
!     y0,y1,y2,y3 is data specified at points 0,1,2,3
!     newx should be a point between 0 and 3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer i,n,xbeg
      real*8 x0,x1,x2,x3,y0,y1,y2,y3,ynew
      real*8 denom0,denom1,denom2,denom3,fact0,fact1,fact2,fact3,newx

      x0 = 0
      x1 = 1
      x2 = 2
      x3 = 3

      denom0 = -6  !(x0-x1)*(x0-x2)*(x0-x3)   ! (-1)(-2)(-3)=-6
      denom1 =  2  !(x1-x0)*(x1-x2)*(x1-x3)   ! ( 1)(-1)(-2)= 2
      denom2 = -2  !(x2-x0)*(x2-x1)*(x2-x3)   ! ( 2)( 1)(-1)=-2
      denom3 =  6  !(x3-x0)*(x3-x1)*(x3-x2)   ! ( 3)( 2)( 1)= 6

      fact0 = (newx-x1)*(newx-x2)*(newx-x3)/denom0
      fact1 = (newx-x0)*(newx-x2)*(newx-x3)/denom1
      fact2 = (newx-x0)*(newx-x1)*(newx-x3)/denom2
      fact3 = (newx-x0)*(newx-x1)*(newx-x2)/denom3
      
      ynew = y0*fact0 + y1*fact1 + y2*fact2 + y3*fact3

      return
      end subroutine





end module
