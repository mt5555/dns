#include "macros.h"
module particles
use params
implicit none
!
!  Particles module
!
!  all cpus have a copy of all particles.  But each particle "belongs"
!  to one pe, and the value on all other pe's will be 0.  
!
!  1. call init_particles with the number of particles desired
!  2. initialize particles::xparticle, yparticle, zparticle 
!  3. call advance_particles(rhs) to advance them one time step
!  
!
!  particle(:,1)  = xcord
!  particle(:,2)  = ycord
!  particle(:,3)  = zcord   (if ndim==3)
!  particle(:,ndim+1) = 'alf', a particle marker used for insertion
!
!  particle_old(:,1:ndim)   used for rk4 time stepping
!  particle_tmp(:,1:ndim)   used for rk4 time stepping
!
!  particle_work(:,1:ndim+1)  used for MPI buffer, other copy operations
!
integer :: nump=0             ! total number of particles
integer :: nump_max=0
real*8,allocatable :: particle(:,:)
real*8,private,allocatable :: particle_old(:,:)
real*8,private,allocatable :: particle_tmp(:,:)
real*8,private,allocatable :: particle_work(:,:)

integer,parameter :: ncross_max=5000
integer :: ncross=0
real*8  :: cross(ncross_max,2)


character(len=80),private :: message


contains

subroutine allocate_particles(in_nump)
implicit none
integer :: in_nump

nump=in_nump
nump_max=2*in_nump


if (allocated(particle)) then
   call abort("allocate_particles(): error: particles allready allocated")
endif

allocate(particle(nump_max,2*ndim+1))  
allocate(particle_work(nump_max,2*ndim+1))  
allocate(particle_old(nump_max,2*ndim))  
allocate(particle_tmp(nump_max,2*ndim))  

particle=-1d100
end subroutine




subroutine enlarge_particles()
!
! double the size of all particle_* arrays.
! preserve data in particle() array ONLY
!
implicit none
integer :: in_nump


nump_max=2*nump_max
call print_message("Increasig size of particle array")
write(message,'(a,i5)') "new size=",nump_max
call print_message(message)


particle_work=particle
deallocate(particle)
allocate(particle(nump_max,2*ndim+1))  
particle=-1d100
particle(1:nump,:)=particle_work(1:nump,:)

deallocate(particle_work)
allocate(particle_work(nump_max,2*ndim+1))  
deallocate(particle_old)
allocate(particle_old(nump_max,2*ndim))  
deallocate(particle_tmp)
allocate(particle_tmp(nump_max,2*ndim))  

end subroutine


subroutine particles_restart(fpe)
!
!  read==1   read particles from restart file
!  read==0   write particles to restart file
!
implicit none
integer :: fpe
real*8 :: time=0

call print_message("Restart particle data from file restart.particle")
call particles_io(1,fpe,'restart.particle')
write(message,'(a,i5)') 'total number of particles: ',nump
call print_message(message)
end subroutine





subroutine particles_save(fpe,time)
!
!  read==1   read particles from restart file
!  read==0   write particles to restart file
!
use params
implicit none
integer :: fpe
real*8 :: time
character(len=280) :: fname

if (nump==0) return

write(message,'(f10.4)') 10000.0000 + time
fname = rundir(1:len_trim(rundir)) // runname(1:len_trim(runname)) // &
     message(2:10) // ".particle"

call particles_io(0,fpe,fname)
end subroutine





subroutine particles_io(read,fpe,fname)
!
!  read==1   read particles from restart file
!  read==0   write particles to restart file
!
use params
use mpi
implicit none

integer :: read,fpe
character(len=*) :: fname

! local
CPOINTER :: fid
integer :: ierr,j,nump_in
real*8 :: x
character,save :: access="0"

if (my_pe==fpe) then

   if (read==1) then
      call copen(fname,"r",fid,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "particle_io(): Error opening file. Error no=",ierr
         call print_message(message)
         call print_message(fname)
         call abort("")
      endif
      call cread8e(fid,x,1,ierr)
      if (ierr/=1) then
         write(message,'(a,i5)') "particle_io(): Error reading file"
         call print_message(message)
         call print_message(fname)
         call abort("")
      endif
      nump_in=x

      call cread8(fid,x,1)
      if (2*ndim+1 /= nint(x)) then
         call abort("Error: ndim in code not the same as ndim in particles restart file")
      endif
      call allocate_particles(nump_in)
      do j=1,2*ndim+1
         call cread8(fid,particle(1,j),nump)
      enddo
      call cclose(fid,ierr)


   else
      call copen(fname,"w",fid,ierr)
      if (ierr/=0) then
         write(message,'(a,i5)') "particle_io(): Error opening file. Error no=",ierr
         call print_message(message)
         call print_message(fname)
         call abort("")
      endif
      x=nump
      call cwrite8(fid,x,1)
      x=2*ndim+1
      call cwrite8(fid,x,1)
      do j=1,2*ndim+1
         call cwrite8(fid,particle(1,j),nump)
      enddo
      call cclose(fid,ierr)

   endif
endif

#ifdef USE_MPI
if (read==1) then
   call MPI_bcast(nump,1,MPI_INTEGER,fpe,comm_3d ,ierr)
   if (my_pe/=fpe) then
      nump_in=nump
      call allocate_particles(nump_in)      
   endif
   call MPI_bcast(particle,(2*ndim+1)*nump_max,MPI_REAL8,fpe,comm_3d ,ierr)
endif
#endif



end subroutine


subroutine particle_advance(ugrid,rk4stage,time)
!
!  input:  psi  stream function
!          rk4stage  = 1,2,3,4 to denote which rk4 stage
!
!  work array: ugrid 
!
! stage 1:
!    particle_old=particle
!    particle_tmp=particle
!    RHS = interpolate psi to position of particle_tmp
!    particle = particle + delt*RHS/6
!    particle_tmp = particle_old + delt*RHS/2
!   
! stage 2:
!    RHS = interpolate pso to poisition of particle_tmp
!    particle = particle + delt*RHS/3
!    particle_tmp = particle_old + delt*RHS/2
!
! stage 3:
!    RHS = interpolate psi to position of particle_tmp
!    particle = particle + delt*RHS/3
!    particle_tmp = particle_old + delt*RHS
!    
! stage 4:
!    RHS = interpolate psi to position of particle_tmp
!    particle = particle + delt*RHS/6
!    particle_tmp = particle_old + 0*RHS (doesn't matter)
!    
!
!
use params
use mpi
use ghost
implicit none
real*8 :: ugrid(nx,ny,nz,2)
integer :: rk4stage
real*8 :: time

!local
integer :: i,j,ii,jj,igrid,jgrid
integer :: k,kk,kgrid
real*8 :: Ulocal(ndim),c1,c2
real*8 :: Q2d(4,4,ndim)
real*8 :: Q1d(4,ndim)
real*8 :: xc,yc,tmx1,tmx2
real*8 :: zc,tmx3
integer :: jc,kc
integer :: ierr

if (nump==0) return
call wallclock(tmx1)


if (rk4stage==1) then
   particle_tmp(:,1:2*ndim)=particle(:,1:2*ndim)
   particle_old(:,1:2*ndim)=particle(:,1:2*ndim)
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


! interpolate u to position in particle_tmp
do i=1,nump
   
   ! find position in global grid:
   igrid = 1 + floor( (particle_tmp(i,1)-g_xcord(1))/delx )
   jgrid = 1 + floor( (particle_tmp(i,2)-g_ycord(1))/dely )
   kgrid = 1 + floor( (particle_tmp(i,3)-g_zcord(1))/delz )
   
   if (1<=igrid .and. igrid+1<o_nx .and. 1<=jgrid .and. jgrid+1<o_ny &
        .and. 1<=kgrid .and. kgrid+1 < o_nz ) then
      
      ! compute a new point in the center of the above cell:
      ! (do this to avoid problems with 2 cpus both claiming a point
      ! on the boundary of a cell)
      xc=.5*(g_xcord(igrid)+g_xcord(igrid+1))
      yc=.5*(g_ycord(jgrid)+g_ycord(jgrid+1))
      zc=.5*(g_zcord(kgrid)+g_zcord(kgrid+1))
      
      ! find cpu which owns the grid point (xc,yc)
      if (xcord(intx1)<xc .and. xc<xcord(intx2)+delx .and. &
          ycord(inty1)<yc .and. yc<ycord(inty2)+dely .and. &
          zcord(intz1)<zc .and. zc<zcord(intz2)+delz ) then
         

         ! find igrid,jgrid,kgrid so that point (xc,yc,zc) is in box:
         ! igrid-1,igrid,igrid+1,igrid+2   and jgrid-1,jgrid,jgrid+1,jgrid+2
         ! and kgrid-1,kgrid,kgrid+1,kgrid+2   
         igrid = intx1 + floor( (xc-xcord(intx1))/delx )
         jgrid = inty1 + floor( (yc-ycord(inty1))/dely )
         kgrid = intz1 + floor( (zc-zcord(intz1))/delz )
         ASSERT("advance_particles(): igrid interp error",igrid<=intx2)
         ASSERT("advance_particles(): jgrid interp error",jgrid<=inty2)
         ASSERT("advance_particles(): kgrid interp error",kgrid<=intz2)
         
         ! interpolate Ulocal 
         ! 3d --> 2d y-z planes --> 1d z-lines --> point
         do kk=1,4
            do jj=1,4
               ! interpolate xcord(igrid-1:igrid+2) to xcord=particle(i,1)
               ! onto 2-d planes in y-z
               ! data  ugrid(igrid-1:igrid+2, jgrid-2+jj,:) 
               xc = 1 + (particle_tmp(i,1)-xcord(igrid))/delx
               jc = jgrid-2+jj
               kc = kgrid-2+kk
               !loop over velocity components
               do j=1,ndim
                  call interp4(ugrid(igrid-1,jc,kc,j),ugrid(igrid,jc,kc,j),&
                       ugrid(igrid+1,jc,kc,j),ugrid(igrid+2,jc,kc,j),&
                       xc,Q2d(jj,kk,j))
               enddo
            enddo
         enddo
         
         do kk=1,4
            ! interpolate ycord(jgrid-1:jgrid+2) to ycord=particle(j,1)
            ! onto lines in z
            ! data  ugrid(igrid-1:igrid+2, jgrid-2+jj,:) 
            yc = 1 + (particle_tmp(i,2)-ycord(jgrid))/dely
            !loop over velocity components
            do j=1,ndim
               call interp4(Q2d(1,kk,j),Q2d(2,kk,j),Q2d(3,kk,j),Q2d(4,kk,j),yc,Q1d(kk,j))
            enddo
         enddo
         
         ! interpolate zcord(kgrid-1:kgrid+2) to zcord=particle(k,2)
         ! data:  Q1d(1:4,j)
!
! BUG???
! JAMAL: should this be particle_tmp(i,3) (the z-coordinate) 
!        in the next line?  
!
!
         zc = 1 + (particle_tmp(i,2)-zcord(kgrid))/delz
         !loop over velocity components
         do j=1,ndim
            call interp4(Q1d(1,j),Q1d(2,j),Q1d(3,j),Q1d(4,j),zc,Ulocal(j))
         enddo

         !for inertial particles, change this step to use appropriate drag law
         ! first ndim data are positions
         ! next ndim data are velocities
         !update velocity first
         ! advance
         do j=1,ndim
            particle(i,j+ndim)=Ulocal(j)
            particle_tmp(i,j+ndim)=Ulocal(j)
            particle(i,j)=particle(i,j)+c1*particle(i,j+ndim)
            particle_tmp(i,j)=particle_old(i,j)+c2*particle(i,j+ndim)
         enddo
      else
         ! this cpu does not own the particle
         ! point does not belong to my_pe, set position to -inf
         do j=1,2*ndim
            particle(i,j)=-1d100
            particle_tmp(i,j)=-1d100
         enddo
      endif
   else
      ! partical has left the domain
      call abort("paritcal has left computational domain")
   endif
enddo



#ifdef USE_MPI
   particle_work=particle
   call MPI_allreduce(particle_work,particle,(2*ndim+1)*nump_max,MPI_REAL8,MPI_MAX,comm_3d,ierr)
   if (rk4stage/=4) then
      ! not necessary on last stage - we no longer need particle_tmp
      particle_work(:,1:2*ndim)=particle_tmp(:,1:2*ndim)
      call MPI_allreduce(particle_work,particle_tmp,2*ndim*nump_max,MPI_REAL8,MPI_MAX,comm_3d,ierr)
   endif
#endif

!
! make sure particle and particle_tmp have not left the domain
do i=1,nump
  if (particle(i,1)>g_xcord(g_nx)) particle(i,1)=particle(i,1)-xscale
  if (particle(i,1)<g_xcord(1)) particle(i,1)=particle(i,1)+xscale
  if (particle(i,2)>g_ycord(g_ny)) particle(i,2)=particle(i,2)-yscale
  if (particle(i,2)<g_ycord(1)) particle(i,2)=particle(i,2)+yscale
  if (particle(i,3)>g_zcord(g_nz)) particle(i,3)=particle(i,3)-zscale
  if (particle(i,3)<g_zcord(1)) particle(i,3)=particle(i,3)+zscale

  if (particle_tmp(i,1)>g_xcord(g_nx)) particle_tmp(i,1)=particle_tmp(i,1)-xscale
  if (particle_tmp(i,1)<g_xcord(1)) particle_tmp(i,1)=particle_tmp(i,1)+xscale
  if (particle_tmp(i,2)>g_ycord(g_ny)) particle_tmp(i,2)=particle_tmp(i,2)-yscale
  if (particle_tmp(i,2)<g_ycord(1)) particle_tmp(i,2)=particle_tmp(i,2)+yscale
  if (particle_tmp(i,3)>g_zcord(g_nz)) particle_tmp(i,3)=particle_tmp(i,3)-zscale
  if (particle_tmp(i,3)<g_zcord(1)) particle_tmp(i,3)=particle_tmp(i,3)+zscale
enddo


if (rk4stage==4) then
   ! xcord set to -1d100 to denote off processor.  
   ! if any particle has left *all* processors, abort:
   if (minval(particle(1:nump,1))<g_xcord(1)) then
      do i=1,nump
         write(*,'(i5,3f12.5)') i,particle(i,1),particle(i,2),particle(i,3)
      enddo
      call abort("particle_advance(): point has left domain") 
   endif



   ! insert points into particle() if necessary:
   if (nump>(3*nump_max/4)) call enlarge_particles()

   ! INSERT PARTICLES HERE IF NEEDED


endif



call wallclock(tmx2)
tims(16)=tims(16)+(tmx2-tmx1)
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
