#include "macros.h"
module subcubes
use params
use mpi
implicit none
#if 0

module to extract subcubes and apply arbritrary rotation

ssize = (input) size (in gridpoints) of the subcubes

nsubcube  = total number of subcubes
subcube_size = length of the subcube
wsubcube_size = length of the larger "window" subcube, centered over the
                original subcube

subcube_corner(3,i)   x,y,z coordinates of the i-th subcubes corner
wsubcube_corner(3,i)   x,y,z coordinates of the i-th window subcube
                

subcube_cords(1:ssize)     base coordinates.  coordinates of ith subcube
                           are given by:
                               x(:)=subcube_corner(1,i)+subcube_cords(:)
                               y(:)=subcube_corner(2,i)+subcube_cords(:)
 


#endif

integer :: nsubcube
real*8,allocatable :: subcube_corner(:,:)
real*8,allocatable :: wsubcube_corner(:,:)
real*8 :: subcube_size
real*8 :: wsubcube_size

real*8,allocatable :: subcube_cords(:)   
real*8,allocatable :: wsubcube_cords(:)   



integer,private   :: init = 0

contains


subroutine setup_subcubes(ssize)
!
!  ssize = 128   cubes of size 128^3
!  diff  = 1/8   spacing between cubes
!
!
integer :: ssize,i
real*8 :: diff,ax,ay,az

subcube_size=ssize*delx
wsubcube_size=2*subcube_size   ! 2x larger
diff = ssize*delx/2  ! overlap 50%


if (allocated(subcube_cords)) deallocate(subcube_cords)
if (allocated(wsubcube_cords)) deallocate(wsubcube_cords)
allocate(subcube_cords(ssize))
allocate(wsubcube_cords(ssize*2))


do i=1,ssize
   subcube_cords(i)=(i-1)*delx
enddo
do i=1,2*ssize
   wsubcube_cords(i)=(i-1)*delx - subcube_size/2
enddo


nsubcube=0
ax=0
do while (ax<1)
   ay=0
   do while (ay<1)
      az=0
      do while (az<1)
         nsubcube=nsubcube+1
         az=az+diff
      enddo
      ay=ay+diff
   enddo
   ax=ax+diff
enddo


if (allocated(subcube_corner)) deallocate(subcube_corner)
if (allocated(wsubcube_corner)) deallocate(wsubcube_corner)

allocate(subcube_corner(3,nsubcube))
allocate(wsubcube_corner(3,nsubcube))

i=0
ax=0
do while (ax<1)
   ay=0
   do while (ay<1)
      az=0
      do while (az<1)
         i=i+1
         subcube_corner(1,i)=ax
         subcube_corner(2,i)=ay
         subcube_corner(3,i)=az
         wsubcube_corner(1,i)=ax - subcube_size/2
         wsubcube_corner(2,i)=ay - subcube_size/2
         wsubcube_corner(3,i)=az - subcube_size/2
         az=az+diff
      enddo
      ay=ay+diff
   enddo
   ax=ax+diff
enddo

end subroutine



subroutine interp_subcube(nsx,nsy,nsz,dx,dy,dz,data,field,irot,rmatrix)
!
! interpolate to the points dx(1:ns),dy(1:ns),dz(1:ns)
!
! if irot=0, then interpolate by finding closest grid point.
!
! if irot=1, apply rotation matrix and then do true interpolation
!
! NOTE: only io_pe will have the correct, interpolated data
!
integer :: nsx,nsy,nsz
real*8 :: dx(nsx),dy(nsy),dz(nsz)
real*8 :: data(nsx,nsy,nsz)
real*8 :: field(nx,ny,nz)
integer :: irot
real*8 :: rmatrix(3,3,3)

real*8 :: xi,yi,zi
integer :: n1(nsx),n2(nsy),n3(nsz),i,j,k,igrid,jgrid,kgrid
#ifdef USE_MPI
real*8,allocatable :: data2(:,:,:)
integer :: ierr
#endif

!if (io_pe==my_pe) then
!   print *,'interpolating to subcube: '
!   print *,dx(1),dy(1),dz(1)
!   print *,dx(nsx),dy(nsy),dz(nsz)
!endif

do i=1,nsx
   xi=dx(i)
   if (xi<g_xcord(1)) xi=xi+1;
   if (xi>g_xcord(g_nx)) xi=xi-1;
   n1(i)=(xi-xcord(nx1))/delx
   n1(i)=n1(i)+nx1
   if (n1(i)<nx1 .or. n1(i)>nx2) n1(i)=0 !flag indicating not on this cpu
enddo
do j=1,nsy
   yi=dy(j)
   if (yi<g_ycord(1)) yi=yi+1;
   if (yi>g_ycord(g_ny)) yi=yi-1;
   n2(j)=(yi-ycord(ny1))/dely
   n2(j)=n2(j)+ny1
   if (n2(j)<ny1 .or. n2(j)>ny2) n2(j)=0 !flag indicating not on this cpu
enddo
do k=1,nsz
   zi=dz(k)
   if (zi<g_zcord(1)) zi=zi+1;
   if (zi>g_zcord(g_nz)) zi=zi-1;
   n3(k)=(zi-zcord(nz1))/delz
   n3(k)=n3(k)+nz1
   if (n3(k)<nz1 .or. n3(k)>nz2) n3(k)=0 !flag indicating not on this cpu
enddo


do k=1,nsz
   do j=1,nsy
      do i=1,nsx

         if (irot==1) then
            xi=dx(i)
            yi=dy(j)
            zi=dz(k)
            !rotate to (xi,yi,zi), and map back into [0,1]
            !call rotate(xi,yi,zi,rmatrix)
            ! find cell containing (xi,yi,zi):

            ! find position in global grid:
            igrid = 1 + floor( (xi-g_xcord(1))/delx )
            jgrid = 1 + floor( (yi-g_ycord(1))/dely )
            kgrid = 1 + floor( (zi-g_zcord(1))/delz )

            ! compute a new point in the center of the above cell:
            ! (do this to avoid problems with 2 cpus both claiming a point
            ! on the boundary of a cell)
            xi=.5*(g_xcord(igrid)+g_xcord(igrid+1))
            yi=.5*(g_ycord(jgrid)+g_ycord(jgrid+1))
            zi=.5*(g_zcord(kgrid)+g_zcord(kgrid+1))

            ! find cpu which owns this cell (xi,yi,zi):
            igrid=(xi-xcord(nx1))/delx
            jgrid=(yi-ycord(ny1))/dely
            kgrid=(zi-zcord(nz1))/delz

            ! interpolate:
            
         else
            if (n1(i)==0 .or. n2(j)==0 .or. n3(k)==0) then
               data(i,j,k)=-9d200
            else
               ! we have this data, copy;
               data(i,j,k)=field(n1(i),n2(j),n3(k))
            endif
         endif
         
      enddo
   enddo
enddo
#ifdef USE_MPI
allocate(data2(nsx,nsy,nsz))
data2=data
call mpi_reduce(data2,data,nsx*nsy*nsz,MPI_REAL8,MPI_MAX,io_pe,comm_3d,ierr)
deallocate(data2)
#endif

#if 0
do k=1,nsz
   do j=1,nsy
      do i=1,nsx

         xi=dx(i)
         yi=dy(j)
         zi=dz(k)
         !if (irot==1) call rotate(xi,yi,zi,rmatrix)

         n3=(zi-zcord(nz1))/delz
         n3=n3+nz1
         n2=(yi-ycord(ny1))/dely
         n2=n2+ny1
         n1=(xi-xcord(nx1))/delx
         n1=n1+nx1

         if (data(i,j,k)<-9d199) then
            print *,my_pe,n1,n2,n3
         endif
      enddo
   enddo
enddo
#endif

end subroutine

end module
