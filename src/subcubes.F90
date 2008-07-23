!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Copyright 2007.  Los Alamos National Security, LLC. This material was
!produced under U.S. Government contract DE-AC52-06NA25396 for Los
!Alamos National Laboratory (LANL), which is operated by Los Alamos
!National Security, LLC for the U.S. Department of Energy. The
!U.S. Government has rights to use, reproduce, and distribute this
!software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
!LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY
!FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
!derivative works, such modified software should be clearly marked, so
!as not to confuse it with the version available from LANL.
!
!Additionally, this program is free software; you can redistribute it
!and/or modify it under the terms of the GNU General Public License as
!published by the Free Software Foundation; either version 2 of the
!License, or (at your option) any later version. Accordingly, this
!program is distributed in the hope that it will be useful, but WITHOUT
!ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
!for more details.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!diff = ssize*delx     ! no overlap
!diff = ssize*delx/2  ! overlap 50%

! 12x12x12 subcubes
i=g_nx/12
if (mod(g_nx,12)/=0) i=i+1
diff = i*delx    


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



subroutine interp_subcube(ssize,k1,k2,x0,e01,e02,e03,data,field,irot)
!
! interpolate to the points x+dx(1:ns),y+dy(1:ns),z+dz(1:ns)
!
! NOTE: only io_pe will have the correct, interpolated data
!
integer :: ssize,k1,k2,irot
real*8 :: data(ssize,ssize,k2-k1+1)
real*8 :: field(nx,ny,nz)
real*8 :: x0(3),e01(3),e02(3),e03(3)


real*8 :: xi(3),err
integer :: i,j,k,ii,jj,kk
#ifdef USE_MPI
real*8,allocatable :: data2(:,:,:)
integer :: ierr
#endif

!if (io_pe==my_pe) then
!   print *,'interpolating to subcube: '
!   print *,dx(1),dy(1),dz(1)
!   print *,dx(nsx),dy(nsy),dz(nsz)
!endif

#define CHECK_IERR


do k=k1,k2
   do j=1,ssize
      do i=1,ssize

         xi = x0 + (i-1)*delx*e01 + (j-1)*dely*e02 + (k-1)*delz*e03
         if (xi(1)>g_xcord(g_nx)) xi(1)=xi(1)-1
	 if (xi(1)<g_xcord(1))    xi(1)=xi(1)+1
	 if (xi(2)>g_ycord(g_ny)) xi(2)=xi(2)-1
	 if (xi(2)<g_ycord(1))    xi(2)=xi(2)+1
	 if (xi(3)>g_zcord(g_nz)) xi(3)=xi(3)-1
	 if (xi(3)<g_zcord(1))    xi(3)=xi(3)+1
         if (irot==1) then
             call interp3d(data(i,j,k-k1+1),field,xi)
         else
            ! no interpolation - just use closes gridpoint to xi
            ii=nint((xi(1)-xcord(nx1))/delx);  
	    jj=nint((xi(2)-ycord(ny1))/dely);  
	    kk=nint((xi(3)-zcord(nz1))/delz);  

#ifdef CHECK_IERR
            err= (ii-(xi(1)-xcord(nx1))/delx)**2
            err = err + (jj-(xi(2)-ycord(ny1))/dely)**2  
	    err = err + (kk-(xi(3)-zcord(nz1))/delz)**2  
#endif

            ii=ii+nx1
            jj=jj+ny1
            kk=kk+nz1

	    if (kk>=nz1 .and. kk<=nz2 .and. &
                ii>=nx1 .and. ii<=nx2 .and. &
                jj>=ny1 .and. jj<=ny2) then
	        ! interpolate 
	        data(i,j,k-k1+1)=field(ii,jj,kk)

#ifdef CHECK_IERR
                if (err>.001**2) then
                   print *,'warning: interpolation error: '
                   print *,'ii',ii-nx1,(xi(1)-xcord(nx1))/delx
                   print *,'jj',jj-ny1,(xi(2)-ycord(ny1))/dely;  
                   print *,'kk',kk-nz1,(xi(3)-zcord(nz1))/delz;  
                endif
#endif
             else	
	        data(i,j,k-k1+1)=-9d99
            endif
         endif
         
      enddo
   enddo
enddo
#ifdef USE_MPI
allocate(data2(ssize,ssize,k2-k1+1))
data2=data
call mpi_reduce(data2,data,ssize*ssize*(k2-k1+1),MPI_REAL8,MPI_MAX,io_pe,comm_3d,ierr)
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







subroutine extract_subcube(nsx,nsy,nsz,dx,dy,dz,data,field)
!
! extract the gridpoints dx(1:ns),dy(1:ns),dz(1:ns) from the global data
!
! NOTE: only io_pe will have the correct, interpolated data
!
integer :: nsx,nsy,nsz
real*8 :: dx(nsx),dy(nsy),dz(nsz)
real*8 :: data(nsx,nsy,nsz)
real*8 :: field(nx,ny,nz)
integer :: irot
real*8 :: rmatrix(3,3,3)

real*8 :: xi,yi,zi,pos(3)
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

         if (n1(i)==0 .or. n2(j)==0 .or. n3(k)==0) then
            data(i,j,k)=-9d200
         else
            ! we have this data, copy;
            data(i,j,k)=field(n1(i),n2(j),n3(k))
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






subroutine interp3d(field_interp,field,pos)
! data:  field(i,j,k)
! interpolation point:  pos(1:3)
! output:  field_interp

! input
real*8 :: field(nx,ny,nz)
real*8 :: pos(3),field_interp

! local
real*8 :: Q2d(4,4)
real*8 :: Q1d(4)
integer :: jj,kk,igrid,jgrid,kgrid,jc,kc
real*8 :: xc,yc,zc


! find position in global grid:
igrid = 1 + floor( (pos(1)-g_xcord(1))/delx )
jgrid = 1 + floor( (pos(2)-g_ycord(1))/dely )
kgrid = 1 + floor( (pos(3)-g_zcord(1))/delz )

! compute a new point in the center of the above cell:
! (do this to avoid problems with 2 cpus both claiming a point
! on the boundary of a cell)
xc=.5*(g_xcord(igrid)+g_xcord(igrid+1))
yc=.5*(g_ycord(jgrid)+g_ycord(jgrid+1))
zc=.5*(g_zcord(kgrid)+g_zcord(kgrid+1))

! find cpu which owns this cell (xi,yi,zi):
igrid=(xc-xcord(nx1))/delx
jgrid=(yc-ycord(ny1))/dely
kgrid=(zc-zcord(nz1))/delz

! find cpu which owns the grid point (xc,yc)
if (xcord(nx1)<xc .and. xc<xcord(nx2)+delx .and. &
     ycord(ny1)<yc .and. yc<ycord(ny2)+dely .and. &
     zcord(nz1)<zc .and. zc<zcord(nz2)+delz ) then
   

   ! 3d --> 2d y-z planes --> 1d z-lines --> point
   do kk=1,4
      do jj=1,4
         ! interpolate xcord(igrid-1:igrid+2) to xcord=particle(i,1)
         ! onto 2-d planes in y-z
         ! data  field(igrid-1:igrid+2, jgrid-2+jj) 
         xc = 1 + (pos(1)-xcord(igrid))/delx
         jc = jgrid-2+jj
         kc = kgrid-2+kk
         call interp4(field(igrid-1,jc,kc),field(igrid,jc,kc),&
              field(igrid+1,jc,kc),field(igrid+2,jc,kc),&
              xc,Q2d(jj,kk))
      enddo
   enddo
   
   do kk=1,4
      ! interpolate ycord(jgrid-1:jgrid+2) to ycord=particle(j,1)
      ! onto lines in z
      ! data  field(igrid-1:igrid+2, jgrid-2+jj) 
      yc = 1 + (pos(2)-ycord(jgrid))/dely
      call interp4(Q2d(1,kk),Q2d(2,kk),Q2d(3,kk),Q2d(4,kk),yc,Q1d(kk))
   enddo
   
   ! interpolate zcord(kgrid-1:kgrid+2) to zcord=particle(k,2)
   ! data:  Q1d(1:4,j)
   zc = 1 + (pos(3)-zcord(kgrid))/delz
   call interp4(Q1d(1),Q1d(2),Q1d(3),Q1d(4),zc,field_interp)
else
   field_interp=-9d100
endif

end subroutine interp3d


subroutine rotcube(ir,jr,x0,x1,x2,x3,c0)
!
! Input:  4 points x0(), x1(), x2(), x3()
!         coordinates of center: c0()
!
!  rotate (about center) by 90 degrees so as to swap the ir and jr axis
! 
!   
integer :: ir,jr
real*8 :: x0(3),x1(3),x2(3),x3(3),c0(3)

real*8 :: xshift(3),tmp

xshift=x0-c0
tmp=xshift(ir); xshift(ir)=xshift(jr); xshift(jr)=-tmp
x0=xshift+c0

xshift=x1-c0
tmp=xshift(ir); xshift(ir)=xshift(jr); xshift(jr)=-tmp
x1=xshift+c0

xshift=x2-c0
tmp=xshift(ir); xshift(ir)=xshift(jr); xshift(jr)=-tmp
x2=xshift+c0

xshift=x3-c0
tmp=xshift(ir); xshift(ir)=xshift(jr); xshift(jr)=-tmp
x3=xshift+c0


end subroutine



   
end module 
