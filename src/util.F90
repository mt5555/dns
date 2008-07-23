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

subroutine f90random_seed(istart)
integer :: istart,k
integer,allocatable :: seed(:)

call random_seed(size=k)
allocate(seed(k))
call random_seed(get=seed)
seed=seed+istart
call random_seed(put=seed)
deallocate(seed)

end subroutine




subroutine abortdns(message)
use mpi
use params
implicit none
integer ierr
character(len=*) message
character(len=15) :: pre="ABORTDNS: "
write(*,'(a)') pre // message

call flush(6)

#ifdef USE_MPI
   call mpi_abort(comm_3d,1,ierr)
#endif
stop
end subroutine


subroutine print_message(message)
use params
use mpi
implicit none
character(len=*) message
integer :: ierr
if (my_pe==io_pe) then
   write(*,'(a)') trim(message)
endif

!done do this: not everyone calls print_message:
! for parallel debugging with print_message():
!call mpi_barrier(comm_3d,ierr)

end subroutine


subroutine wallclock(tmx)

#ifdef USE_MPI
use params
use mpi
implicit none
real*8 tmx
#ifdef MPI_UNDERSCORE
real*8 mpi_wtime
#endif
tmx = mpi_wtime()

#else
use params
implicit none
integer count,count_rate,count_max
real*8 tmx
call system_clock(count,count_rate,count_max)
tmx=real(count,r8kind)/real(count_rate,r8kind)

#endif

end subroutine





subroutine logplotASCII(spectrum,n,title)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ASCII plot of spectrum(0:n) on a log-log scale
!
! Y-axis scale is  10^-10 ... 10^0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use params
implicit none
integer :: n
real*8  ::  spectrum(0:n)
character(len=*) :: title


! local variables
integer,parameter :: numx=20,numy=13
integer i,j
character :: plot(0:numx,0:numy)
real*8 cspec(0:numx)
real*8 :: delta,ireal
integer ix,iy


if (my_pe==io_pe) then
cspec=0
plot=" "

#if 0
! bin all values
do i=0,n
   if (i==0) then
      ix=0
   else
      ix=1 + ( (numx-1)*log10(real(i))/log10(real(n)) )     ! ranges from 1..numx
   endif
   if (ix<0) ix=0
   if (ix>numx) ix=numx
   cspec(ix)=cspec(ix)+spectrum(i)
enddo
#endif

! interpolate
cspec(0)=spectrum(0)
do ix=1,numx
   ireal = 10**( log10(real(n)) * (ix-1)/real(numx-1) )
   i=floor(ireal)
   delta=ireal-i
   if (i>=0 .and. i<=n-1) then
      cspec(ix)=spectrum(i)*(1-delta)+spectrum(i+1)*delta
   endif
enddo


! scale from 0..numy, log scale
do i=0,numx
   iy = -(numy/10.0) * log10(1d-200+cspec(i))  ! ranges from 0..numy
   if (iy<0) iy=0
   if (iy>numy) iy=numy
   cspec(i)=iy
enddo

do i=0,numx
    j=cspec(i)
    plot(i,j)="*"
enddo

print *
print *,"1E0   |",plot(:,0),title
do i=1,numy-1
   print *,"      |",plot(:,i)
enddo
print *,"1E-10 |",plot(:,numy)
print *,"      +---------------------+"
write (*,'(a,i4)') "      k=0                 k=",n

endif
end subroutine




subroutine plotASCII(spectrum,n,title,ymin,ymax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ASCII plot of spectrum(0:n) on a log-log scale
!
! Y-axis scale is  10^-10 ... 10^0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use params
implicit none
integer :: n
real*8  ::  spectrum(0:n),ymin,ymax
character(len=*) :: title


! local variables
integer,parameter :: numx=20,numy=13
integer i,j
character :: plot(0:numx,0:numy)
real*8 cspec(0:numx)
real*8 :: delta,ireal
integer ix,iy


if (my_pe==io_pe) then
cspec=0
plot=" "

! interpolate
cspec(0)=spectrum(0)
do ix=1,numx
   ireal = 10**( log10(real(n)) * (ix-1)/real(numx-1) )
   !ireal = n * (ix-1)/real(numx-1) 
   i=floor(ireal)
   delta=ireal-i
   if (i>=0 .and. i<=n-1) then
      cspec(ix)=spectrum(i)*(1-delta)+spectrum(i+1)*delta
   endif
enddo


! scale from 0..numy, log scale
do i=0,numx
   iy = numy*(cspec(i) - ymax)/(ymin-ymax)
   if (iy<0) iy=0
   if (iy>numy) iy=numy
   cspec(i)=iy
enddo

do i=0,numx
    j=cspec(i)
    plot(i,j)="*"
enddo

print *
write(*,'(1x,f6.3,99a)') ymax,"|",plot(:,0),title
do i=1,numy-1
   print *,"      |",plot(:,i)
enddo
!print *,"1E-10 |",plot(:,numy)
write(*,'(1x,f6.3,99a)') ymin,"|",plot(:,numy)
print *,"      +---------------------+"
write (*,'(a,i4)') "      k=0                 k=",n

endif
end subroutine



subroutine random_data(buf,n)
implicit none
integer n
real*8 :: buf(n)
!call random_number(buf)
call gaussian(buf,n)
end subroutine




!.... Gaussian Random number generator ran1 from Numerical recipes
subroutine gaussian(buf,n)
implicit none
      
integer n
real*8 :: buf(n)
real*8 ::  fac,rsq,ran1(2)
integer i,j1,j2

! n=10:  i=1..5   j1max=9, j2max=10 
! n=11:  i=1..6   j1max=11, j2max=12(ignored)
do i=1,(n+1)/2
   j1=2*i-1
   j2=j1+1
   do 
      call random_number(ran1)
      ran1=2*ran1-1
      rsq = ran1(1)**2 + ran1(2)**2
      if (rsq<1 .and. rsq>0) exit
   enddo
   fac = Sqrt(-2.0 * Log(rsq)/rsq)
   buf(j1) = ran1(1) * fac
   if (j2<=n) buf(j2)= ran1(2) * fac
enddo

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




integer function zerosign(i)
integer i
if (i==0) then
   zerosign=0
else if (i<0) then
   zerosign=-1
else
   zerosign=1
endif
end function
