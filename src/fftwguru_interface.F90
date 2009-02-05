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
#if 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! our wrapper for ECMWF FFT99 
! 
!
! provides public interfaces:
!   fft_interface_init               call this before using any other routines
!   fft1                             fft along first dimension of 3D array
!   ifft1	                     ifft along first dimension of 3D array
! 
! Routines work on data of the form:  p(n1d,n2d,n3d)
! Size of the grid point data         p(1:n1,1:n2,1:n3)
! Size of fourier coefficients        p(1:n1+2,1:n2+2,1:n3+2)
!

FFT data representation:

sum over m=1..n/2:

   f = fhat(1)  +  2 fhat(2*m+1) cos(m*2pi*x) - 2*fhat(2m+2) sin(m*2pi*x)


     if isign = +1, and m coefficient vectors are supplied
     each containing the sequence:

     a(0),b(0),a(1),b(1),...,a(n/2),b(n/2)  (n+2 values)

     then the result consists of m data vectors each
     containing the corresponding n+2 gridpoint values:

     x(0), x(1), x(2),...,x(n-1),0,0.

     note: the fact that the gridpoint values x(j) are real
     implies that b(0)=b(n/2)=0.  for a call with isign=+1,
     it is not actually necessary to supply these zeros.


In otherwords:

 grid space data:    1 2 3 4 5 6 7 8 * *
 
 fourier space:      0 0 1 1 2 2 3 3 4 4
 rearranged:         0 4 1 1 2 2 3 3 * *


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif


module fft_interface
use params
use transpose
implicit none
include 'fftw3.f'

integer :: init=0
integer, parameter ::  num_fftsizes=8
type fftdata_d
   CPOINTER :: plan
   integer :: size
   integer :: howmany
   integer :: direction
   integer :: stride
end type
type(fftdata_d) :: fftdata(num_fftsizes)


contains 





subroutine fft_interface_init(f,work)
real*8,optional :: f(nx,ny,nz)
real*8,optional :: work(nx,ny,nz)
integer :: i,n1,n1d,n2,n2d,n3,n3d,index

init=1
do i=1,num_fftsizes
   fftdata(i)%size = 0	
   fftdata(i)%howmany = 0	
   fftdata(i)%stride = 0
   fftdata(i)%direction = 0	
enddo


if (present(f)) then

   call print_message("Initializing FFTW Plans...")
   ! force initialization of plans.
   f=0
   call transpose_to_z(f,work,n1,n1d,n2,n2d,n3,n3d)
   call getindex(n1,n2*n3,1,1,index,work)
   call getindex(n1,n2*n3,1,-1,index,work)
   
   call transpose_to_y(f,work,n1,n1d,n2,n2d,n3,n3d)
   call getindex(n1,n2*n3,1,1,index,work)
   call getindex(n1,n2*n3,1,-1,index,work)
   
   call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d)
   call getindex(n1,n2*n3,1,1,index,work)
   call getindex(n1,n2*n3,1,-1,index,work)

   call getindex(ny2,nx2,nx,-1,index,work)
   ! for howmany = 1 debug case:
   ! call getindex(ny2,1,nx,-1,index,work)
   call print_message("Done initializing FFTW Plans.")
endif


end subroutine





subroutine fftinit(n,howmany,stride,direction,index,f,n1d)
integer n,index,direction,howmany,n1d,stride
real*8 :: f(n1d,howmany)
character(len=80) message

if (init==0) call abortdns("fftw_interface.F90: call fft_interface_init to initialize first!")
if (n>1000000) call abortdns("fftw_interface.F90: n>1 million")

fftdata(index)%size=n
fftdata(index)%howmany=howmany
fftdata(index)%stride=stride
fftdata(index)%direction=direction

call print_message('Starting FFTW_MEASURE:')
write(message,'(a,i6,i3,i12,i3)') 'size, direction, howmany, stride=',n,direction,howmany,stride
call print_message(message)


fftdata(index)%plan=1

if (stride==1) then
   if (direction==1) then
      call dfftw_plan_guru_dft_r2c(fftdata(index)%plan,1,n,1,1, &
           1, howmany,n1d,n1d/2,f,f,FFTW_MEASURE)
   endif
   if (direction==-1) then
      call dfftw_plan_guru_dft_c2r(fftdata(index)%plan,1,n,1,1, &
           1, howmany,n1d/2,n1d,f,f,FFTW_MEASURE)
   endif
else
   if (direction==1) then
      call dfftw_plan_guru_split_dft_r2c(fftdata(index)%plan,1,&
           n,stride,stride*2,   1,    howmany, 1,1,&
           f,f(1,1),f(1,2),FFTW_MEASURE)
   endif
   if (direction==-1) then
      call dfftw_plan_guru_split_dft_c2r(fftdata(index)%plan,1,n,stride*2,stride, &
           1, howmany,1,1,f(1,1),f(1,2),f,FFTW_MEASURE)
   endif
endif

if (fftdata(index)%plan==0) then
   call abortdns("fftwguru_interface.F90:  fftw_plan_guru failed")
endif

end subroutine




subroutine getindex(n1,howmany,stride,direction,i,f)
!
!  call this routine with optional last argument to compute plans
!  with FFT_MEASURE
!
!  if plans have been computed, call this routine to retrieve index
!  of the prevously computed plan
!
!  if no plans have been computed, it will return i=0, and the FFTs
!  will use the FFTW advanced interface (which computes a new plan each time)
!
!
integer :: n1,direction,howmany,stride
real*8,optional :: f(*)
character(len=180) message_str
integer i,k


i=0
do 
   i=i+1
   if (i>num_fftsizes) then
      write(message_str,'(a,i10)') "fftwguru_interface.F90:  Failed initializing an fft of size =",n1
      call abortdns(message_str)
   endif

   if (fftdata(i)%size==0) then  ! we didnt find a plan - make a new one
      if (present(f)) then
         call fftinit(n1,howmany,stride,direction,i,f,n1+2)      
      else
         if (init.eq.1) then
            write(message_str,'(a,i10)') "NOTE: FFTW Guru interface not initialized.  Reverting for n=",n1
            call print_message(message_str)
         endif
         init=init+1
         i=0  ! flag 
      endif
      exit 
   endif
   if (n1==fftdata(i)%size .and. direction==fftdata(i)%direction &
          .and. howmany==fftdata(i)%howmany &
          .and. stride==fftdata(i)%stride .and. 0/=fftdata(i)%plan) exit 
enddo
end subroutine





subroutine fft_get_mcord(mcord,n)
!
!  i=1   0 cosine mode             mcord=0
!  i=2   n/2 cosine mode           mcord=  n/2
!  i=3   1 cosine mode             mcord =  1
!  i=4   1 sine mode               mcord = -1
!  i=5   2 cosine mode             mcord =  2
!  i=6   2 sine mode               mcord=  -2
!  etc...
!
integer n,mcord(:)
integer i
do i=1,n
   mcord(i)=(i-1)/2	
   if (mod(i,2)==0) mcord(i)=-mcord(i)
   if (i==2) mcord(i)=n/2
enddo
end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute 3D in-place iFFT of p
! FFT taken along first direction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ifft1(p,n1,n1d,n2,n2d,n3,n3d)
integer n1,n1d,n2,n2d,n3,n3d
real*8 p(n1d,n2d,n3d)

real*8 :: tmx1,tmx2
character(len=80) message_str
integer j,k,index
CPOINTER :: plan

call wallclock(tmx1)
if (tims(19)==0) ncalls(19)=0  ! timer was reset, so reset counter too

if (n1==1) return
ASSERT("ifft1: dimension too small ",n1+2<=n1d)
ASSERT("ifft1: dimension too small ",n2==n2d)
call getindex(n1,n2*n3,1,-1,index)


! number of FFTs to compute:  n2*n3
! stride:  n1d

do k=1,n3
   do j=1,n2
      ! move the last cosine mode back into correct location:
      p(n1+1,j,k)=p(2,j,k)
      p(2,j,k)=0             ! not needed?
      p(n1+2,j,k)=0          ! not needed?
   enddo
enddo

if (index==0) then
   call dfftw_plan_many_dft_c2r(plan, 1, n1,n2*n3,&
          p, n1d/2, 1, n1d/2, &
          p, n1d, 1, n1d, FFTW_ESTIMATE)
   call dfftw_execute(plan)
   call dfftw_destroy_plan(plan)
else
   call dfftw_execute_dft_c2r(fftdata(index)%plan,p,p)
endif


call wallclock(tmx2) 
tims(19)=tims(19)+(tmx2-tmx1)          
ncalls(19)=ncalls(19)+1
end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute 3D in-place iFFT of p
! FFT taken along second direction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ifft1_dim2(p,n1,n1d,n2,n2d,n3,n3d)
integer n1,n1d,n2,n2d,n3,n3d
real*8 p(n1d,n2d,n3d)

real*8 :: tmx1,tmx2
character(len=80) message_str
integer i,j,k,index
CPOINTER :: plan

call wallclock(tmx1)
if (tims(19)==0) ncalls(19)=0  ! timer was reset, so reset counter too

if (n2==1) return
ASSERT("ifft1: dimension too small ",n2+2<=n2d)
! do we need this???
ASSERT("ifft1: dimension too small ",n1==n1d)
call getindex(n2,n1,n1d,-1,index)

!used when doing howmany=1
!call getindex(n2,1,n1d,-1,index)


if (index==0) then
   stop 'error: cant do this with advanced interface'
   ! because of how are array is stored, we have to call FFTW n3 times
endif

! number of FFTs to compute:  n2*n3
! stride:  n1d

do k=1,n3
   do i=1,n1
      ! move the last cosine mode back into correct location:
      p(i,n2+1,k)=p(i,2,k)
      p(i,2,k)=0             ! not needed?
      P(i,n2+2,k)=0          ! not needed?
!      call dfftw_execute_split_dft_c2r(fftdata(index)%plan,&
!           p(i,1,k),p(i,2,k),p(i,1,k))
   enddo
   call dfftw_execute_split_dft_c2r(fftdata(index)%plan,p(1,1,k),p(1,2,k),p(1,1,k))
enddo



call wallclock(tmx2) 
tims(19)=tims(19)+(tmx2-tmx1)          
ncalls(19)=ncalls(19)+1
end subroutine






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute 3D in-place FFT of p
! FFT taken along first direction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft1(p,n1,n1d,n2,n2d,n3,n3d)
integer n1,n1d,n2,n2d,n3,n3d
real*8 :: p(n1d,n2d,n3d)
real*8 :: tmx1,tmx2
integer j,k,index
CPOINTER :: plan

call wallclock(tmx1)
if (tims(18)==0) ncalls(18)=0  ! timer was reset, so reset counter too


if (n1==1) return
ASSERT("fft1: dimension too small ",n1+2<=n1d)
ASSERT("fft1: dimension too small ",n2==n2d)
call getindex(n1,n2*n3,1,1,index)

! number of FFTs to compute:  n2*n3
! stride:  n1d
if (index==0) then
   call dfftw_plan_many_dft_r2c(plan, 1, n1,n2*n3,&
      p, n1d, 1, n1d, p, n1d/2, 1, n1d/2,  FFTW_ESTIMATE)
   call dfftw_execute(plan)
   call dfftw_destroy_plan(plan)
else
   call dfftw_execute_dft_r2c(fftdata(index)%plan,p,p)
endif

do k=1,n3
   !   do j=1,n2
   !         p(n1+1,jj,k)=0
   !         p(n1+2,jj,k)=0
   !   enddo
   !     move the last cosine mode into slot of first sine mode:
   do j=1,n2
      p(2,j,k)=p(n1+1,j,k)
      p(1:n1,j,k)=p(1:n1,j,k)/n1 
   enddo
enddo
call wallclock(tmx2) 
tims(18)=tims(18)+(tmx2-tmx1)          
ncalls(18)=ncalls(18)+1
end subroutine




end ! module mod_fft_interface


