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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Fully resolved compressible Navier Stokes code
! copyright 2001 by Hal Marshal, Mark Taylor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program conv
use params
use transpose
implicit none
real*8,save :: Q(nx,ny,nz,n_var)
real*8,save :: p(nx,ny,nz)
real*8,save :: work1(nx,ny,nz)
real*8,save :: work2(nx,ny,nz)
CPOINTER :: fid
real*8 :: time
character(len=80) message,fname
integer ierr,i,j,k,n


call init_mpi       
call init_mpi_comm3d()
call init_model      

fname="temp0001.0000.u"
call singlefile_io(time,Q(1,1,1,1),fname,work1,work2,1,io_pe)
fname="temp0001.0000.v"
call singlefile_io(time,Q(1,1,1,2),fname,work1,work2,1,io_pe)
fname="temp0001.0000.w"
call singlefile_io(time,Q(1,1,1,3),fname,work1,work2,1,io_pe)


if (my_pe==io_pe) then
fname="iogrid"
call copen(fname,"w",fid,ierr)
call cwrite8(fid,o_nx,1)
call cwrite8(fid,o_ny,1)
call cwrite8(fid,o_nz,1)
endif

do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   p(i,j,k)=xcord(i)
enddo
enddo
enddo

call output1(p,work1,work2,fid,io_pe)

do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   p(i,j,k)=ycord(i)
enddo
enddo
enddo

call output1(p,work1,work2,fid,io_pe)

do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   p(i,j,k)=zcord(i)
enddo
enddo
enddo

call output1(p,work1,work2,fid,io_pe)

p=1

if (my_pe==io_pe) then
fname="ioflow"
call copen(fname,"w",fid,ierr)
call cwrite8(fid,o_nx,1)
call cwrite8(fid,o_ny,1)
call cwrite8(fid,o_nz,1)
call cwrite8(fid,p,1)
call cwrite8(fid,p,1)
call cwrite8(fid,p,1)
endif

!call output1(p,work1,work2,fid,io_pe)
call output1(Q(1,1,1,1),work1,work2,fid,io_pe)
call output1(Q(1,1,1,2),work1,work2,fid,io_pe)
call output1(Q(1,1,1,3),work1,work2,fid,io_pe)
call output1(p,work1,work2,fid,io_pe)
	


call close_mpi


end program conv

