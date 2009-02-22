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

subroutine init_mpi
use params
use mpi

implicit none
integer ierr1,ierr2,ierr3,rank
character(len=80) message

my_pe=0
mpicoords=0
mpidims(1)=ncpu_x
mpidims(2)=ncpu_y
mpidims(3)=ncpu_z
io_pe=0
initial_live_procs=1


my_x=0
my_y=0
my_z=0


#ifdef USE_MPI

call mpi_init(ierr1)
if (ierr1/=0) call abortdns("mpi_init failure")
call mpi_comm_rank(MPI_COMM_WORLD,my_world_pe,ierr2)
if (ierr2/=0) call abortdns("mpi_comm_rank failure")
call mpi_comm_size(MPI_COMM_WORLD,initial_live_procs,ierr3)
if (ierr3/=0) call abortdns("mpi_comm_size failure")

#else
! make sure we are running with MPI!
if (ncpu_x * ncpu_y * ncpu_z /= 1) then
   call abortdns("code was not compiled with MPI support!")
endif
#endif


ncpus=initial_live_procs

end subroutine





subroutine init_mpi_comm3d
use params
use mpi

implicit none
logical isperiodic(3),reorder
integer ierr1,ierr2,ierr3,rank,key,color,my_pe2
character(len=80) message

if (ncpu_x * ncpu_y * ncpu_z /= ncpus) then
   call print_message("Error: incorrect number of cpus");
   write(message,'(a,i5,a,i3,a,i3,a,i3)') "Parallel decomposition requested: ncpus= ", &
      ncpu_x*ncpu_y*ncpu_z," = ",ncpu_x," x",ncpu_y," x",ncpu_z
   call print_message(message)
   write(message,'(a,i5)') "Total cpus = ",initial_live_procs
   write(message,'(a,i5)') "cpus for cartesian communicator = ",ncpus
   call print_message(message)
   call abortdns("Terminating.")
endif

#ifdef USE_P3DFFT
if (ncpu_x==1) p3dfft_cartmap=.true.
#endif


#ifdef USE_MPI
isperiodic(1)=.false.
isperiodic(2)=.false.
isperiodic(3)=.false.

if (p3dfft_cartmap) then
   ! to match P3DFFT cartesian communicator, we need to call 
   ! mpi_cart_create exactly as P3DFFT does:
   reorder=.false.
   mpidims(1)=ncpu_z
   mpidims(2)=ncpu_y
   call mpi_cart_create(MPI_COMM_WORLD,2,mpidims,isperiodic,reorder,comm_3d,ierr1)
else
   reorder=.true.
   call mpi_cart_create(MPI_COMM_WORLD,3,mpidims,isperiodic,reorder,comm_3d,ierr1)
endif
call mpi_comm_rank(comm_3d,my_pe,ierr2)

if (ierr1/=0) call abortdns("mpi_cart_create failure")

call cart_get(comm_3d,3,mpidims,isperiodic,mpicoords,ierr1)
!print *,'me=',my_pe," mpi coords: ",mpicoords(1),mpicoords(2),mpicoords(3)

! get my_pe (sanity check?)
call cart_rank(comm_3d,mpicoords,my_pe2,ierr2)
if (my_pe /= my_pe2) then
   call abortdns("Error: mpicoords -> cart_rank mapping error")
endif

! create a communicator with just io_pe in it:
color = 0
if (my_pe==io_pe) color=1
key=0

! everyone with ntot>0 joins a new group, comm_sforcing
call mpi_comm_split(comm_3d,color,key,comm_1,ierr2);
if (color==0) call mpi_comm_free(comm_1,ierr2)




write(message,'(a,i5,a,i4,a,i4,a,i4)') "Running parallel.  NCPUS=", &
   ncpus," = ",ncpu_x," x",ncpu_y," x",ncpu_z
call print_message(message)
#else
call print_message("Running single threaded")
#endif



end subroutine




subroutine close_mpi
use params
use mpi
integer ierr

#ifdef USE_MPI
call mpi_finalize(ierr)
#endif

end subroutine


#ifdef USE_MPI
subroutine cart_get(comm1,size,dims,isperiodic,coords,ierr1)
use params
implicit none
integer :: size
logical :: isperiodic(size)
integer :: dims(size),coords(size)
integer :: comm1,ierr1,itmp
integer :: coords_map(2)

if (p3dfft_cartmap) then
   call mpi_cart_get(comm1,2,dims,isperiodic,coords_map,ierr1)
   coords(1)=0
   coords(2)=coords_map(2)
   coords(3)=coords_map(1)
else
   call mpi_cart_get(comm1,3,dims,isperiodic,coords,ierr1)
endif
if (ierr1/=0) call abortdns("mpi_cart_get failure")


end subroutine


subroutine cart_rank(comm1,coords,pe1,ierr1)
use params
implicit none
integer :: coords(*)
integer :: comm1,ierr1,pe1
integer :: coords_map(2)

if (p3dfft_cartmap) then
   coords_map(1)=coords(3)
   coords_map(2)=coords(2)
   call mpi_cart_rank(comm1,coords_map,pe1,ierr1)
else
   call mpi_cart_rank(comm1,coords,pe1,ierr1)
endif
if (ierr1/=0) call abortdns("mpi_cart_rank failure")

end subroutine
#endif

