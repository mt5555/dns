program slab_alltoall_test
implicit none
include "mpif.h"

integer,parameter :: nx=2048,ny=2048,nz=2048
integer,parameter :: nproc=1024
integer,parameter :: nz_proc=nz/nproc

integer :: ierr,my_rank,nproc2,i,k,j,n,count,scount
real*8 :: u(nx,ny,nz_proc)
real*8 :: spack(nx*ny*2)
real*8 :: rpack(nx*ny*2)
real*8 :: uz(nz_proc,ny,nz)    ! the transposed array
real*8 :: tmx2,tmx1,tmx_max,tmx_min

call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD,my_rank,ierr)
call mpi_comm_size(MPI_COMM_WORLD,nproc2,ierr)
if (nproc2 /= nproc) then
   print *,'error: nproc,nproc2=',nproc,nproc2
   call mpi_abort(MPI_COMM_WORLD,1,ierr)
endif




! fill with some data
u=1

tmx1 = MPI_Wtime()
! send u(2n:2n+1,1:ny,1:2) to process n   
count=0
do n=1,nproc
   do i=1,nz_proc
      do j=1,ny
         do k=1,nz_proc
            count=count+1
            spack(count)=u(i,j,k)
         enddo
      enddo
   enddo
enddo
scount=count/nproc

call MPI_Alltoall(spack,scount,MPI_DOUBLE_PRECISION,rpack,scount,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
! we received the pencil of data u(2n:2n+1,1:ny,2) from process n.  unpack:

count=0
do n=1,nproc
   do i=1,nz_proc
      do j=1,ny
         do k=1,nz_proc
            count=count+1
            ! rpack(count) = u(i,j,k) from process n.  process n 
            ! has the n'th slab of data   
            uz(i,j,nz_proc*(n-1)+k)=rpack(count)
         enddo
      enddo
   enddo
enddo

tmx2 = MPI_Wtime()
tmx2=tmx2-tmx1
call MPI_Reduce(tmx2,tmx_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(tmx2,tmx_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)

if(my_rank==0) then
   print *,'Min alltoall + copy time: (in minutes)',tmx_max/60
   print *,'Max alltoall + copy time: (in minutes)',tmx_max/60
endif


end program





