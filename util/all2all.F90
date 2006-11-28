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
real*8 :: tmx3,tmx4,tmx2,tmx1,tmx_max,tmx_min
real*8 :: a2a_max,a2a_min

call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD,my_rank,ierr)
call mpi_comm_size(MPI_COMM_WORLD,nproc2,ierr)
if (nproc2 /= nproc) then
   print *,'error: nproc,nproc2=',nproc,nproc2
   call mpi_abort(MPI_COMM_WORLD,1,ierr)
endif



! fill with some data
u=my_rank

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

tmx2 = MPI_Wtime()
call MPI_Alltoall(spack,scount,MPI_DOUBLE_PRECISION,rpack,scount,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
! we received the pencil of data u(2n:2n+1,1:ny,2) from process n.  unpack:
tmx3 = MPI_Wtime()

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

tmx4 = MPI_Wtime()

tmx4=tmx4-tmx1
call MPI_Reduce(tmx4,tmx_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(tmx4,tmx_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)

tmx4=tmx3-tmx2
call MPI_Reduce(tmx4,a2a_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
call MPI_Reduce(tmx4,a2a_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)


if(my_rank==0) then
   print *,'Min alltoall : (in minutes)',a2a_min/60
   print *,'Max alltoall : (in minutes)',a2a_max/60
   print *,'Min alltoall + copy time: (in minutes)',tmx_min/60
   print *,'Max alltoall + copy time: (in minutes)',tmx_max/60
endif

! test data:
count=0
do i=1,nz_proc
do j=1,ny
do k=1,nz
   if (uz(i,j,k) /= (k-1)/2 ) then
      count=count+1
      if (count<100) print *,'Error: ',i,j,k,uz(i,j,k)
   endif
enddo
enddo
enddo
call MPI_Finalize(ierr)
end program





