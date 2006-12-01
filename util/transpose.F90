#if 0
compile with:

ftn -target=catamount -byteswapio -fast -r8 -Mpreprocess -DMESHSIZE=512
-DPROCY=4 -DPROCZ=4 -o ./transpose.x transpose.f90

on Redstorm:
 
4096^3 on 8192 cores:  (I dont have a run like this)
mpif90  -fast -r8 -Mpreprocess -DMESHSIZE=4096 \
-DPROCY=2 -DPROCZ=4096 -o ./transpose.x transpose.F90
yod -VN -sz 8192  transpose.x

4096^3 on 16384 cores: 

mpif90  -fast -r8 -Mpreprocess -DMESHSIZE=4096 \
-DPROCY=128 -DPROCZ=128 -o ./t2.x transpose.F90


yod -VN -sz 16384  transpose.x 


 



"bandwidth"   8*MESHSIZE^3 / NPROC / time_in_seconds


#endif

program test_transpose
implicit none

! This code is for pencils, in y and z directions. 
#define PROCX 1 

!Problem size and processor counts. To be #defined
!No divisibility tests performed here
integer, parameter :: nxg=MESHSIZE, nyg=MESHSIZE, nzg=MESHSIZE
integer, parameter :: npx = PROCX, npy=PROCY, npz=PROCZ
integer, parameter :: nproc = npx*npy*npz
integer, parameter :: nx = nxg/npx, ny = nyg/npy, nz=nzg/npz

include 'mpif.h'
integer my_rank
integer ycomm, zcomm
integer yid, zid

integer ierr, nproc2

real *8 :: ux(nxg, ny, nz)
real *8 :: uy(nxg/PROCY, nyg, nz)

real*8 :: spack(nxg*ny*nz), rpack(nxg*ny*nz)

integer bsize
integer i, j, k, n

real*8 :: tmx3,tmx4,tmx2,tmx1, tme
real*8 :: a2a_max=0.0, a2a_min=1e10, a2a_sum=0.0
real*8 :: trs_max=0.0, trs_min=1e10, trs_sum=0.0
real*8 :: tme_g
integer :: iter, iter_max=1000

real*8 :: error, error_max, error_max_g

call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD,my_rank,ierr)
call mpi_comm_size(MPI_COMM_WORLD,nproc2,ierr)
if (nproc2 /= nproc) then
   if(my_rank .eq. 0) print *,'error: nproc,nproc2=',nproc,nproc2
   call mpi_abort(MPI_COMM_WORLD,1,ierr)
endif

!Form y and z communicators
zid = my_rank/npy
yid = mod(my_rank, npy)

call mpi_comm_split(mpi_comm_world, zid, my_rank, ycomm, ierr)
call mpi_comm_split(mpi_comm_world, yid, my_rank, zcomm, ierr)

!Fill with data that will allow me to debug the transpose operation
do k = 1, nz
  do j = 1, ny
    do i = 1, nx
      ux(i,j,k) = 1e8*real(i)+1e4*real(yid*ny+j)+1e0*real(zid*nz+k)
    end do
  end do
end do

iter: do iter = 1, iter_max

  tmx1 = MPI_Wtime()
  
  !Pack data for alltoall in y-direction
  bsize = nxg*ny*nz/npy
  do n = 0, npy-1
    spack(n*bsize+1:(n+1)*bsize) = & 
      reshape(ux(n*nxg/npy+1:(n+1)*nxg/npy, 1:ny, 1:nz), (/bsize/))
  end do
  
  tmx2 = MPI_Wtime()
  
  call MPI_Alltoall(spack,bsize,MPI_REAL8,rpack,bsize,MPI_REAL8,ycomm,ierr)
  
  tmx3 = MPI_Wtime()
  
  !Unpack the data
  do n = 0, npy-1
    uy(1:nxg/PROCY, n*ny+1:(n+1)*ny, 1:nz) = &
      reshape(rpack(n*bsize+1:(n+1)*bsize), (/nxg/PROCY, ny, nz/) ) 
  end do
  
  tmx4 = MPI_Wtime()

  tme = tmx4-tmx1
  if(tme>trs_max) trs_max = tme
  if(tme<trs_min) trs_min = tme
  trs_sum = trs_sum+tme

  tme = tmx3-tmx2
  if(tme>a2a_max) a2a_max = tme
  if(tme<a2a_min) a2a_min = tme
  a2a_sum = a2a_sum+tme
end do iter

!DEBUG -  if transpose was done correct
!error_max = 0.0
!do k = 1, nz
!  do j = 1, nyg
!    do i = 1, nxg/PROCY
!      error = abs(uy(i,j,k)-1e8*real(yid*nxg/PROCY+i)-1e4*real(j)-1e0*real(zid*nz+k))
!      if(error .gt. error_max) error_max = error
!    end do
!  end do
!end do 
!call MPI_Reduce(error_max,error_max_g,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
!if(my_rank .eq. 0) print *, 'maximum error is', error_max_g


if(my_rank==0) print *, 'All times in seconds'

call MPI_Reduce(trs_max,tme_g,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
if(my_rank==0) print *, 'Max transpose: ', tme_g 
call MPI_Reduce(trs_min,tme_g,1,MPI_REAL8,MPI_MIN,0,MPI_COMM_WORLD,ierr)
if(my_rank==0) print *, 'Min transpose: ', tme_g 
call MPI_Reduce(trs_sum,tme_g,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
if(my_rank==0) print *, 'Ave transpose: ', tme_g/real(iter_max)/real(nproc)

call MPI_Reduce(a2a_max,tme_g,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
if(my_rank==0) print *, 'Max All_to_all: ', tme_g 
call MPI_Reduce(a2a_min,tme_g,1,MPI_REAL8,MPI_MIN,0,MPI_COMM_WORLD,ierr)
if(my_rank==0) print *, 'Min All_to_all: ', tme_g 
call MPI_Reduce(a2a_sum,tme_g,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
if(my_rank==0) print *, 'Ave All_to_all: ', tme_g/real(iter_max)/real(nproc)

call MPI_Finalize(ierr)

end program test_transpose

