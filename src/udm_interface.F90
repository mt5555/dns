#include "macros.h"
subroutine udm_write_uvw(fname,time,Q,Qhat,work1,bufudm)
!
! 
!
use params
use mpi
use fft_interface
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 bufudm(nx2-nx1+1, ny2-ny1+1, nz2-nz1+1)    
real*8 :: time
character(len=80) fname


!local
character(len=80) message
integer :: n,i,j,k


integer len1, len2, fidudm, dsidudm, infoid, ierr
integer attrdim, attrcnt
integer i1, j1, k1
integer*8 sizesudm(3), gsizesudm(3), offsetsudm(3)
integer*8 cpusizesudm(3)
integer*8 cpulistudm(ncpu_x * ncpu_y * ncpu_z)  
integer*4 cpulocudm(3) 


character(len=80) fnameudm
character(len=80) dsnameudm 
character(len=80) attrname 
character*(2) dotudm 

! character*(6) dbgname

integer*8 fidmpi,infoin

#ifdef HAVE_UDM

#include "UDMFCInterfaceF.h"

!  open udm file
   write(message,'(f10.4)') 10000.0000 + time

   len1 = len_trim(rundir)
   len2 = len_trim(runname)
   fnameudm = rundir(1:len1) // runname(1:len2) // message(2:10) //'.h5'// char(0) 

   if (io_pe==my_pe) then
   call print_message("opening MPI file")
   call MPI_Info_create(infoin,ierr)
   call MPI_Info_set(infoin, "striping_factor", "128",ierr) 	
   call MPI_Info_set(infoin, "striping_unit",  "8388608",ierr)
   call MPI_File_open(comm_1,fnameudm, &
           MPI_MODE_WRONLY + MPI_MODE_CREATE ,&
           infoin, fidmpi,ierr)
   call print_message("closing MPI file")
   call MPI_File_close(fidmpi,ierr)
   call print_message("done")
   endif

   call print_message("opening file")
   call UDM_FILE_OPEN(fnameudm, UDM_HDF_OPEN_CREATE, fidudm, ierr)
   call print_message("done with UDM open")



   dotudm = '.'//char(0)
   attrname = 'time'//char(0)
   attrdim = 1
   attrcnt = 1
   call UDM_ATTR_WRITE(attrname,dotudm,fidudm,UDM_ATOMIC_DOUBLE, 1, 1, time,ierr)  

!   i = ncpu_x *(nx2 - nx1 + 1)  
!   j = ncpu_y *(ny2 - ny1 + 1)  
!   k = ncpu_z *(nz2 - nz1 + 1)  
!   attrname = 'nx'//char(0) 
!   call UDM_ATTR_WRITE(attrname,dotudm,fidudm,UDM_ATOMIC_INT,attrdim,attrcnt,i,ierr)
!   attrname = 'ny'//char(0)     
!   call UDM_ATTR_WRITE(attrname,dotudm,fidudm,UDM_ATOMIC_INT,attrdim,attrcnt,j,ierr)
!   attrname = 'nz'//char(0)     
!   call UDM_ATTR_WRITE(attrname,dotudm,fidudm,UDM_ATOMIC_INT,attrdim,attrcnt,k,ierr)
!
!  open a group for coordinates
!
!   gnameudm = 'coords'//char(0)
!   call UDM_GROUP_OPEN(gnameudm, fidudm, UDM_HDF_OPEN_CREATE, gidudm, ierr)
!   attrname = 'delx'//char(0)
!   call UDM_ATTR_WRITE(attrname,dotudm,fidudm,UDM_ATOMIC_DOUBLE,attrdim,attrcnt,delx,ierr)
!   attrname = 'dely'//char(0)
!   call UDM_ATTR_WRITE(attrname,dotudm,fidudm,UDM_ATOMIC_DOUBLE,attrdim,attrcnt,dely,ierr)
!   attrname = 'delz'//char(0)
!   call UDM_ATTR_WRITE(attrname,dotudm,fidudm,UDM_ATOMIC_DOUBLE,attrdim,attrcnt,delz,ierr)

!  create a empty structure for a dataset

   dsnameudm(1:1) = char(0)
   call UDM_INFO_CREATE(dsnameudm, dsnameudm, -1, UDM_INFO_DATA_DIST, infoid, ierr)
   call UDM_INFO_ITEM_SET(infoid, UDM_INFO_NDIMS, 1, ierr)

!  only my_pe = 0 writes corrdinates

   gsizesudm(1) = o_nx 
   offsetsudm(1) = 0
   if (my_pe .eq. 0) then  
      sizesudm(1) = gsizesudm(1)
   else 
      sizesudm(1) = 0
   endif
   call UDM_INFO_ITEM_SET(infoid, UDM_INFO_DIMS, sizesudm, ierr)
   call UDM_INFO_ITEM_SET(infoid, UDM_INFO_DIMS_TOTAL, gsizesudm, ierr)
   call UDM_INFO_ITEM_SET(infoid, UDM_INFO_STARTS, offsetsudm, ierr)  
   dsnameudm = 'xcord'//char(0)
   call UDM_DATASET_WRITE(dsnameudm,fidudm,UDM_ATOMIC_DOUBLE,infoid,g_xcord(1),ierr)

   gsizesudm(1) = o_ny
   if (my_pe .eq. 0) then
      sizesudm(1) = gsizesudm(1)
   endif 
   call UDM_INFO_ITEM_SET(infoid, UDM_INFO_DIMS,  sizesudm,  ierr)
   call UDM_INFO_ITEM_SET(infoid, UDM_INFO_DIMS_TOTAL, gsizesudm, ierr)
   dsnameudm = 'ycord'//char(0)
   call UDM_DATASET_WRITE(dsnameudm,fidudm,UDM_ATOMIC_DOUBLE,infoid,g_ycord(1),ierr)

   gsizesudm(1) = o_nz
   if (my_pe .eq. 0) then
      sizesudm(1) = gsizesudm(1)
   endif 
   call UDM_INFO_ITEM_SET(infoid, UDM_INFO_DIMS,  sizesudm,  ierr)
   call UDM_INFO_ITEM_SET(infoid, UDM_INFO_DIMS_TOTAL, gsizesudm, ierr)
   dsnameudm = 'zcord'//char(0) 
   call UDM_DATASET_WRITE(dsnameudm,fidudm,UDM_ATOMIC_DOUBLE,infoid,g_zcord(1),ierr)

   call UDM_INFO_FREE(1,infoid,ierr)







!     since the current fortran interface treats the rightmost index first, as in C, 
!     I have to reverse the orders of indices 

!      gnameudm = 'ns_uvw'//char(0)
!      call UDM_GROUP_OPEN(gnameudm, fidudm, UDM_HDF_OPEN_CREATE, gidudm, ierr)

      dsnameudm(1:1) = char(0)
      call UDM_INFO_CREATE(dsnameudm, dsnameudm, -1, UDM_INFO_DATA_DIST, infoid, ierr)
      call UDM_INFO_ITEM_SET(infoid, UDM_INFO_NDIMS, 3, ierr)

      sizesudm(3) = nx2
      sizesudm(2) = ny2
      sizesudm(1) = nz2
      call UDM_INFO_ITEM_SET(infoid, UDM_INFO_DIMS, sizesudm, ierr);

      cpusizesudm(3) = ncpu_x
      cpusizesudm(2) = ncpu_y
      cpusizesudm(1) = ncpu_z

      cpulocudm(3) = my_x
      cpulocudm(2) = my_y
      cpulocudm(1) = my_z 

      call UDM_INFO_ITEM_SET(infoid, UDM_INFO_PGRID_DIMS, cpusizesudm, ierr)
      call FIO_FGET_PLIST(cpusizesudm, cpulocudm, 3, cpulistudm, ierr)
      call UDM_INFO_ITEM_SET(infoid, UDM_INFO_PGRID_ORDER, cpulistudm, ierr)

      dsnameudm = 'u'//char(0)
      do k = nz1, nz2 
      k1 = k - nz1 + 1  
      do j = ny1, ny2 
      j1 = j - ny1 + 1 
      do i = nx1, nx2 
      i1 = i - nx1 + 1  
      bufudm(i1,j1,k1) = Q(i,j,k,1) 
      enddo
      enddo
      enddo 
      call UDM_DATASET_WRITE(dsnameudm,fidudm,UDM_ATOMIC_DOUBLE,infoid,bufudm,ierr)

      dsnameudm = 'v'//char(0)
      do k = nz1, nz2        
      k1 = k - nz1 + 1
      do j = ny1, ny2        
      j1 = j - ny1 + 1
      do i = nx1, nx2        
      i1 = i - nx1 + 1 
      bufudm(i1,j1,k1) = Q(i,j,k,2)
      enddo
      enddo
      enddo
      call UDM_DATASET_WRITE(dsnameudm,fidudm,UDM_ATOMIC_DOUBLE,infoid,bufudm,ierr)

      dsnameudm = 'w'//char(0)
      do k = nz1, nz2        
      k1 = k - nz1 + 1
      do j = ny1, ny2        
      j1 = j - ny1 + 1
      do i = nx1, nx2        
      i1 = i - nx1 + 1 
      bufudm(i1,j1,k1) = Q(i,j,k,3)
      enddo
      enddo
      enddo
      call UDM_DATASET_WRITE(dsnameudm,fidudm,UDM_ATOMIC_DOUBLE,infoid,bufudm,ierr)

      CALL UDM_INFO_FREE(1, infoid, ierr)

   call UDM_FILE_CLOSE(fidudm, ierr)

#else
   call abort("udm_write_uvw: Error: UDM support not compiled in")
#endif
end subroutine




subroutine udm_read_uvw(Q,Qhat,work1,bufudm)
!
! 
!
use params
use mpi
use fft_interface
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 bufudm(nx2-nx1+1, ny2-ny1+1, nz2-nz1+1)


!local
character(len=80) message
character(len=80) fname
integer :: n

integer fidudm, dsidudm, infoid, ierr
integer i, j, k, i1, j1, k1
integer dimsudm 
integer*8 sizesudm(3), gsizesudm(3), offsetsudm(3)

character(len=80) fnameudm
character(len=80) attrname 
character(len=80) dsnameudm
character(len=80) dotudm
character*2 nochar 
! character*(5) dbgname 

#ifdef HAVE_UDM

#include "UDMFCInterfaceF.h"


Q=0

   dotudm = '.'//char(0)
   attrname = 'time'//char(0)
   nochar(1:1) = char(0)
   fnameudm = rundir(1:len_trim(rundir)) // 'restart.h5'//char(0)
   call print_message("Restarting from UDM file:")
   call print_message(fnameudm)
   call UDM_FILE_OPEN(fnameudm, UDM_HDF_OPEN_READ, fidudm, ierr)

   dotudm = '.'//char(0)
   attrname = 'time'//char(0)
   
   call UDM_ATTR_READ(attrname,dotudm,fidudm,UDM_ATOMIC_DOUBLE,time_initial,ierr)

!   attrname = 'nx'//char(0)
!   call UDM_ATTR_READ(attrname,dotudm,fidudm,UDM_ATOMIC_INT, nxstored, ierr)
!   attrname = 'ny'//char(0)
!   call UDM_ATTR_READ(attrname,dotudm,fidudm,UDM_ATOMIC_INT, nystored, ierr)
!   attrname = 'nz'//char(0)
!   call UDM_ATTR_READ(attrname,dotudm,fidudm,UDM_ATOMIC_INT, nzstored, ierr)
!
!   if (nxstored .ne. ncpu_x *(nx2 - nx1 + 1))  &   
!       call abort("ERROR: data file nx <> nx set in params.h")
!   if (nystored .ne. ncpu_y *(ny2 - ny1 + 1))  & 
!       call abort("ERROR: data file ny <> ny set in params.h")
!   if (nzstored .ne. ncpu_z *(nz2 - nz1 + 1))  & 
!       call abort("ERROR: data file nz <> nz set in params.h")

!  read g_xcord, g_ycor and g_zcord

!   dsnameudm(1:12) = 'coords/xcord'
!   dsnameudm(13:13) = char(0)
!   call FIO_INFO_CREATE(fidudm,dsnameudm,UDM_INFOTYPE_DATASET,infoid,ierr) 
!   call FIO_INFO_ITEM_GET(infoid, UDM_INFO_DIMS,   dimsudm,   ierr)
!   if (dimsudm .ne. 1) call abort("ERROR: datafile g_xcord wrong")
!   call FIO_INFO_ITEM_GET(infoid, UDM_INFO_GSIZES, sizesudm, ierr)
!   offsetsudm(1) = 0
!   call FIO_INFO_ITEM_SET(infoid, UDM_INFO_SIZES, sizesudm, ierr) 
!   call FIO_INFO_ITEM_SET(infoid, UDM_INFO_OFFSETS, offsetsudm, ierr) 
!   call FIO_READ(fidudm, dsnameudm, infoid, UDM_DOUBLE_DATASET, g_xcord(1), ierr)
!   call FIO_INFO_FREE(1, infoid, ierr)
!
!   dsnameudm(1:12) = 'coords/ycord'
!   call FIO_INFO_CREATE(fidudm,dsnameudm,UDM_INFOTYPE_DATASET,infoid,ierr)
!   call FIO_INFO_ITEM_GET(infoid, UDM_INFO_DIMS,   dimsudm,   ierr)
!   if (dimsudm .ne. 1) call abort("ERROR: datafile g_ycord wrong")
!   call FIO_INFO_ITEM_GET(infoid, UDM_INFO_GSIZES, sizesudm, ierr)
!   call FIO_INFO_ITEM_SET(infoid, UDM_INFO_SIZES, sizesudm, ierr)  
!   call FIO_INFO_ITEM_SET(infoid, UDM_INFO_OFFSETS, offsetsudm, ierr)  
!   call FIO_READ(fidudm, dsnameudm, infoid, UDM_DOUBLE_DATASET, g_ycord(1), ierr)
!   call FIO_INFO_FREE(1, infoid, ierr)
!
!   dsnameudm(1:12) = 'coords/zcord'
!   call FIO_INFO_CREATE(fidudm,dsnameudm,UDM_INFOTYPE_DATASET,infoid,ierr)
!   call FIO_INFO_ITEM_GET(infoid, UDM_INFO_DIMS,   dimsudm,   ierr)
!   if (dimsudm .ne. 1) call abort("ERROR: datafile g_zcord wrong")
!   call FIO_INFO_ITEM_GET(infoid, UDM_INFO_GSIZES, sizesudm, ierr)
!   call FIO_INFO_ITEM_SET(infoid, UDM_INFO_SIZES, sizesudm, ierr)  
!   call FIO_INFO_ITEM_SET(infoid, UDM_INFO_OFFSETS, offsetsudm, ierr)  
!   call FIO_READ(fidudm, dsnameudm, infoid, UDM_DOUBLE_DATASET, g_zcord(1), ierr)
!   call FIO_INFO_FREE(1, infoid, ierr)

   dsnameudm = 'u'//char(0)
   call UDM_INFO_CREATE(nochar, dsnameudm,fidudm,UDM_INFO_DATA_DIST, infoid, ierr)
   CALL UDM_INFO_ITEM_GET(infoid, UDM_INFO_NDIMS, dimsudm, ierr)
   if (dimsudm .ne. 3) call abort("ERROR: datafile is wrong")
   CALL UDM_INFO_ITEM_GET(infoid, UDM_INFO_DIMS_TOTAL, gsizesudm, ierr)

   if (gsizesudm(1) .ne. (ncpu_z *(nz2 - nz1 + 1)))    &
       call abort("ERROR: nz, ncpu_z in param.h not consistent the datafile")
   if (gsizesudm(2) .ne. (ncpu_y *(ny2 - ny1 + 1)))    &
       call abort("ERROR: ny, ncpu_y in param.h not consistent the datafile")
   if (gsizesudm(3) .ne. (ncpu_x *(nx2 - nx1 + 1)))    &
       call abort("ERROR: nx, ncpu_x in param.h not consistent the datafile")

   sizesudm(1) = gsizesudm(1) / ncpu_z  
   sizesudm(2) = gsizesudm(2) / ncpu_y
   sizesudm(3) = gsizesudm(3) / ncpu_x
   if ( (gsizesudm(1) - sizesudm(1) * ncpu_z .ne. 0) .or.   &
        (gsizesudm(2) - sizesudm(2) * ncpu_y .ne. 0) .or.   &
        (gsizesudm(3) - sizesudm(3) * ncpu_x .ne. 0) )      &
        call abort("ERROR: mod(nx, ncpu_x) != 0")
   offsetsudm(1) = my_z * sizesudm(1)
   offsetsudm(2) = my_y * sizesudm(2)
   offsetsudm(3) = my_x * sizesudm(3)

!   call UDM_INFO_ITEM_SET(infoid, UDM_INFO_NDIMS, dimsudm, ierr)
   call UDM_INFO_ITEM_SET(infoid, UDM_INFO_DIMS, sizesudm, ierr)
   call UDM_INFO_ITEM_SET(infoid, UDM_INFO_STARTS, offsetsudm, ierr)

   CALL UDM_DATASET_READ(dsnameudm,fidudm,UDM_ATOMIC_DOUBLE,infoid,bufudm,ierr)

   do k = nz1, nz2  
   k1 = k - nz1 + 1  
   do j = ny1, ny2   
   j1 = j - ny1 + 1   
   do i = nx1, nx2 
   i1 = i - nx1 + 1  
   Q(i,j,k,1) = bufudm(i1,j1,k1)
   enddo
   enddo
   enddo
   
   dsnameudm = 'v'//char(0)
   CALL UDM_DATASET_READ(dsnameudm,fidudm,UDM_ATOMIC_DOUBLE,infoid,bufudm,ierr)
   do k = nz1, nz2       
   k1 = k - nz1 + 1 
   do j = ny1, ny2      
   j1 = j - ny1 + 1  
   do i = nx1, nx2      
   i1 = i - nx1 + 1 
   Q(i,j,k,2) = bufudm(i1,j1,k1)
   enddo
   enddo
   enddo
 
   dsnameudm = 'w'//char(0)
   CALL UDM_DATASET_READ(dsnameudm,fidudm,UDM_ATOMIC_DOUBLE,infoid,bufudm,ierr)
   do k = nz1, nz2       
   k1 = k - nz1 + 1 
   do j = ny1, ny2      
   j1 = j - ny1 + 1  
   do i = nx1, nx2      
   i1 = i - nx1 + 1 
   Q(i,j,k,3) = bufudm(i1,j1,k1)
   enddo
   enddo
   enddo

   call UDM_INFO_FREE(1,infoid,ierr)
   call UDM_FILE_CLOSE(fidudm, ierr)


#else
   call abort("udm_write_uvw: Error: UDM support not compiled in")
#endif

end subroutine











subroutine udm_read_uvw_old(Q,Qhat,work1,bufudm)
!
! 
!
use params
use mpi
use fft_interface
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: Qhat(nx,ny,nz,n_var)
real*8 :: work1(nx,ny,nz)
real*8 bufudm(nx2-nx1+1, ny2-ny1+1, nz2-nz1+1)


!local
character(len=80) message
character(len=80) fname
integer :: n

integer fidudm, dsidudm, infoid, ierr
integer i, j, k, i1, j1, k1
integer nxstored, nystored, nzstored
integer dimsudm 
integer*8 sizesudm(3), gsizesudm(3), offsetsudm(3)

character(len=80) fnameudm
character(len=80) attrname 
character(len=80) dsnameudm
character(len=80) dotudm
character*2 nochar 
! character*(5) dbgname 

#ifdef HAVE_UDM

#include "UDMFCInterfaceF.h"


Q=0

   dotudm = '.'//char(0)
   attrname = 'time'//char(0)
   nochar(1:1) = char(0)
   fnameudm = rundir(1:len_trim(rundir)) // 'restart.h5'//char(0)
   call print_message("Restarting from UDM file:")
   call print_message(fnameudm)
   call UDM_FILE_OPEN(fnameudm, UDM_HDF_OPEN_READ, fidudm, ierr)

   dotudm = '.'//char(0)
   attrname = 'time'//char(0)
   
   call UDM_ATTR_READ(attrname,dotudm,fidudm,UDM_ATOMIC_DOUBLE,time_initial,ierr)

   attrname = 'nx'//char(0)
   call UDM_ATTR_READ(attrname,dotudm,fidudm,UDM_ATOMIC_INT, nxstored, ierr)
   attrname = 'ny'//char(0)
   call UDM_ATTR_READ(attrname,dotudm,fidudm,UDM_ATOMIC_INT, nystored, ierr)
   attrname = 'nz'//char(0)
   call UDM_ATTR_READ(attrname,dotudm,fidudm,UDM_ATOMIC_INT, nzstored, ierr)

   if (nxstored .ne. ncpu_x *(nx2 - nx1 + 1))  &   
       call abort("ERROR: data file nx <> nx set in params.h")
   if (nystored .ne. ncpu_y *(ny2 - ny1 + 1))  & 
       call abort("ERROR: data file ny <> ny set in params.h")
   if (nzstored .ne. ncpu_z *(nz2 - nz1 + 1))  & 
       call abort("ERROR: data file nz <> nz set in params.h")

!  read g_xcord, g_ycor and g_zcord

!   dsnameudm(1:12) = 'coords/xcord'
!   dsnameudm(13:13) = char(0)
!   call FIO_INFO_CREATE(fidudm,dsnameudm,UDM_INFOTYPE_DATASET,infoid,ierr) 
!   call FIO_INFO_ITEM_GET(infoid, UDM_INFO_DIMS,   dimsudm,   ierr)
!   if (dimsudm .ne. 1) call abort("ERROR: datafile g_xcord wrong")
!   call FIO_INFO_ITEM_GET(infoid, UDM_INFO_GSIZES, sizesudm, ierr)
!   offsetsudm(1) = 0
!   call FIO_INFO_ITEM_SET(infoid, UDM_INFO_SIZES, sizesudm, ierr) 
!   call FIO_INFO_ITEM_SET(infoid, UDM_INFO_OFFSETS, offsetsudm, ierr) 
!   call FIO_READ(fidudm, dsnameudm, infoid, UDM_DOUBLE_DATASET, g_xcord(1), ierr)
!   call FIO_INFO_FREE(1, infoid, ierr)
!
!   dsnameudm(1:12) = 'coords/ycord'
!   call FIO_INFO_CREATE(fidudm,dsnameudm,UDM_INFOTYPE_DATASET,infoid,ierr)
!   call FIO_INFO_ITEM_GET(infoid, UDM_INFO_DIMS,   dimsudm,   ierr)
!   if (dimsudm .ne. 1) call abort("ERROR: datafile g_ycord wrong")
!   call FIO_INFO_ITEM_GET(infoid, UDM_INFO_GSIZES, sizesudm, ierr)
!   call FIO_INFO_ITEM_SET(infoid, UDM_INFO_SIZES, sizesudm, ierr)  
!   call FIO_INFO_ITEM_SET(infoid, UDM_INFO_OFFSETS, offsetsudm, ierr)  
!   call FIO_READ(fidudm, dsnameudm, infoid, UDM_DOUBLE_DATASET, g_ycord(1), ierr)
!   call FIO_INFO_FREE(1, infoid, ierr)
!
!   dsnameudm(1:12) = 'coords/zcord'
!   call FIO_INFO_CREATE(fidudm,dsnameudm,UDM_INFOTYPE_DATASET,infoid,ierr)
!   call FIO_INFO_ITEM_GET(infoid, UDM_INFO_DIMS,   dimsudm,   ierr)
!   if (dimsudm .ne. 1) call abort("ERROR: datafile g_zcord wrong")
!   call FIO_INFO_ITEM_GET(infoid, UDM_INFO_GSIZES, sizesudm, ierr)
!   call FIO_INFO_ITEM_SET(infoid, UDM_INFO_SIZES, sizesudm, ierr)  
!   call FIO_INFO_ITEM_SET(infoid, UDM_INFO_OFFSETS, offsetsudm, ierr)  
!   call FIO_READ(fidudm, dsnameudm, infoid, UDM_DOUBLE_DATASET, g_zcord(1), ierr)
!   call FIO_INFO_FREE(1, infoid, ierr)

   dsnameudm = 'ns_uvw/u'//char(0)
   call UDM_INFO_CREATE(nochar, dsnameudm,fidudm,UDM_INFO_DATA_DIST, infoid, ierr)
   CALL UDM_INFO_ITEM_GET(infoid, UDM_INFO_NDIMS, dimsudm, ierr)
   if (dimsudm .ne. 3) call abort("ERROR: datafile is wrong")
   CALL UDM_INFO_ITEM_GET(infoid, UDM_INFO_DIMS_TOTAL, gsizesudm, ierr)

   if (gsizesudm(1) .ne. (ncpu_z *(nz2 - nz1 + 1)))    &
       call abort("ERROR: nz, ncpu_z in param.h not consistent the datafile")
   if (gsizesudm(2) .ne. (ncpu_y *(ny2 - ny1 + 1)))    &
       call abort("ERROR: ny, ncpu_y in param.h not consistent the datafile")
   if (gsizesudm(3) .ne. (ncpu_x *(nx2 - nx1 + 1)))    &
       call abort("ERROR: nx, ncpu_x in param.h not consistent the datafile")

   sizesudm(1) = gsizesudm(1) / ncpu_z  
   sizesudm(2) = gsizesudm(2) / ncpu_y
   sizesudm(3) = gsizesudm(3) / ncpu_x
   if ( (gsizesudm(1) - sizesudm(1) * ncpu_z .ne. 0) .or.   &
        (gsizesudm(2) - sizesudm(2) * ncpu_y .ne. 0) .or.   &
        (gsizesudm(3) - sizesudm(3) * ncpu_x .ne. 0) )      &
        call abort("ERROR: mod(nx, ncpu_x) != 0")
   offsetsudm(1) = my_z * sizesudm(1)
   offsetsudm(2) = my_y * sizesudm(2)
   offsetsudm(3) = my_x * sizesudm(3)

!   call UDM_INFO_ITEM_SET(infoid, UDM_INFO_NDIMS, dimsudm, ierr)
   call UDM_INFO_ITEM_SET(infoid, UDM_INFO_DIMS, sizesudm, ierr)
   call UDM_INFO_ITEM_SET(infoid, UDM_INFO_STARTS, offsetsudm, ierr)

   CALL UDM_DATASET_READ(dsnameudm,fidudm,UDM_ATOMIC_DOUBLE,infoid,bufudm,ierr)

   do k = nz1, nz2  
   k1 = k - nz1 + 1  
   do j = ny1, ny2   
   j1 = j - ny1 + 1   
   do i = nx1, nx2 
   i1 = i - nx1 + 1  
   Q(i,j,k,1) = bufudm(i1,j1,k1)
   enddo
   enddo
   enddo
   
   dsnameudm = 'ns_uvw/v'//char(0)
   CALL UDM_DATASET_READ(dsnameudm,fidudm,UDM_ATOMIC_DOUBLE,infoid,bufudm,ierr)
   do k = nz1, nz2       
   k1 = k - nz1 + 1 
   do j = ny1, ny2      
   j1 = j - ny1 + 1  
   do i = nx1, nx2      
   i1 = i - nx1 + 1 
   Q(i,j,k,2) = bufudm(i1,j1,k1)
   enddo
   enddo
   enddo
 
   dsnameudm = 'ns_uvw/w'//char(0)
   CALL UDM_DATASET_READ(dsnameudm,fidudm,UDM_ATOMIC_DOUBLE,infoid,bufudm,ierr)
   do k = nz1, nz2       
   k1 = k - nz1 + 1 
   do j = ny1, ny2      
   j1 = j - ny1 + 1  
   do i = nx1, nx2      
   i1 = i - nx1 + 1 
   Q(i,j,k,3) = bufudm(i1,j1,k1)
   enddo
   enddo
   enddo

   call UDM_INFO_FREE(1,infoid,ierr)
   call UDM_FILE_CLOSE(fidudm, ierr)


#else
   call abort("udm_write_uvw: Error: UDM support not compiled in")
#endif

end subroutine
