#include "macros.h"

module structf
implicit none

integer,parameter :: delta_num=10
integer           :: delta_val(delta_num)

! max range:  10 ... 10
integer,parameter :: pdf_max_bin=1000
real*8            :: pdf_bin_size = .01







type pdf_structure_function
!
! pdf%data(-bin:bin,delta_num,npdfs)
!
real*8,allocatable :: pdf(:,:)
integer             :: nbin         ! number of bins
integer             :: ncalls       ! number of instances PDF has been computed
                               ! (i.e. we may compute the PDF for 10 time steps
end type
 
!
!  Our structure functions
!  SF(NUM_SF,3) :  structure functions, in the x,y and z directions
!  
!
integer,parameter :: NUM_SF=3
type(pdf_structure_function) ::  SF(NUM_SF,3)


contains





subroutine init_pdf_module()
!
! Initialize this module.  Set the values of delta to use
! for structure fuctions  <u(x+delta)-u(x)>
!
! 1, 2, 4, 8 16 32 64 128 256 512 
!
implicit none

integer idel,i
do idel=1,delta_num
   delta_val(idel)=2**(idel-1)
enddo

do i=1,NUM_SF
   call init_pdf(SF(i,1),100)
   call init_pdf(SF(i,2),100)
   call init_pdf(SF(i,3),100)
enddo
end subroutine





subroutine init_pdf(str,bin)
!
! call this to initialize a pdf_structure_function
!
implicit none
type(pdf_structure_function) :: str
integer :: bin

if (allocated(str%pdf)) deallocate(str%pdf)

str%nbin=bin
allocate(str%pdf(-bin:bin,delta_num))
str%pdf=0
str%ncalls=0

end subroutine





subroutine resize_pdf(str,bin)
!
! resize a pdf_structure_function
!
implicit none
type(pdf_structure_function) :: str
integer :: bin

!local vars
type(pdf_structure_function) :: str2
integer n

n=str%nbin
if (bin<n) call abort("resize_pdf: attempting to shrink pdf");
if (bin==n) return;

! make a copy of data
allocate(str2%pdf(-n:n,delta_num))
str2%pdf=str%pdf
str2%ncalls=str%ncalls

! create a larger structure function
deallocate(str%pdf)
allocate(str%pdf(-bin:bin,delta_num))

! copy data into new, larger structure function
str%pdf=0
str%pdf(-n:n,:)=str2%pdf(-n:n,:)
str%ncalls=str2%ncalls
str%nbin=bin

! delete the copy of original data
deallocate(str2%pdf)

end subroutine








subroutine outputSF(time,fid)
use params
implicit none
real*8 time
CPOINTER fid

integer i,ierr
character*80 message
real*8 x


do i=1,NUM_SF
   call mpisum_pdf(SF(i,1))  
   call mpisum_pdf(SF(i,2))
   call mpisum_pdf(SF(i,3))
enddo

if (my_pe==io_pe) then
   call cwrite8(fid,time,1)
   ! number of structure functions
   x=3 ; call cwrite8(fid,x,1)
   do i=1,NUM_SF
      call normalize_and_write_pdf(fid,SF(i,1),SF(i,1)%nbin)   
      call normalize_and_write_pdf(fid,SF(i,2),SF(i,2)%nbin)   
      call normalize_and_write_pdf(fid,SF(i,3),SF(i,3)%nbin)   
      
      ! reset PDF's
      SF(i,1)%ncalls=0
      SF(i,1)%pdf=0
      SF(i,2)%ncalls=0
      SF(i,2)%pdf=0
      SF(i,3)%ncalls=0
      SF(i,3)%pdf=0
      
   enddo
endif


end subroutine






subroutine normalize_and_write_pdf(fid,str,nbin)
use params
implicit none
integer i,ierr,nbin
CPOINTER :: fid
type(pdf_structure_function) :: str
real*8 x,pdfdata(-nbin:nbin,delta_num)

! the delta values
x=delta_num; call cwrite8(fid,x,1)
do i=1,delta_num
   x=delta_val(i); call cwrite8(fid,x,1)
enddo

! PDF data
call cwrite8(fid,pdf_bin_size,1)
x=str%nbin; call cwrite8(fid,x,1)
x=str%ncalls; call cwrite8(fid,x,1)

! normalize
pdfdata=str%pdf
x=max(str%ncalls,1)
pdfdata=pdfdata / x;
pdfdata=pdfdata/g_nx/g_ny/g_nz

call cwrite8(fid,pdfdata,(2*nbin+1)*delta_num)
end subroutine










subroutine mpisum_pdf(str)
!
! sum PDF's over all CPUs into one PDF on cpu = io_pe
!
use params
use mpi
implicit none
type(pdf_structure_function) :: str


! local variables
type(pdf_structure_function) :: totstr
real*8,allocatable :: pdfdata(:,:)
integer :: bin,ierr,n



#ifdef USE_MPI
! find maximum value of bin  MPI max
call MPI_allreduce(str%nbin,bin,1,MPI_INTEGER,MPI_MAX,comm_3d,ierr)

! resize all str's to size bin
call resize_pdf(str,bin)

! MPI sum into totstr.
! str%pdf(-bin:bin) = 2*bin+1 data points

n=(2*bin+1)*delta_num
if (my_pe == io_pe) then
   allocate(pdfdata(-bin:bin,delta_num))
   call MPI_reduce(str%pdf,pdfdata,n,MPI_REAL8,MPI_MAX,io_pe,comm_3d,ierr)
   str%pdf=pdfdata
   deallocate(pdfdata)
else
   call MPI_reduce(str%pdf,0,n,MPI_REAL8,MPI_MAX,io_pe,comm_3d,ierr)
endif

#endif

end subroutine







subroutine compute_pdf(Q,n1,n1d,n2,n2d,n3,n3d,str)
!
! compute a pdf_structure function along the first dimension of Q
! for all the values of delta given by delta_val(:)
!
implicit none
integer :: n1,n1d,n2,n2d,n3,n3d
real*8 :: Q(n1d,n2d,n3d)
type(pdf_structure_function) :: str

! local variables
real*8  :: del
integer :: bin,idel,i,j,k,i2

do k=1,n3
   do j=1,n2
      do idel=1,delta_num
         if (delta_val(idel) < n1/2) then
            do i=1,n1
               i2 = i + delta_val(idel)
               if (i2>n1) i2=i2-n1
               
               del = Q(i,j,k)-Q(i2,j,k)
               del = del/pdf_bin_size
               bin = nint(del)

               if (abs(bin)>pdf_max_bin) then
                  print *,"Warning pdf bin overflow"
                  if (bin<0) then
                     bin=-pdf_max_bin
                  else
                     bin=pdf_max_bin
                  endif
               endif

               ! increase the size of our PDF function
               if (abs(bin)>str%nbin) call resize_pdf(str,abs(bin)) 
               str%pdf(bin,idel)=str%pdf(bin,idel)+1
            enddo
         endif
      enddo
   enddo
enddo
str%ncalls=str%ncalls+1
end subroutine





subroutine z_ifft3d_str(fin,f,strid,compy_instead_of_x,compz)
!
!  compute inverse fft 3d of fin, return in f
!  fin and f can overlap in memory
!  Also compute structure functions.
!  strid=1    U
!        2    V
!        3    W
!
!  This routine can compute PDF's of structure functions in the x direction
!  or y direction (but not both) and in the z direction.
! 
!  compy_instead_of_x    compute PDF in y direction instead of x
!  compz                 also compute PDF in z direction
!
!
use params
use fft_interface
use transpose
implicit none
real*8 fin(g_nz2,nslabx,ny_2dz)  ! input
real*8 f(nx,ny,nz)    ! output
real*8 work(nx,ny,nz) ! work array
real*8 work2(g_nz2,nslabx,ny_2dz) ! work array
logical :: compy_instead_of_x,compz
integer strid



!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k

n1=g_nz
n1d=g_nz2   	
n2=nslabx
n2d=nslabx
n3=ny_2dz
n3d=ny_2dz

work2=fin
call ifft1(work2,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(work2,f,n1,n1d,n2,n2d,n3,n3d)

if (compy_instead_of_x) then

   call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d)
   call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
   call transpose_from_x(work,f,n1,n1d,n2,n2d,n3,n3d)

   call transpose_to_y(f,work,n1,n1d,n2,n2d,n3,n3d)
   call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
   call compute_pdf(work,n1,n1d,n2,n2d,n3,n3d,SF(strid,2))
   call transpose_from_y(work,f,n1,n1d,n2,n2d,n3,n3d)
   
else

   call transpose_to_y(f,work,n1,n1d,n2,n2d,n3,n3d)
   call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
   call transpose_from_y(work,f,n1,n1d,n2,n2d,n3,n3d)
   
   
   call transpose_to_x(f,work,n1,n1d,n2,n2d,n3,n3d)
   call ifft1(work,n1,n1d,n2,n2d,n3,n3d)
   call compute_pdf(work,n1,n1d,n2,n2d,n3,n3d,SF(strid,1))
   call transpose_from_x(work,f,n1,n1d,n2,n2d,n3,n3d)
endif

if (compz) then
   call transpose_to_z(f,work,n1,n1d,n2,n2d,n3,n3d)
   call compute_pdf(work,n1,n1d,n2,n2d,n3,n3d,SF(strid,3))
endif


end









end module

