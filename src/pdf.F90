#include "macros.h"

module structf
implicit none

integer           :: structf_init=0
integer,parameter :: delta_num=10
integer           :: delta_val(delta_num)

! max range:  10 ... 10
integer,parameter :: pdf_max_bin=1000
real*8            :: pdf_bin_size = .01





type pdf_structure_function
!
! pdf%data(-bin:bin,delta_num,npdfs)
!
!SGI f90 does not allow allocatable arrays in a struct
real*8,pointer      :: pdf(:,:)
integer             :: nbin         ! number of bins
integer             :: ncalls       ! number of instances PDF has been computed
                               ! (i.e. we may compute the PDF for 10 time steps
end type
 
!
!  Our structure functions
!  SF(NUM_SF,3) :  structure functions, in the x,y and z directions
!  
!
integer,parameter :: NUM_SF=6
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

str%nbin=bin
allocate(str%pdf(-bin:bin,delta_num))
str%pdf=0
str%ncalls=0

end subroutine





subroutine resize_pdf(str,bin)
!
! resize a pdf_structure_function so it can hold value = bin
! bin should be out of the range:   -str%nbin : str%nbin
! 
! if bin is inside that range, we abort.  (cant shrink str)
! if bin is one of the end points of that range, we do nothing
! if bin is outside the range:
!      if abs(bin) > pdf_max_bin:  print warning, set bin to end point
!      else: resize str to contain bin.  
!      
!
use params
implicit none
type(pdf_structure_function) :: str
integer :: bin

!local vars
type(pdf_structure_function) :: str2
integer n

n=str%nbin
if (bin<n) call abort("resize_pdf: attempting to shrink pdf");
if (bin==n) return;


if (bin>pdf_max_bin) then
   print *,"Warning pdf bin overflow on pe=",my_pe
endif



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
character(len=80) message
real*8 x

if (structf_init==0) then
   call abort("Error: outputSF() called, but structf module not initialized")
endif

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
   call MPI_reduce(str%pdf,pdfdata,n,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   str%pdf=pdfdata
   deallocate(pdfdata)
else
   call MPI_reduce(str%pdf,0,n,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
endif

#endif

end subroutine







subroutine compute_pdf(Q,n1,n1d,n2,n2d,n3,n3d,str)
!
! compute a pdf_structure function along the first dimension of Q
! for all the values of delta given by delta_val(:)
!
use params
implicit none
integer :: n1,n1d,n2,n2d,n3,n3d,n
real*8 :: Q(n1d,n2d,n3d,3)
type(pdf_structure_function) :: str(NUM_SF)

! local variables
real*8  :: del,delv(3),delq
real*8  :: one_third = (1d0/3d0)
real*8  :: tmx1,tmx2
integer :: bin,idel,i,j,k,i2,nsf

if (structf_init==0) then
   structf_init=1
   call init_pdf_module()
endif
call wallclock(tmx1)

do k=1,n3
   do j=1,n2
      do idel=1,delta_num
         if (delta_val(idel) < n1/2) then
            do i=1,n1
               ! compute structure functions for U,V,W 
               do n=1,3
                  i2 = i + delta_val(idel)
                  if (i2>n1) i2=i2-n1
                  
                  del = Q(i2,j,k,n)-Q(i,j,k,n)
                  delv(n)=del
                  del = del/pdf_bin_size
                  bin = nint(del)
                  
                  ! increase the size of our PDF function
                  if (abs(bin)>str(n)%nbin) call resize_pdf(str(n),abs(bin)+10) 
                  if (bin>pdf_max_bin) bin=pdf_max_bin
                  if (bin<-pdf_max_bin) bin=-pdf_max_bin
                  str(n)%pdf(bin,idel)=str(n)%pdf(bin,idel)+1
               enddo

               ! compute structure functions for U(U**2+V**2+W**2)
               delq = delv(1)**2 + delv(2)**2 + delv(3)**2
               do n=1,3
                  nsf=n+3
                  del=delv(n)*delq
                  if (del>=0) then
                     del=del**one_third
                  else
                     del=-(-del)**one_third
                  endif
                  del = del/pdf_bin_size
                  bin=nint(del)
                  if (abs(bin)>str(nsf)%nbin) call resize_pdf(str(nsf),abs(bin)+10) 
                  if (bin>pdf_max_bin) bin=pdf_max_bin
                  if (bin<-pdf_max_bin) bin=-pdf_max_bin
                  str(nsf)%pdf(bin,idel)=str(nsf)%pdf(bin,idel)+1
               enddo
            enddo
         endif
      enddo
   enddo
enddo


do n=1,NUM_SF
str(n)%ncalls=str(n)%ncalls+1
enddo
call wallclock(tmx2)
tims(12)=tims(12)+(tmx2-tmx1)
end subroutine





subroutine z_ifft3d_str(fin,f,w1,w2,compx,compy,compz)
!
!  compute inverse fft 3d of fin, return in f
!  fin and f can overlap in memory
!  Also compute structure functions.
!
!  This routine can compute PDF's of structure functions in the x direction
!  or y direction (but not both) and in the z direction.
! 
!  compx   compute PDF in x direction
!  compy   compute PDF in y direction
!  compz   compute PDF in z direction
!
! only one of compx or compy can be computed per call.
! compx or compy can be computed w/o any additional transposes.
! compz requires an additional z-transpose.
!
use params
use fft_interface
use transpose
implicit none
real*8 fin(g_nz2,nslabx,ny_2dz,n_var)  ! input
real*8 f(nx,ny,nz,n_var)               ! output

! overlaped in memory:
real*8 w1(g_nz2,nslabx,ny_2dz,n_var)
real*8 w2(nx,ny,nz,n_var)    
logical :: compx,compy,compz


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k,n


do n=1,3
   n1=g_nz
   n1d=g_nz2   	
   n2=nslabx
   n2d=nslabx
   n3=ny_2dz
   n3d=ny_2dz
   w1(:,:,:,1)=fin(:,:,:,n)
   call ifft1(w1,n1,n1d,n2,n2d,n3,n3d)
   call transpose_from_z(w1,f(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo

if (compx) then
   do n=1,3
      call transpose_to_y(f(1,1,1,n),w2(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      call ifft1(w2(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      call transpose_from_y(w2(1,1,1,n),f(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      
      call transpose_to_x(f(1,1,1,n),w2(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      call ifft1(w2(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      call transpose_from_x(w2(1,1,1,n),f(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
   enddo
   call compute_pdf(w2,n1,n1d,n2,n2d,n3,n3d,SF(1,1))
   
else ! compy, or default (compute no structure functions)
   do n=1,3
      call transpose_to_x(f(1,1,1,n),w2(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      call ifft1(w2(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      call transpose_from_x(w2(1,1,1,n),f(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      
      call transpose_to_y(f(1,1,1,n),w2(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      call ifft1(w2(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      call transpose_from_y(w2(1,1,1,n),f(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
   enddo
   if (compy) call compute_pdf(w2,n1,n1d,n2,n2d,n3,n3d,SF(1,2))
endif



if (compz) then
   do n=1,3
      call transpose_to_z(f(1,1,1,n),w2(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
   enddo
   call compute_pdf(w2,n1,n1d,n2,n2d,n3,n3d,SF(1,3))
endif


end subroutine





end module

