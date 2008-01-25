#include "macros.h"

module pdf
use params
implicit none

!
!
! compute velocity and epsilon PDFs  and joint PDFs
! 
! 
!
logical ::  compute_uvw_pdfs=.true.
logical ::  compute_uvw_jpdfs=.false.
logical ::  compute_passive_pdfs=.true.

real*8 :: uscale=.01             ! bin size for vel increment
real*8 :: epsscale=.01           ! bin size for epsilon increment
real*8 :: pscale=.0025           ! bin size for scalar increment
                                 ! for t-mix problem, scalar will 
                                 ! range from 0..1


real*8 :: core_data(g_nx,n_var)      ! used to take x direction cores thru data
real*8 :: core_ddata(g_nx,3,n_var)   ! used to take x direction cores of derivativs
integer :: core_num                  ! number of cores taken

integer           :: structf_init=0
integer,parameter :: delta_num_max=100
integer           :: delta_val(delta_num_max)
integer :: numx=0,numy=0,numz=0          ! actual number of delta's in each direction
                                         ! allowed by the grid

! max range:  10 ... 10
integer,parameter :: pdf_max_bin=4000
integer,parameter :: jpdf_max_bin=200


type pdf_structure_function
!
! pdf%data(-bin:bin,delta_num)
!
!SGI f90 does not allow allocatable arrays in a struct
real*8,pointer      :: pdf(:,:)
integer             :: delta_num    ! number of different values of delta
integer             :: nbin         ! number of bins
real*8              :: pdf_bin_size ! size of bin
integer             :: ncalls       ! number of instances PDF has been computed
                               ! (i.e. we may compute the PDF for 10 time steps
end type

type jpdf_structure_function
!
! JOINT PDF's
! pdf%data(-bin:bin,-bin:bin,delta_num)
!
!SGI f90 does not allow allocatable arrays in a struct
real*8,pointer      :: pdf(:,:,:)
integer             :: delta_num    ! number of different values of delta
integer             :: nbin         ! number of bins
real*8              :: pdf_bin_size ! size of bin
integer             :: ncalls       ! number of instances PDF has been computed
                               ! (i.e. we may compute the PDF for 10 time steps
end type
 

!
!  Our structure functions
!  SF(NUM_SF,3) :  structure functions, in the x,y and z directions
!  
!  j=1..3
!  i=1..3
!  SF(i,j)     PDF of delta_j(u_i)     i'th component of u, j'th direction
!
!  3rd order:
!  SF(i+3,j)   i<>j:   delta_j(u_j) delta_j(u_i**2)                     
!              i==j:   delta_j(u_j) delta_j(u_1**2 + u_2**2 + u_3**2)   
!
! 4th order
!  SF(i+6)    (i,j)
!             (1,1)     delta_j(u_1**2) delta_j(u_2**2)                 
!             (2,1)     delta_j(u_1**2) delta_j(u_3**2)                 
!             (1,2)     delta_j(u_2**2) delta_j(u_1**2)           
!             (2,2)     delta_j(u_2**2) delta_j(u_3**2)
!             (1,3)     delta_j(u_3**2) delta_j(u_1**2)
!             (2,3)     delta_j(u_3**2) delta_j(u_2**2)
!              
! Note: the PDF's are further normalized so that they all have the same
! units.  3rd order SF are raised to the 1/3 power, 4th order SF
! are raised to the 1/4 power before binning up the PDF.  
!
!integer,parameter :: NUM_SF=8
integer,parameter :: NUM_SF=3   ! just compute first 3 pdfs
type(pdf_structure_function) ::  SF(NUM_SF,3)
type(pdf_structure_function) ::  epsilon
! number of scalar structure functions:  params.F90::npassive 
! which is never more than n_var-2
type(pdf_structure_function) ::  SCALARS(n_var-2)  

! pdfs of a generic scalar
! convert program will use these to compute PDFs of any quantity of interest
type(pdf_structure_function),allocatable ::  cpdf(:)
integer :: number_of_cpdf = 0

integer :: overflow=0    ! count the number of overflow messages
integer :: joverflow=0    ! count the number of overflow messages
real*8  :: one_third = (1d0/3d0)


!
!  joint PDF's 
!
!     delx(u) delx(v)
!     delx(u) delx(w)
!     delx(v) delx(w)
!
!     dely(u) dely(v)
!     dely(u) dely(w)
!     dely(v) dely(w)
!
!     delz(u) delz(v)
!     delz(u) delz(w)
!     delz(v) delz(w)
!
integer,parameter :: NUM_JPDF=3
type(jpdf_structure_function) ::  jpdf_v(NUM_JPDF,3)



contains




subroutine init_pdf_module()
!
! Initialize this module.  Set the values of delta to use
! for structure fuctions  <u(x+delta)-u(x)>
!
! 1, 2, 4, 8 16 32 64 128 256 512 1024 2048  4096 
!
use params
implicit none

integer idel,i,j
integer :: nmax,max_delta
real*8  :: xtmp
! we can compute them up to g_nx/2, but no reason to go that far.
nmax=max(g_nx/3,g_ny/3,g_nz/3)


delta_val=100000

! compute all deltas up to ndelta_max
! (but we only use deltas up to ndelta)
! 1..16
do i=1,16
   delta_val(i)=i
enddo
j=16

! 18..32  (start with 2*9)
do i=9,16
   j=j+1
   delta_val(j)=2*i
enddo

! 36..64  (start with 4*9)
do i=9,16
   j=j+1
   delta_val(j)=4*i
enddo

! 72..128 (start with 8*9)
do i=9,16
   j=j+1
   delta_val(j)=8*i
enddo

! 144..256  (start with 16*9)
do i=9,16
   j=j+1
   delta_val(j)=16*i
enddo

! 288..512  (start with 32*9)
do i=9,16
   j=j+1
   delta_val(j)=32*i
enddo

! 576..1024  (start with 64*9)
do i=9,16
   j=j+1
   delta_val(j)=64*i
enddo

! 1152..2048
do i=9,16
   j=j+1
   delta_val(j)=128*i
enddo
! determine maximum value of delta to use for this grid
if (j > delta_num_max) then
   call abortdns("structf init: j > delta_num_max")
endif

max_delta = g_nmin/2
if (user_specified_isodel>0) then
   max_delta=min(max_delta,user_specified_isodel)
endif

do idel=1,delta_num_max
   if (delta_val(idel) < max_delta) then
      numx=idel
   endif
   if (delta_val(idel) < max_delta) then
      numy=idel
   endif
   if (delta_val(idel) < max_delta) then
      numz=idel
   endif
enddo

if (compute_uvw_pdfs) then
do i=1,NUM_SF
   call init_pdf(SF(i,1),100,uscale,numx)
   call init_pdf(SF(i,2),100,uscale,numy)
   call init_pdf(SF(i,3),100,uscale,numz)
enddo
call init_pdf(epsilon,100,epsscale,1)
endif

if (npassive==0) compute_passive_pdfs=.false.
if (compute_passive_pdfs) then
do i=1,npassive
   call init_pdf(SCALARS(i),100,pscale,1)
enddo
endif

if (number_of_cpdf > 0) then
   allocate(cpdf(number_of_cpdf))
   xtmp = 0 ! binsize will be set later
   do i=1,number_of_cpdf
      call init_pdf(cpdf(i),100, xtmp, 1)
   enddo
endif



if (compute_uvw_jpdfs) then
do i=1,NUM_JPDF
   call init_jpdf(jpdf_v(i,1),100,.1d0, min(numx,numy))
   call init_jpdf(jpdf_v(i,2),100,.1d0, min(numx,numz))
   call init_jpdf(jpdf_v(i,3),100,.1d0, min(numy,numz))
enddo
endif


core_num=0
end subroutine






subroutine init_pdf(str,bin,binsize,ndelta)
!
! call this to initialize a pdf_structure_function
!
implicit none
type(pdf_structure_function) :: str
integer :: bin,ndelta
real*8 binsize

str%nbin=bin
str%pdf_bin_size=binsize
allocate(str%pdf(-bin:bin,ndelta))
str%pdf=0
str%ncalls=0
str%delta_num=ndelta

end subroutine







subroutine init_jpdf(str,bin,binsize,ndelta)
!
! call this to initialize a JOINT pdf_structure_function
!
implicit none
type(jpdf_structure_function) :: str
integer :: bin,ndelta
real*8 binsize

str%nbin=bin
str%pdf_bin_size=binsize
allocate(str%pdf(-bin:bin,-bin:bin,ndelta))
str%pdf=0
str%ncalls=0
str%delta_num=ndelta

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
real*8, allocatable :: pdfdata(:,:)
integer n,ndelta

ndelta=str%delta_num

n=str%nbin
if (bin<n) call abortdns("resize_pdf: attempting to shrink pdf");
if (bin==n) return;


if (bin>pdf_max_bin) then
   overflow=overflow+1
   if (overflow<10) print *,"Warning pdf bin overflow on pe=",my_pe
   if (overflow==10) print *,"disabling bin overflow messages"
   bin=pdf_max_bin
endif



! make a copy of data
allocate(pdfdata(-n:n,ndelta))
pdfdata=str%pdf

! create a larger structure function
deallocate(str%pdf)
allocate(str%pdf(-bin:bin,ndelta))

! copy data into new, larger structure function
str%pdf=0
str%pdf(-n:n,:)=pdfdata(-n:n,:)
str%nbin=bin

! delete the copy of original data
deallocate(pdfdata)

end subroutine









subroutine resize_jpdf(str,bin)
!
! resize a jpdf_structure_function so it can hold value = bin
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
type(jpdf_structure_function) :: str
integer :: bin

!local vars
real*8, allocatable :: pdfdata(:,:,:)
integer n,ndelta

ndelta=str%delta_num

n=str%nbin
if (bin<n) call abortdns("resize_pdf: attempting to shrink pdf");
if (bin==n) return;


if (bin>jpdf_max_bin) then
   joverflow=joverflow+1
   if (joverflow<10) print *,"Warning jpdf bin overflow on pe=",my_pe
   if (joverflow==10) print *,"disabling bin overflow messages"
   bin=pdf_max_bin
endif



! make a copy of data
allocate(pdfdata(-n:n,-n:n,ndelta))
pdfdata=str%pdf

! create a larger structure function
deallocate(str%pdf)
allocate(str%pdf(-bin:bin,-bin:bin,ndelta))

! copy data into new, larger structure function
str%pdf=0
str%pdf(-n:n,-n:n,:)=pdfdata(-n:n,-n:n,:)
str%nbin=bin

! delete the copy of original data
deallocate(pdfdata)

end subroutine








subroutine output_pdf(time,fid,fidj,fidS,fidC,fidcore)
use params
implicit none
real*8 time
CPOINTER fid   ! velcoity PDFs
CPOINTER fidj  ! joint velocity PDFs
CPOINTER fidS  ! passive scalar PDFs
CPOINTER fidC  ! cpdf scalar PDFs
CPOINTER fidcore  ! x-direction core data

integer i,j,ierr,ndelta,numm,n
character(len=80) message
real*8 x

if (structf_init==0) then
   call abortdns("Error: output_pdf() called, but structf module not initialized")
endif

if (core_num>0) then
   ! output cores:
   if (my_pe==io_pe) then
      x=g_nx;  call cwrite8(fidcore,x,1)
      x=core_num;  call cwrite8(fidcore,x,1)
      call cwrite8(fidcore,time,1)
      do n=1,core_num
         call cwrite8(fidcore,core_data(1,n),g_nx)
      enddo
      do n=1,core_num
         do i=1,3
            call cwrite8(fidcore,core_ddata(1,i,n),g_nx)
         enddo
      enddo
   endif
   core_num=0
endif

if (compute_uvw_pdfs) then
do j=1,3
do i=1,NUM_SF
   call mpisum_pdf(SF(i,j))  
enddo
enddo
call mpisum_pdf(epsilon)
endif

if (compute_passive_pdfs) then
do i=1,npassive
   call mpisum_pdf(SCALARS(i))
enddo	
endif

if (number_of_cpdf>0) then
   do i=1,number_of_cpdf
      call mpisum_pdf(cpdf(i))
   enddo
endif

if (compute_uvw_jpdfs) then
do j=1,3
do i=1,NUM_JPDF
   call mpisum_jpdf(jpdf_v(i,j))
enddo
enddo
endif



if (my_pe==io_pe) then
   if (compute_uvw_pdfs) then
      call cwrite8(fid,time,1)
      ! number of structure functions
      x=NUM_SF ; call cwrite8(fid,x,1)
      do i=1,NUM_SF
         do j=1,3
            call normalize_and_write_pdf(fid,SF(i,j),SF(i,j)%nbin)   
         enddo
      enddo
      ! number of plain pdf's
      x=1 ; call cwrite8(fid,x,1)
      call normalize_and_write_pdf(fid,epsilon,epsilon%nbin)   
   endif

   if (compute_passive_pdfs) then
      call cwrite8(fidS,time,1)
      ! number of structure functions
      x=npassive ; call cwrite8(fidS,x,1)
      do i=1,npassive
         call normalize_and_write_pdf(fidS,SCALARS(i),SCALARS(i)%nbin)   
      enddo	
   endif

   if (number_of_cpdf>0) then
      call cwrite8(fidC,time,1)
      ! number of structure functions
      x=number_of_cpdf ; call cwrite8(fidC,x,1)
      do i=1,number_of_cpdf
         call normalize_and_write_pdf(fidC,cpdf(i),cpdf(i)%nbin)   
      enddo	
   endif
   

   if (compute_uvw_jpdfs) then
      call cwrite8(fidj,time,1)
      ! number of structure functions
      x=NUM_SF ; call cwrite8(fidj,x,1)
      do i=1,NUM_JPDF
      do j=1,3
         call normalize_and_write_jpdf(fidj,jpdf_v(i,j),jpdf_v(i,j)%nbin)
      enddo
      enddo
   endif
endif


if (compute_uvw_pdfs) then
   do i=1,NUM_SF
      do j=1,3
         ! reset PDF's
         SF(i,j)%ncalls=0
         SF(i,j)%pdf=0
      enddo
   enddo
   epsilon%ncalls=0
   epsilon%pdf=0
endif

if (compute_passive_pdfs) then
   do i=1,npassive
      SCALARS(i)%ncalls=0
      SCALARS(i)%pdf=0
   enddo	
endif

if (number_of_cpdf>0) then
   do i=1,number_of_cpdf
      cpdf(i)%ncalls=0
      cpdf(i)%pdf=0
   enddo
endif




if (compute_uvw_jpdfs) then
do i=1,NUM_JPDF
do j=1,3
   ! reset JPDF's
   jpdf_v(i,j)%ncalls=0
   jpdf_v(i,j)%pdf=0
enddo
enddo
endif


end subroutine






subroutine normalize_and_write_pdf(fid,str,nbin)
use params
implicit none
integer j,i,ierr,nbin
CPOINTER :: fid
type(pdf_structure_function) :: str
real*8 x
real*8,allocatable :: pdfdata(:,:)
integer ndelta


ndelta=str%delta_num
allocate(pdfdata(-nbin:nbin,ndelta))


! the delta values
x=ndelta; call cwrite8(fid,x,1)
do i=1,ndelta
   x=delta_val(i); call cwrite8(fid,x,1)
enddo

! PDF data
call cwrite8(fid,str%pdf_bin_size,1)
x=str%nbin; call cwrite8(fid,x,1)
x=str%ncalls; call cwrite8(fid,x,1)

! normalize
pdfdata=str%pdf
x=max(str%ncalls,1)
pdfdata=pdfdata / x;
pdfdata=pdfdata/g_nx/g_ny/g_nz


! consistency check:
x=1
if (str%ncalls==0) x=0
do j=1,ndelta
   if (1e-9<abs(x - sum(pdfdata(:,j)))) then 
      print *,'ERROR in PDF normalization: '
      print *,'ndelta = ',j
      print *,'sum:     ',sum(pdfdata(:,j))
      print *,'un-norm: ',sum(str%pdf),g_nx*g_ny*real(g_nz)  ! g_n^3 can overflow integers
      print *,'ncalls:  ',str%ncalls 	
      print *,'nbin:    ',str%nbin,nbin
   endif
enddo


call cwrite8(fid,pdfdata,(2*nbin+1)*ndelta)
deallocate(pdfdata)
end subroutine








subroutine normalize_and_write_jpdf(fid,str,nbin)
use params
implicit none
integer j,i,ierr,nbin
CPOINTER :: fid
type(jpdf_structure_function) :: str
real*8 x
real*8,allocatable :: pdfdata(:,:,:)
integer ndelta


ndelta=str%delta_num
allocate(pdfdata(-nbin:nbin,-nbin:nbin,ndelta))


! the delta values
x=ndelta; call cwrite8(fid,x,1)
do i=1,ndelta
   x=delta_val(i); call cwrite8(fid,x,1)
enddo

! PDF data
call cwrite8(fid,str%pdf_bin_size,1)
x=str%nbin; call cwrite8(fid,x,1)
x=str%ncalls; call cwrite8(fid,x,1)

! normalize
pdfdata=str%pdf
x=max(str%ncalls,1)
pdfdata=pdfdata / x;
pdfdata=pdfdata/g_nx/g_ny/g_nz


! consistency check:
x=1
if (str%ncalls==0) x=0
do j=1,ndelta
   ASSERT("n_and_w(): bad normalization",1e-9>abs(x - sum(pdfdata(:,:,j)))  )
enddo


call cwrite8(fid,pdfdata,(2*nbin+1)*(2*nbin+1)*ndelta)
deallocate(pdfdata)
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
real*8,allocatable :: pdfdata(:,:)
integer :: bin,ierr,n,ndelta,ncalls2


#ifdef USE_MPI
! find maximum value of bin  MPI max
call mpi_allreduce(str%nbin,bin,1,MPI_INTEGER,MPI_MAX,comm_3d,ierr)
call mpi_allreduce(str%ncalls,ncalls2,1,MPI_INTEGER,MPI_MAX,comm_3d,ierr)
str%ncalls=ncalls2  ! in case io_pe has ncalls=0

! resize all str's to size bin
call resize_pdf(str,bin)

!
! MPI sum into pdfdata
!
ndelta=str%delta_num
n=(2*bin+1)*ndelta
if (my_pe == io_pe) then
   allocate(pdfdata(-bin:bin,ndelta))
   call mpi_reduce(str%pdf,pdfdata,n,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   str%pdf=pdfdata
   deallocate(pdfdata)
else
   call mpi_reduce(str%pdf,0,n,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
endif

#endif

end subroutine











subroutine mpisum_jpdf(str)
!
! sum PDF's over all CPUs into one PDF on cpu = io_pe
!
use params
use mpi
implicit none
type(jpdf_structure_function) :: str


! local variables
real*8,allocatable :: pdfdata(:,:,:)
integer :: bin,ierr,n,ndelta,ncalls2


#ifdef USE_MPI
! find maximum value of bin  MPI max
call mpi_allreduce(str%nbin,bin,1,MPI_INTEGER,MPI_MAX,comm_3d,ierr)
call mpi_allreduce(str%ncalls,ncalls2,1,MPI_INTEGER,MPI_MAX,comm_3d,ierr)
str%ncalls=ncalls2  ! in case io_pe has ncalls=0

! resize all str's to size bin
call resize_jpdf(str,bin)

!
! MPI sum into pdfdata
!
ndelta=str%delta_num
n=(2*bin+1)*(2*bin+1)*ndelta
if (my_pe == io_pe) then
   allocate(pdfdata(-bin:bin,-bin:bin,ndelta))
   call mpi_reduce(str%pdf,pdfdata,n,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   str%pdf=pdfdata
   deallocate(pdfdata)
else
   call mpi_reduce(str%pdf,0,n,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
endif

#endif

end subroutine



subroutine compute_cores(pt,n1,n1d,n2,n2d,n3,n3d,core,n)
use params
implicit none
integer :: n1,n1d,n2,n2d,n3,n3d,n
real*8 :: pt(n1d,n2d,n3d)
real*8 :: core(n1)
if (n1>g_nx) then
   call abortdns("Error: compute_core() called on non traspose-x data?")
endif
core(1:n1)=pt(1:n1,1,1)
core_num=max(core_num,n)
end subroutine




subroutine compute_pdf(u,v,w,n1,n1d,n2,n2d,n3,n3d,str,ncomp)
!
! compute a pdf_structure function along the first dimension of Q
! for all the values of delta given by delta_val(:)
!
!
! ncomp=1,2 or 3  we are computing them in the x,y or z direction
!
use params
implicit none
integer :: n1,n1d,n2,n2d,n3,n3d,n,ncomp
real*8 :: u(n1d,n2d,n3d)
real*8 :: v(n1d,n2d,n3d)
real*8 :: w(n1d,n2d,n3d)
type(pdf_structure_function) :: str(NUM_SF)

! local variables
real*8  :: del,delv(3)
integer :: bin,idel,i,j,k,i2,nsf,ndelta


if (structf_init==0) then
   structf_init=1
   call init_pdf_module()
endif


ndelta=str(1)%delta_num
do j=2,NUM_SF
ASSERT("ndelta must be the same for all U structure functions",ndelta==str(j)%delta_num)
enddo



do k=1,n3
   do j=1,n2
      do idel=1,ndelta
         if (delta_val(idel) < n1/2) then
            do i=1,n1

               ! compute structure functions for U,V,W 
               do n=1,3
                  i2 = i + delta_val(idel)
                  if (i2>n1) i2=i2-n1
                  if (n==1) then
                     del = u(i2,j,k)-u(i,j,k)
                  else if (n==2) then
                     del = v(i2,j,k)-v(i,j,k)
                  else
                     del = w(i2,j,k)-w(i,j,k)
                  endif

                  delv(n)=del
                  del = del/str(n)%pdf_bin_size
                  bin = nint(del)

                  ! increase the size of our PDF function
                  if (abs(bin)>str(n)%nbin) call resize_pdf(str(n),abs(bin)+10) 
                  if (bin>pdf_max_bin) bin=pdf_max_bin
                  if (bin<-pdf_max_bin) bin=-pdf_max_bin
                  str(n)%pdf(bin,idel)=str(n)%pdf(bin,idel)+1
               enddo

               if (NUM_SF>=6) then
               ! compute structure functions for U U**2, U V**2 and U W**2
               ! but U * U**2 was already computed above, so replace
               ! this one with U (U**2 + V**2 + W**2)
               do n=1,3
                  nsf=n+3
                  if (n==ncomp) then
                     del=delv(ncomp)*(delv(1)**2+delv(2)**2+delv(3)**2)
                  else
                     del=delv(ncomp)*delv(n)**2
                  endif

                  if (del>=0) then
                     del=del**one_third
                  else
                     del=-((-del)**one_third)
                  endif
                  del = del/str(nsf)%pdf_bin_size
                  bin=nint(del)

                  if (abs(bin)>str(nsf)%nbin) call resize_pdf(str(nsf),abs(bin)+10) 
                  if (bin>pdf_max_bin) bin=pdf_max_bin
                  if (bin<-pdf_max_bin) bin=-pdf_max_bin
                  str(nsf)%pdf(bin,idel)=str(nsf)%pdf(bin,idel)+1
               enddo
               endif


               ! compute structure functions for U**2 V**2 and U**2  W**2
               if (NUM_SF>=8) then
               nsf=6
               do n=1,3
                  if (n==ncomp) then
                     ! skip diagonal term
                  else
                     nsf=nsf+1
                     ! (delv(ncomp)**2 delv(n)**2  ) ** 1/4
                     del=sqrt(abs(delv(ncomp)*delv(n)))
                     del = del/str(nsf)%pdf_bin_size
                     bin=nint(del)
                     
                     if (abs(bin)>str(nsf)%nbin) call resize_pdf(str(nsf),abs(bin)+10) 
                     if (bin>pdf_max_bin) bin=pdf_max_bin
                     if (bin<-pdf_max_bin) bin=-pdf_max_bin
                     str(nsf)%pdf(bin,idel)=str(nsf)%pdf(bin,idel)+1
                  endif
               enddo
               endif



            enddo
         endif
      enddo
   enddo
enddo

do n=1,NUM_SF
str(n)%ncalls=str(n)%ncalls+1
enddo

end subroutine




subroutine compute_jpdf(u,v,w,n1,n1d,n2,n2d,n3,n3d,str,ncomp)
!
! compute a pdf_structure function along the first dimension of Q
! for all the values of delta given by delta_val(:)
!
!
! ncomp=1,2 or 3  we are computing them in the x,y or z direction
!
use params
implicit none
integer :: n1,n1d,n2,n2d,n3,n3d,n,ncomp
real*8 :: u(n1d,n2d,n3d)
real*8 :: v(n1d,n2d,n3d)
real*8 :: w(n1d,n2d,n3d)
type(jpdf_structure_function) :: str(NUM_JPDF)

! local variables
real*8  :: del1,del2,delv(3)
integer :: bin1,bin2,bin,idel,i,j,k,i2,nsf,ndelta


if (structf_init==0) then
   structf_init=1
   call init_pdf_module()
endif


ndelta=str(1)%delta_num
do j=2,NUM_JPDF
ASSERT("ndelta must be the same for all U structure functions",ndelta==str(j)%delta_num)
enddo



do k=1,n3
   do j=1,n2
      do idel=1,ndelta
         if (delta_val(idel) < n1/2) then
            do i=1,n1
               ! compute structure functions for U,V,W 

               i2 = i + delta_val(idel)
               if (i2>n1) i2=i2-n1
               do n=1,3
                  if (n==1) then
                     del1=u(i2,j,k)-u(i,j,k)
                     del2=v(i2,j,k)-v(i,j,k)
                  endif
                  if (n==2) then
                     del1=u(i2,j,k)-u(i,j,k)
                     del2=w(i2,j,k)-v(i,j,k)
                  endif
                  if (n==3) then
                     del1=v(i2,j,k)-u(i,j,k)
                     del2=w(i2,j,k)-v(i,j,k)
                  endif

                  del1 = del1/str(n)%pdf_bin_size
                  del2 = del2/str(n)%pdf_bin_size
                  bin1 = nint(del1)
                  bin2 = nint(del2)
                  bin=max(abs(bin1),abs(bin2))

                  ! increase the size of our PDF function
                  if (bin>str(n)%nbin) call resize_jpdf(str(n),abs(bin)+10) 
                  if (bin1>jpdf_max_bin) bin1=jpdf_max_bin
                  if (bin1<-jpdf_max_bin) bin1=-jpdf_max_bin
                  if (bin2>jpdf_max_bin) bin2=jpdf_max_bin
                  if (bin2<-jpdf_max_bin) bin2=-jpdf_max_bin

                  str(n)%pdf(bin1,bin2,idel)=str(n)%pdf(bin1,bin2,idel)+1
               enddo
            enddo
         endif
      enddo
   enddo
enddo

do n=1,NUM_JPDF
str(n)%ncalls=str(n)%ncalls+1
enddo

end subroutine





subroutine compute_pdf_scalar(ux,pdfdata,binsize)
!
! compute a pdf of a scalar, 1 point corrolation
! 
!
use params
implicit none
real*8 :: ux(nx,ny,nz)
type(pdf_structure_function) ::  pdfdata
real*8,optional :: binsize

! local variables
real*8  :: del
integer :: bin,idel,i,j,k,i2,nsf

if (structf_init==0) then
   structf_init=1
   call init_pdf_module()
endif

if (present(binsize)) then
   print *,'here a',pdfdata%ncalls
   if (pdfdata%ncalls/=0) then
      call abortdns("compute_pdf_scalar(): ERROR:  cant change PDF binsize unless ncalls=0")
   endif
   pdfdata%pdf_bin_size=binsize
endif

print *,'done'
if (pdfdata%pdf_bin_size == 0) then
   call abortdns("compute_pdf_scalar(): ERROR:  binsize not initialized")
endif



do k=nz1,nz2
   do j=ny1,ny2
      do i=nx1,nx2
         ! compute structure functions for U,V,W 
         print *,i,j,k
         del=ux(i,j,k)
         del = del/pdfdata%pdf_bin_size
         bin = nint(del)

         ! increase the size of our PDF function
         if (abs(bin)>pdfdata%nbin) call resize_pdf(pdfdata,abs(bin)+10) 
         if (bin>pdf_max_bin) bin=pdf_max_bin
         if (bin<-pdf_max_bin) bin=-pdf_max_bin
         pdfdata%pdf(bin,1)=pdfdata%pdf(bin,1)+1
      enddo
   enddo
enddo

pdfdata%ncalls=pdfdata%ncalls+1
end subroutine








subroutine compute_all_pdfs(Q,gradu)
!
!  Compute the velocity increment PDFs 
!  compute the epsilon PDF
!  compute the passive scalar PDFs
!
use params
use fft_interface
use transpose
implicit none
integer :: ns
real*8 Q(nx,ny,nz,n_var)    
real*8 gradu(nx,ny,nz,3)    

!local
integer n1,n1d,n2,n2d,n3,n3d,ierr
integer i,j,k,n,m1,m2
real*8 :: vor(3),uij,uji
real*8 :: dummy(1),ensave
real*8 :: tmx1,tmx2


if (compute_uvw_pdfs) then

call wallclock(tmx1)
call print_message("computing x direction pdfs...")
do n=1,3
   call transpose_to_x(Q(1,1,1,n),gradu(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
   call compute_cores(gradu(1,1,1,n),n1,n1d,n2,n2d,n3,n3d,core_data(1,n),n)
enddo
call compute_pdf(gradu(1,1,1,1),gradu(1,1,1,2),gradu(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,SF(1,1),1)
if (compute_uvw_jpdfs) then
   call compute_jpdf(gradu(1,1,1,1),gradu(1,1,1,2),gradu(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,jpdf_v(1,1),1)
endif


call print_message("computing y direction pdfs...")
do n=1,3
   call transpose_to_y(Q(1,1,1,n),gradu(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo
call compute_pdf(gradu(1,1,1,1),gradu(1,1,1,2),gradu(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,SF(1,2),2)
if (compute_uvw_jpdfs) then
   call compute_jpdf(gradu(1,1,1,1),gradu(1,1,1,2),gradu(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,jpdf_v(1,2),2)
endif


call print_message("computing z direction pdfs...")
do n=1,3
   call transpose_to_z(Q(1,1,1,n),gradu(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo
call compute_pdf(gradu(1,1,1,1),gradu(1,1,1,2),gradu(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,SF(1,3),3)
if (compute_uvw_jpdfs) then
   call compute_jpdf(gradu(1,1,1,1),gradu(1,1,1,2),gradu(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,jpdf_v(1,3),3)
endif

call print_message("computing eps pdfs...")
gradu=0
do n=1,3
   ! u_x, u_y, u_z
   call der(Q(1,1,1,1),gradu(1,1,1,3),dummy,gradu(1,1,1,2),DX_ONLY,n)
   gradu(:,:,:,1)=gradu(:,:,:,1)+gradu(:,:,:,3)**2	

   call transpose_to_x(gradu(1,1,1,3),gradu(1,1,1,2),n1,n1d,n2,n2d,n3,n3d)
   call compute_cores(gradu(1,1,1,2),n1,n1d,n2,n2d,n3,n3d,core_ddata(1,n,1),1)



   ! v_x, v_y, v_z
   call der(Q(1,1,1,2),gradu(1,1,1,3),dummy,gradu(1,1,1,2),DX_ONLY,n)
   gradu(:,:,:,1)=gradu(:,:,:,1)+gradu(:,:,:,3)**2	

   call transpose_to_x(gradu(1,1,1,3),gradu(1,1,1,2),n1,n1d,n2,n2d,n3,n3d)
   call compute_cores(gradu(1,1,1,2),n1,n1d,n2,n2d,n3,n3d,core_ddata(1,n,2),2)


   ! w_x, w_y, w_z
   call der(Q(1,1,1,3),gradu(1,1,1,3),dummy,gradu(1,1,1,2),DX_ONLY,n)
   gradu(:,:,:,1)=gradu(:,:,:,1)+gradu(:,:,:,3)**2	

   call transpose_to_x(gradu(1,1,1,3),gradu(1,1,1,2),n1,n1d,n2,n2d,n3,n3d)
   call compute_cores(gradu(1,1,1,2),n1,n1d,n2,n2d,n3,n3d,core_ddata(1,n,3),3)

enddo
gradu(:,:,:,1)=mu*gradu(:,:,:,1); 
gradu(:,:,:,1)=gradu(:,:,:,1)**one_third
call compute_pdf_scalar(gradu,epsilon)
endif




if (compute_passive_pdfs) then
   call print_message("computing passive scalar pdfs...")
   do i=np1,np2
      call compute_pdf_scalar(Q(1,1,1,i),SCALARS(i-np1+1)) 

      call transpose_to_x(Q(1,1,1,i),gradu,n1,n1d,n2,n2d,n3,n3d)
      call compute_cores(gradu,n1,n1d,n2,n2d,n3,n3d,core_data(1,i),i)

      do n=1,3
         call der(Q(1,1,1,i),gradu,dummy,gradu(1,1,1,2),DX_ONLY,n)
         call transpose_to_x(gradu,gradu(1,1,1,2),n1,n1d,n2,n2d,n3,n3d)
         call compute_cores(gradu(1,1,1,2),n1,n1d,n2,n2d,n3,n3d,core_ddata(1,n,i),i)
      enddo
   enddo
endif
call print_message("done with pdfs.")




call wallclock(tmx2)
tims(12)=tims(12)+(tmx2-tmx1)

end subroutine






end module






