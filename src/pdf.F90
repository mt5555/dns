#include "macros.h"
#undef COMP_JPDF

module structf
implicit none

! compute structure functions every "struct_n*" time steps
! x and y structure functions can be computed with no additional
! transforms at most every other timestep. (struct_nx,y >= 2)
!
! struct_nz controls z structure functions, and the epsilon 1 point
! structure functions.  This requires 2 extra z-tranpose plus
! an x and y tranpose.   
!
! structure functions are output every 'diag_dt' time interval.  
! Set struct_n* is very large to only compute structure functions
! once per diag_dt time interval.  
! 
! 
!
real*8 :: uscale=.01             ! bin size for vel increment
real*8 :: epsscale=.01           ! bin size for epsilon increment


integer :: struct_nx=0 
integer :: struct_ny=0          ! every "struct_n*" time steps
integer :: struct_nz=0          !
integer :: countx=-1,county=-1,countz=-1


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
integer,parameter :: NUM_SF=8
type(pdf_structure_function) ::  SF(NUM_SF,3),epsilon

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
   call abort("structf init: j > delta_num_max")
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


do i=1,NUM_SF
   call init_pdf(SF(i,1),100,uscale,numx)
   call init_pdf(SF(i,2),100,uscale,numy)
   call init_pdf(SF(i,3),100,uscale,numz)
enddo
call init_pdf(epsilon,100,epsscale,1)

#ifdef COMP_JPDF
do i=1,NUM_JPDF
   call init_jpdf(jpdf_v(i,1),100,.1d0, min(numx,numy))
   call init_jpdf(jpdf_v(i,2),100,.1d0, min(numx,numz))
   call init_jpdf(jpdf_v(i,3),100,.1d0, min(numy,numz))
enddo
#endif
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
if (bin<n) call abort("resize_pdf: attempting to shrink pdf");
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
if (bin<n) call abort("resize_pdf: attempting to shrink pdf");
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








subroutine output_pdf(time,fid,fidj,fidS)
use params
implicit none
real*8 time
CPOINTER fid,fidj,fidS

integer i,j,ierr,ndelta,numm
character(len=80) message
real*8 x

if (structf_init==0) then
   call abort("Error: output_pdf() called, but structf module not initialized")
endif

! reset structure function counters
countx=-1
county=-1
countz=-1


do j=1,3
do i=1,NUM_SF
   call mpisum_pdf(SF(i,j))  
enddo
enddo
call mpisum_pdf(epsilon)

#ifdef COMP_JPDF
do j=1,3
do i=1,NUM_JPDF
   call mpisum_jpdf(jpdf_v(i,j))
enddo
enddo
#endif



if (my_pe==io_pe) then
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

#ifdef COMP_JPDF
   call cwrite8(fidj,time,1)
   ! number of structure functions
   x=NUM_SF ; call cwrite8(fidj,x,1)
   do i=1,NUM_JPDF
   do j=1,3
      call normalize_and_write_jpdf(fidj,jpdf_v(i,j),jpdf_v(i,j)%nbin)
   enddo
   enddo
#endif

endif



do i=1,NUM_SF
do j=1,3
   ! reset PDF's
   SF(i,j)%ncalls=0
   SF(i,j)%pdf=0
enddo
enddo
epsilon%ncalls=0
epsilon%pdf=0

#ifdef COMP_JPDF
do i=1,NUM_JPDF
do j=1,3
   ! reset JPDF's
   jpdf_v(i,j)%ncalls=0
   jpdf_v(i,j)%pdf=0
enddo
enddo
#endif

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
   ASSERT("n_and_w(): bad normalization",1e-9>abs(x - sum(pdfdata(:,j)))  )
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
integer :: bin,ierr,n,ndelta,ncalls


#ifdef USE_MPI
! find maximum value of bin  MPI max
call MPI_allreduce(str%nbin,bin,1,MPI_INTEGER,MPI_MAX,comm_3d,ierr)
call MPI_allreduce(str%ncalls,ncalls,1,MPI_INTEGER,MPI_MAX,comm_3d,ierr)
str%ncalls=ncalls  ! in case io_pe has ncalls=0

! resize all str's to size bin
call resize_pdf(str,bin)

!
! MPI sum into pdfdata
!
ndelta=str%delta_num
n=(2*bin+1)*ndelta
if (my_pe == io_pe) then
   allocate(pdfdata(-bin:bin,ndelta))
   call MPI_reduce(str%pdf,pdfdata,n,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   str%pdf=pdfdata
   deallocate(pdfdata)
else
   call MPI_reduce(str%pdf,0,n,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
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
integer :: bin,ierr,n,ndelta,ncalls


#ifdef USE_MPI
! find maximum value of bin  MPI max
call MPI_allreduce(str%nbin,bin,1,MPI_INTEGER,MPI_MAX,comm_3d,ierr)
call MPI_allreduce(str%ncalls,ncalls,1,MPI_INTEGER,MPI_MAX,comm_3d,ierr)
str%ncalls=ncalls  ! in case io_pe has ncalls=0

! resize all str's to size bin
call resize_jpdf(str,bin)

!
! MPI sum into pdfdata
!
ndelta=str%delta_num
n=(2*bin+1)*(2*bin+1)*ndelta
if (my_pe == io_pe) then
   allocate(pdfdata(-bin:bin,-bin:bin,ndelta))
   call MPI_reduce(str%pdf,pdfdata,n,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
   str%pdf=pdfdata
   deallocate(pdfdata)
else
   call MPI_reduce(str%pdf,0,n,MPI_REAL8,MPI_SUM,io_pe,comm_3d,ierr)
endif

#endif

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


               ! compute structure functions for U**2 V**2 and U**2  W**2
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





subroutine compute_pdf_epsilon(ux)
!
! compute a pdf of epsilon.  < epsilon > = KE diffusion. 
! ux = norm(grad(u))**2
!
use params
implicit none
real*8 :: ux(nx,ny,nz)

! local variables
real*8  :: del
integer :: bin,idel,i,j,k,i2,nsf

if (structf_init==0) then
   structf_init=1
   call init_pdf_module()
endif



do k=nz1,nz2
   do j=ny1,ny2
      do i=nx1,nx2
         ! compute structure functions for U,V,W 
         del=ux(i,j,k)**one_third
         del = del/epsilon%pdf_bin_size
         bin = nint(del)

         ! increase the size of our PDF function
         if (abs(bin)>epsilon%nbin) call resize_pdf(epsilon,abs(bin)+10) 
         if (bin>pdf_max_bin) bin=pdf_max_bin
         if (bin<-pdf_max_bin) bin=-pdf_max_bin
         epsilon%pdf(bin,1)=epsilon%pdf(bin,1)+1
      enddo
   enddo
enddo

epsilon%ncalls=epsilon%ncalls+1
end subroutine



#if 0

subroutine z_ifft3d_str(fin,f,w1,Qt,works,work)
!
!  compute inverse fft 3d of fin, return in f
!  fin and f can overlap in memory
!  Also compute structure functions.
!
!  This routine can compute PDF's of structure functions in the x direction
!  or y direction (but not both) and in the z direction.
! 
!  compx   compute PDF in x direction of U
!  compy   compute PDF in y direction of U
!  compz   compute PDF in z direction of U and epsilon = gradU**2
!
! only one of compx or compy can be computed per call.
! compx or compy can be computed w/o any additional transposes.
! compz requires an additional z-transpose.
!
use params
use fft_interface
use transpose
implicit none
real*8 fin(g_nz2,nslabx,ny_2dz,3)  ! input
real*8 f(nx,ny,nz,3)               ! output
real*8 w1(nx,ny,nz,3)
real*8 Qt(nx,ny,nz,3)    

! overlapped in memory:
real*8 work(nx,ny,nz)
real*8 works(g_nz2,nslabx,ny_2dz)


logical :: compx,compy,compz


!local
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k,n
real*8 dummy(1)
real*8 :: tmx1,tmx2

call wallclock(tmx1)

if (struct_nx>0) countx=mod(countx+1,struct_nx)  
if (struct_ny>0) county=mod(county+1,struct_ny)  
if (struct_nz>0) countz=mod(countz+1,struct_nz)  
compx=(countx==0)
compy=(county==1)
compz=(countz==0)




n1=g_nz
n1d=g_nz2   	
n2=nslabx
n2d=nslabx
n3=ny_2dz
n3d=ny_2dz

do n=1,3
   works=fin(:,:,:,n)
   call ifft1(works,n1,n1d,n2,n2d,n3,n3d)
   call transpose_from_z(works,f(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo

if (compx) then
   do n=1,3
      call transpose_to_y(f(1,1,1,n),Qt(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      call ifft1(Qt(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      call transpose_from_y(Qt(1,1,1,n),f(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      
      call transpose_to_x(f(1,1,1,n),Qt(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      call ifft1(Qt(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      call transpose_from_x(Qt(1,1,1,n),f(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
   enddo
   call compute_pdf(Qt(1,1,1,1),Qt(1,1,1,2),Qt(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,SF(1,1),1)
   
else ! compy, or default (compute no structure functions)
   do n=1,3
      call transpose_to_x(f(1,1,1,n),Qt(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      call ifft1(Qt(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      call transpose_from_x(Qt(1,1,1,n),f(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      
      call transpose_to_y(f(1,1,1,n),Qt(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      call ifft1(Qt(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
      call transpose_from_y(Qt(1,1,1,n),f(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
   enddo
   if (compy) call compute_pdf(Qt(1,1,1,1),Qt(1,1,1,2),Qt(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,SF(1,2),2)
endif

if (compz) then

   ! compute:
   ! 1. epsilon PDF calculations
   ! 2. z-direction Q PDF
   !
   !   Qt = x or y transpose of (u,v,w)
   !
   work=0
   do n=1,3
      if (compx) then
         ! Qt = x transpoe of (u,v,w), computed above
         ! compute x derivative
         call fft_derivatives(Qt(1,1,1,n),dummy,1,n1,n1d,n2,n2d,n3,n3d)
         call transpose_from_x(Qt(1,1,1,n),w1,n1,n1d,n2,n2d,n3,n3d)
         work=work+w1(:,:,:,1)**2

         ! compute y derivative
         call der(f(1,1,1,n),w1,dummy,w1(1,1,1,2),DX_ONLY,2)
         work=work+w1(:,:,:,1)**2

      else
         ! Qt = y transpoe of (u,v,w), computed above
         call fft_derivatives(Qt(1,1,1,n),dummy,1,n1,n1d,n2,n2d,n3,n3d)
         call transpose_from_y(Qt(1,1,1,n),w1,n1,n1d,n2,n2d,n3,n3d)
         work=work+w1(:,:,:,1)**2

         ! compute x derivative
         call der(f(1,1,1,n),w1,dummy,w1(1,1,1,2),DX_ONLY,1)
         work=work+w1(:,:,:,1)**2
      endif
   enddo

   ! do the Z component
   do n=1,3
      call transpose_to_z(f(1,1,1,n),Qt(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
   enddo
   ! Qt = z transpose of (u,v,w)
   ! compute regular structure functions:
    call compute_pdf(Qt(1,1,1,1),Qt(1,1,1,2),Qt(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,SF(1,3),3)

   ! add in epsilon component
   do n=1,3
      call fft_derivatives(Qt(1,1,1,n),dummy,1,n1,n1d,n2,n2d,n3,n3d)
      call transpose_from_z(Qt(1,1,1,n),w1,n1,n1d,n2,n2d,n3,n3d)
      work=work+w1(:,:,:,1)**2
   enddo
   work=mu*work
   call compute_pdf_epsilon(work)

endif

call wallclock(tmx2)
tims(12)=tims(12)+(tmx2-tmx1)

end subroutine

#endif







subroutine compute_all_pdfs(Q,gradu)
!
!
use params
use fft_interface
use transpose
implicit none
integer :: ns
real*8 Q(nx,ny,nz,3)    
real*8 gradu(nx,ny,nz,3)    

!local
integer n1,n1d,n2,n2d,n3,n3d,ierr
integer i,j,k,n,m1,m2
real*8 :: vor(3),uij,uji
real*8 :: dummy(1),ensave
real*8 :: tmx1,tmx2



call wallclock(tmx1)
call print_message("computing u pdfs...")
do n=1,3
   call transpose_to_x(Q(1,1,1,n),gradu(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo
call compute_pdf(gradu(1,1,1,1),gradu(1,1,1,2),gradu(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,SF(1,1),1)
#ifdef COMP_JPDF
call compute_jpdf(gradu(1,1,1,1),gradu(1,1,1,2),gradu(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,jpdf_v(1,1),1)
#endif

call print_message("computing v pdfs...")
do n=1,3
   call transpose_to_y(Q(1,1,1,n),gradu(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo
call compute_pdf(gradu(1,1,1,1),gradu(1,1,1,2),gradu(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,SF(1,2),2)
#ifdef COMP_JPDF
call compute_jpdf(gradu(1,1,1,1),gradu(1,1,1,2),gradu(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,jpdf_v(1,2),2)
#endif

call print_message("computing w pdfs...")
do n=1,3
   call transpose_to_z(Q(1,1,1,n),gradu(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo
call compute_pdf(gradu(1,1,1,1),gradu(1,1,1,2),gradu(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,SF(1,3),3)
#ifdef COMP_JPDF
call compute_jpdf(gradu(1,1,1,1),gradu(1,1,1,2),gradu(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,jpdf_v(1,3),3)
#endif

call print_message("computing eps pdfs...")
gradu=0
do n=1,3
   ! u_x, u_y, u_z
   call der(Q(1,1,1,1),gradu(1,1,1,3),dummy,gradu(1,1,1,2),DX_ONLY,n)
   gradu(:,:,:,1)=gradu(:,:,:,1)+gradu(:,:,:,3)**2	

   ! v_x, v_y, v_z
   call der(Q(1,1,1,2),gradu(1,1,1,3),dummy,gradu(1,1,1,2),DX_ONLY,n)
   gradu(:,:,:,1)=gradu(:,:,:,1)+gradu(:,:,:,3)**2	

   ! w_x, w_y, w_z
   call der(Q(1,1,1,3),gradu(1,1,1,3),dummy,gradu(1,1,1,2),DX_ONLY,n)
   gradu(:,:,:,1)=gradu(:,:,:,1)+gradu(:,:,:,3)**2	
enddo
gradu=mu*gradu
call compute_pdf_epsilon(gradu)
call print_message("done with pdfs.")




call wallclock(tmx2)
tims(12)=tims(12)+(tmx2-tmx1)

end subroutine






end module















