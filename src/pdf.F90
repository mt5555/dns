#include "macros.h"

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
integer :: struct_nx=0 
integer :: struct_ny=0          ! every "struct_n*" time steps
integer :: struct_nz=0          !
integer :: countx=-1,county=-1,countz=-1
integer :: compute_struct=0


integer           :: structf_init=0
integer,parameter :: delta_num_max=16
integer           :: delta_val(delta_num_max)

! max range:  10 ... 10
integer,parameter :: pdf_max_bin=2000





type pdf_structure_function
!
! pdf%data(-bin:bin,delta_num,npdfs)
!
!SGI f90 does not allow allocatable arrays in a struct
real*8,pointer      :: pdf(:,:)
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
!
integer,parameter :: NUM_SF=6
type(pdf_structure_function) ::  SF(NUM_SF,3),epsilon

integer :: overflow=0    ! count the number of overflow messages
real*8  :: one_third = (1d0/3d0)

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

integer idel,i
integer :: numx=0,numy=0,numz=0

delta_val=99999
delta_val(1)=1
delta_val(2)=2
delta_val(3)=3
delta_val(4)=4
delta_val(5)=6
delta_val(6)=8
delta_val(7)=11
delta_val(8)=16
delta_val(9)=23
delta_val(10)=32
delta_val(11)=64
delta_val(12)=128
delta_val(13)=256
delta_val(14)=512
delta_val(15)=1024
delta_val(16)=2048
ASSERT("delta_num_max to small. ",16==delta_num_max)


do idel=1,delta_num_max
   if (delta_val(idel) < g_nx/2) then
      numx=idel
   endif
   if (delta_val(idel) < g_ny/2) then
      numy=idel
   endif
   if (delta_val(idel) < g_nz/2) then
      numz=idel
   endif
enddo



do i=1,NUM_SF
   call init_pdf(SF(i,1),100,.01d0,numx)
   call init_pdf(SF(i,2),100,.01d0,numy)
   call init_pdf(SF(i,3),100,.01d0,numz)
enddo
call init_pdf(epsilon,100,.01d0,1)
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








subroutine outputSF(time,fid)
use params
implicit none
real*8 time
CPOINTER fid

integer i,j,ierr,ndelta
character(len=80) message
real*8 x

if (structf_init==0) then
   call abort("Error: outputSF() called, but structf module not initialized")
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
            enddo
         endif
      enddo
   enddo
enddo

do n=1,NUM_SF
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
real*8 fin(g_nz2,nslabx,ny_2dz,n_var)  ! input
real*8 f(nx,ny,nz,n_var)               ! output
real*8 w1(nx,ny,nz,n_var)
real*8 Qt(nx,ny,nz,n_var)    

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









subroutine compute_all_pdfs(Q,gradu,gradv,gradw,work,scalars,ns)
!
!
use params
use fft_interface
use transpose
implicit none
integer :: ns
real*8 :: scalars(ns)
real*8 Q(nx,ny,nz,n_var)    
real*8 work(nx,ny,nz)
real*8 gradu(nx,ny,nz,n_var)    
real*8 gradv(nx,ny,nz,n_var)    
real*8 gradw(nx,ny,nz,n_var)    

!local
real*8 :: scalars2(ns)
integer n1,n1d,n2,n2d,n3,n3d,ierr
integer i,j,k,n,m1,m2
real*8 :: vor(3),Sw(3),wS(3),Sww,ux2(3),ux3(3),ux4(3),uij,uji,u2(3),S2sum
real*8 :: S(3,3)
real*8 dummy(1)
real*8 :: tmx1,tmx2


call wallclock(tmx1)

do n=1,3
   call transpose_to_x(Q(1,1,1,n),gradu(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo
call compute_pdf(gradu(1,1,1,1),gradu(1,1,1,2),gradu(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,SF(1,1),1)

do n=1,3
   call transpose_to_y(Q(1,1,1,n),gradu(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo
call compute_pdf(gradu(1,1,1,1),gradu(1,1,1,2),gradu(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,SF(1,2),2)

do n=1,3
   call transpose_to_z(Q(1,1,1,n),gradu(1,1,1,n),n1,n1d,n2,n2d,n3,n3d)
enddo
call compute_pdf(gradu(1,1,1,1),gradu(1,1,1,2),gradu(1,1,1,3),n1,n1d,n2,n2d,n3,n3d,SF(1,3),3)


do n=1,3
   call der(Q(1,1,1,1),gradu(1,1,1,n),dummy,work,DX_ONLY,n)
   call der(Q(1,1,1,2),gradv(1,1,1,n),dummy,work,DX_ONLY,n)
   call der(Q(1,1,1,3),gradw(1,1,1,n),dummy,work,DX_ONLY,n)
enddo
work=0
do n=1,3
   work=work+gradu(:,:,:,n)**2
   work=work+gradv(:,:,:,n)**2
   work=work+gradw(:,:,:,n)**2
enddo
work=mu*work
call compute_pdf_epsilon(work)


#if 0
! cj structure functions
! 
! note: fix S2 below.  it is wrong
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   vor(1)=gradw(i,j,k,2)-gradv(i,j,k,3)
   vor(2)=gradu(i,j,k,3)-gradw(i,j,k,1)
   vor(3)=gradv(i,j,k,1)-gradu(i,j,k,2)
   v2(i,j,k)=vor(1)**2+vor(2)**2+vor(3)**2
   S(m1,m2)=      
enddo
enddo
enddo

! compute integrals of:
! <v2(x),v2(x+r)>
! <S2(x),S2(x+r)>
! <S2(x),v2(x+r)>
!
#endif



! scalars
S2sum=0
Sww=0
ux2=0
ux3=0
ux4=0
u2=0


do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   do n=1,3
      u2(n)=u2(n)+Q(i,j,k,n)**2
   enddo

   vor(1)=gradw(i,j,k,2)-gradv(i,j,k,3)
   vor(2)=gradu(i,j,k,3)-gradw(i,j,k,1)
   vor(3)=gradv(i,j,k,1)-gradu(i,j,k,2)

   ! compute Sw = Sij*wj
   Sw=0
   !wS=0
   do m1=1,3
      do m2=1,3
         if (m1==1) uij=gradu(i,j,k,m2)
         if (m1==2) uij=gradv(i,j,k,m2)
         if (m1==3) uij=gradw(i,j,k,m2)
         if (m2==1) uji=gradu(i,j,k,m1)
         if (m2==2) uji=gradv(i,j,k,m1)
         if (m2==3) uji=gradw(i,j,k,m1)
         ! S(m1,m2) = .5*(uij_uji)
         Sw(m1)=Sw(m1)+.5*(uij+uji)*vor(m2)
         !wS(m2)=wS(m2)+.5*(uij+uji)*vor(m1)
      enddo
   enddo
   ! compute Sww = wi*(Sij*wj)
   Sww=Sww+Sw(1)*vor(1)+Sw(2)*vor(2)+Sw(3)*vor(3)

   ! if we use gradu(i,j,k,1)**3, do we preserve the sign?  
   ! lets not put f90 to that test!
   uij=gradu(i,j,k,1)**2
   ux2(1)=ux2(1)+uij
   ux3(1)=ux3(1)+uij*gradu(i,j,k,1)
   ux4(1)=ux4(1)+uij*uij

   uij=gradv(i,j,k,2)**2
   ux2(2)=ux2(2)+uij
   ux3(2)=ux3(2)+uij*gradv(i,j,k,2)
   ux4(2)=ux4(2)+uij*uij

   uij=gradw(i,j,k,3)**2
   ux2(3)=ux2(3)+uij
   ux3(3)=ux3(3)+uij*gradw(i,j,k,3)
   ux4(3)=ux4(3)+uij*uij

enddo
enddo
enddo

S2sum=S2sum/g_nx/g_ny/g_nz
Sww=Sww/g_nx/g_ny/g_nz
ux2=ux2/g_nx/g_ny/g_nz
ux3=ux3/g_nx/g_ny/g_nz
ux4=ux4/g_nx/g_ny/g_nz
u2=u2/g_nx/g_ny/g_nz

ASSERT("compute_all_pdfs: ns too small ",ns>=14)
do n=1,3
scalars(n)=ux2(n)
scalars(n+3)=ux3(n)
scalars(n+6)=ux4(n)
enddo
scalars(10)=Sww
do n=1,3
scalars(10+n)=u2(n)
enddo
scalars(14)=S2sum


#ifdef USE_MPI
   scalars2=scalars
   call MPI_allreduce(scalars2,scalars,ns,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif


call wallclock(tmx2)
tims(12)=tims(12)+(tmx2-tmx1)

end subroutine


end module

