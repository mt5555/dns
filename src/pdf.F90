#include "macros.h"

module structf
implicit none

integer,parameter :: delta_num=10
integer,parameter :: pdf_max_bin=1000
real*8            :: pdf_bin_size = .1

integer           :: delta_val(delta_num)





type pdf_structure_function
!
! pdf%data(-bin:bin,delta_num,npdfs)
!
integer,allocatable :: pdf(:,:)
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

str%nbin=bin
allocate(str%pdf(-bin:bin,delta_num))
str%pdf=0

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

! make a copy of data
allocate(str2%pdf(-n:n,delta_num))
str2%pdf=str%pdf

! create a larger structure function
deallocate(str%pdf)
allocate(str%pdf(-bin:bin,delta_num))

! copy data into new, larger structure function
str%pdf(-n:n,:)=str2%pdf

! delete the copy of original data
deallocate(str2%pdf)

end subroutine





subroutine compute_pdf(Q,n1,n1d,n2,n2d,n3,n3d,str)
!
! compute a pdf_structure function along the first dimension of Q
! for all the values of delta given by delta_val(:)
!
implicit none
real*8 :: n1,n1d,n2,n2d,n3,n3d
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
end subroutine


end module
