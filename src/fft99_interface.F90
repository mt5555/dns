!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! our wrapper for ECMWF FFT99.
!
! provides public interfaces:
!   fft_interface_init               call this before using any other routines
!   fft
!   ifft
!   fft_derivatives
!   fft_laplace_inverse
! 
! Routines work on data of the form:  p(n1d,n2d,n3d)
! Size of the grid point data         p(1:n1,1:n2,1:n3)
! Size of fourier coefficients        p(1:n1+2,1:n2+2,1:n3+2)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
module fft_interface

implicit none
integer, parameter ::  num_fftsizes=3

integer :: init=0
type fftdata_d
   real*8,dimension(:),pointer :: trigs
   integer :: ifax(13)
   integer :: size
end type
type(fftdata_d) :: fftdata(num_fftsizes)


integer, parameter ::  fftblocks=100   ! do fft's in blocks of size fftblocks
                                       ! set very large to disable
private :: fftinit, getindex
contains 


#if 0
output:  
     if isign = +1, and m coefficient vectors are supplied
     each containing the sequence:

     a(0),b(0),a(1),b(1),...,a(n/2),b(n/2)  (n+2 values)

     then the result consists of m data vectors each
     containing the corresponding n+2 gridpoint values:

     x(0), x(1), x(2),...,x(n-1),0,0.

     note: the fact that the gridpoint values x(j) are real
     implies that b(0)=b(n/2)=0.  for a call with isign=+1,
     it is not actually necessary to supply these zeros.
#endif


subroutine fft_interface_init()
integer i
do i=1,num_fftsizes
   fftdata(i).size=0	
enddo
init=1
end subroutine


subroutine fftinit(n,index)
integer n,index
real*8,allocatable,target  :: trigdata(:)

if (init==0) call abort("fft99_interface.F90: call fft_interface_init to initialize first!");
if (n>1000000) call abort("fft99_interface.F90: n>1 million")

allocate(trigdata(3*n/2+1))
fftdata(index).size=n
call set99(trigdata,fftdata(index).ifax,n)
fftdata(index).trigs => trigdata


end subroutine




subroutine getindex(n1,index)
integer n1,index

character*80 message_str
integer i

i=0
do 
   i=i+1
   if (i>num_fftsizes) then
      write(message_str,'(a,i10)') "fft_interface.F90:  Failed initializing an fft of size =",n1
      call abort(message_str)
   endif

   if (fftdata(i).size==0) then
      call fftinit(n1,fftdata(i).size)      
      exit
   endif
   if (n1==fftdata(i).size) exit
enddo
end subroutine






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute 3D in-place iFFT of p
! FFT taken along first direction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ifft(p,n1,n1d,n2,n2d,n3,n3d)
real*8 p(n1d,n2d,n3d)
real*8 w(min(fftblocks,n2)*(n1+1))
integer n1,n1d,n2,n2d,n3,n3d
character*80 message_str

integer index,j,k,numffts
call getindex(n1,index)

j=0  ! j=number of fft's computed for each k
do k=1,n3
   j=0  ! j=number of fft's computed for each k
   do while (j<n2)
      numffts=min(fftblocks,n2-j)	
      call fft991(p(1,1,k),w,fftdata(index).trigs,fftdata(index).ifax,1,n1d,n1,numffts,1)
      j=j+numffts
   enddo
enddo

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Compute 3D in-place FFT of p
! FFT taken along first direction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft(p,n1,n1d,n2,n2d,n3,n3d)
real*8 p(n1d,n2d,n3d)
real*8 w(min(fftblocks,n2)*(n1+1))
integer n1,n1d,n2,n2d,n3,n3d

integer index,j,k,numffts

call getindex(n1,index)
do k=1,n3
   j=0  ! j=number of fft's computed for each k
   do while (j<n2)
      numffts=min(fftblocks,n2-j)	
      call fft991(p(1,1,k),w,fftdata(index).trigs,fftdata(index).ifax,1,n1d,n1,n2,-1)
      j=j+numffts
   enddo
enddo
end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Compute FFT derivates.
!
!  input: px 
!  output:
!     if numder=1   return d/dx along first direction in px
!                   (and pxx is not accessed) 
!     if numder=2   return d2/dx2 along first direction in pxx
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_derivatives(px,pxx,numder,n1,n1d,n2,n2d,n3,n3d)
real*8 px(n1d,n2d,n3d)
real*8 pxx(n1d,n2d,n3d)
integer numder,n1,n1d,n2,n2d,n3,n3d

integer i,j,k,m
real*8 temp


call fft(px,n1,n1d,n2,n2d,n3,n3d)

if (numder>=2) then
   do k=1,n3
   do j=1,n2
   do m = 0, n1/2
      i = 2*m+1
      pxx(i,j,k) = -m*m * px(i,j,k)
      pxx(i+1,j,k) = -m*m * px(i+1,j,k)
   enddo
   enddo
   enddo
   call ifft(pxx,n1,n1d,n2,n2d,n3,n3d)
endif

if (numder>=1) then
   do k=1,n3
   do j=1,n2
   do m = 0, n1/2
      i = 2*m+1
      temp =  m * px(i,j,k)
      px(i,j,k) = -m * px(i+1,j,k)
      px(i+1,j,k) = temp
   enddo
   enddo
   enddo
   call ifft(px,n1,n1d,n2,n2d,n3,n3d)
endif
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Solve  [alpha + beta*Laplacian] p = rhs
!
! on input,  p = fourier coefficients of rhs
! on output, p = fourier coefficients of solution
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_laplace_inverse(p,n1,n1d,n2,n2d,n3,n3d,alpha,beta)
real*8 p(n1d,n2d,n3d)
real*8 alpha,beta
integer n1,n1d,n2,n2d,n3,n3d

integer i,j,k,im,jm,km
real*8 xfac

   do k=1,n3+2
      km=(k-1)/2
      do j=1,n2+2
         jm=(j-1)/2
         do i=1,n2+2
            im=(i-1)/2
            xfac= alpha + beta*(-im*im -km*km - jm*jm)      
            if (xfac<>0) xfac = 1/xfac
            p(i,j,k)=p(i,j,k)*xfac
         enddo
      enddo
   enddo

end subroutine




end ! module mod_fft_interface
