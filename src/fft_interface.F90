module mod_fft_interface

real*8 trigs(3*maxn/2+1,3)
integer ifax(13,3)




subroutine fftinit()
implicit none
use params

integer ndims(3)  ! number of fft initilizations


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
endif

! check that we have enough padding for FFT
ASSERT(nx+2<=nxd)
ASSERT(ny+2<=nyd)
ASSERT(nz+2<=nzd)


i=i+1
ndims(i)=nx
if (ny<>nx) then
   i=i+1
   ndims(i)=ny
endif
if (nz<>nx & nz<>ny) then
   i=i+1
   ndims(i)=nz
endif

do j=1,i
   call set99(trigs(1,j),ifax(1,j),ndims(j))
enddo
end


subroutine ifft(p,n1,n1d,n2,n2d,n3,n3d)
implicit none
real*8 p(n1d,n2d,n3d)
real*8 w(n2*(n1+1))

find i so that n1==ndims(i) 

m = n2  ! number of tranforms in each z slab
do k=1,n3
   fft991(p(1,1,k),w,trigs(1,i),ifax(1,i),1,n1d,n1,m,1)
end
end


subroutine fft(p,n1,n1d,n2,n2d,n3,n3d)
implicit none
real*8 p(n1d,n2d,n3d)
real*8 w(n2*(n1+1))

find i so that n1==ndims(i) 

m = n2  ! number of tranforms in each z slab
do k=1,n3
   fft991(p(1,1,k),w,trigs(1,i),ifax(1,i),1,n1d,n1,m,-1)
end
end




fft_derivatives(px,pxx,numder,n1,n1d,n2,n2d,n3,n3d)
!
!  input is original given in px.  
!  if numder=1   compute d/dx, return in px
!  if numder=2   compute d2/dx2, return in pxx
!
implicit none

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
end



subroutine fft_laplace_inverse3d()
implicit none
integer i,j,k
real*8 xfac

   do k=1,n3+2
   do j=1,n2+2
   do i=1,n2+2
      im=(i-1)/2
      jm=(j-1)/2
      km=(k-1)/2
      xfac= -im*im -km*km - jm*jm      
      if (xfac<0) xfac = 1/xfac
      p(i,j,k)=p(i,j,k)*xfac
   enddo
   enddo
   enddo

end




end module mod_fft_interface