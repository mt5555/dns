module mod_fft_interface

use mod_params, only: g_maxn
implicit none
integer, parameter ::  num_fftsizes=3
integer fftsizes(num_fftsizes)
real*8  trigs(3*g_maxn/2+1,num_fftsizes)
integer ifax(13,num_fftsizes)


private :: fftinit, getindex


contains 



subroutine fftinit(n,index)
integer n,index

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

fftsizes(index)=n
call set99(trigs(1,index),ifax(1,index),fftsizes(index))
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

   if (fftsizes(i)==0) then
      call fftinit(n1,i)      
      exit
   endif
   if (n1==fftsizes(i)) exit
enddo
end subroutine






subroutine ifft(p,n1,n1d,n2,n2d,n3,n3d)
real*8 p(n1d,n2d,n3d)
real*8 w(n2*(n1+1))
integer n1,n1d,n2,n2d,n3,n3d
character*80 message_str

integer i,k
call getindex(n1,i)
do k=1,n3
   call fft991(p(1,1,k),w,trigs(1,i),ifax(1,i),1,n1d,n1,n2,1)
enddo
end subroutine



subroutine fft(p,n1,n1d,n2,n2d,n3,n3d)
real*8 p(n1d,n2d,n3d)
real*8 w(n2*(n1+1))
integer n1,n1d,n2,n2d,n3,n3d

integer i,k

call getindex(n1,i)
do k=1,n3
   call fft991(p(1,1,k),w,trigs(1,i),ifax(1,i),1,n1d,n1,n2,-1)
enddo
end subroutine




subroutine fft_derivatives(px,pxx,numder,n1,n1d,n2,n2d,n3,n3d)
!
!  input is original given in px.  
!  if numder=1   compute d/dx, return in px
!  if numder=2   compute d2/dx2, return in pxx
!
real*8 px(n1d,n2d,n3d)
real*8 pxx(n1d,n2d,n3d)
integer numder,n1,n1d,n2,n2d,n3,n3d

integer i,j,k,m


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




end ! module mod_fft_interface
