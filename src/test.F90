#include "macros.h"
subroutine test

call test_fft
call test_poisson
call test_divfree

end subroutine






subroutine test_divfree
use params
use fft_interface
implicit none

real*8 work(nx,ny,nz)
real*8 input(nx,ny,nz,3)
real*8 p(nx,ny,nz)
real*8 d1(nx,ny,nz)
real*8 dummy
character*80 message

integer i,j,k,dim

input=0
do dim=1,3
do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2

if (dim==1) input(i,j,k,dim)= xcord(i)*(xcord(i)-1)*ycord(j)*(ycord(j)-1)  

enddo
enddo
enddo
enddo

! remove that pesky highest cosine mode
call fft3d(input,work)
call fft_filter(input)
call ifft3d(input,work)


call divfree(input)

! compute p = div(u)
i=1
call der(input(1,1,1,i),d1,dummy,work,DX_ONLY,i)
p = d1
i=2
call der(input(1,1,1,i),d1,dummy,work,DX_ONLY,i)
p = p + d1
i=3
call der(input(1,1,1,i),d1,dummy,work,DX_ONLY,i)
p = p + d1

write(message,'(a,e10.5)') 'maxval of divergence = ',maxval(p(nx1:nx2,ny1:ny2,nz1:nz2))
call print_message(message)
end subroutine







subroutine test_poisson
use params
implicit none

real*8 work(nx,ny,nz)
real*8 rhs(nx,ny,nz)
real*8 input(nx,ny,nz)
integer i,j,k
real*8 cf1,alpha,beta
real*8 error
character*80 message

alpha=1
beta=.5

rhs = 0
input=.33

if (alpha==0) input=0  ! in this case there is an arbritrary constant

do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2

   cf1=2*pi
   input(i,j,k)=         2*sin(cf1*xcord(i))  + input(i,j,k)
   rhs(i,j,k)=  -cf1*cf1*2*sin(cf1*xcord(i))  + rhs(i,j,k)

   cf1=3*2*pi
   input(i,j,k)=           3*cos(cf1*xcord(i)) + input(i,j,k)
   rhs(i,j,k)=    -cf1*cf1*3*cos(cf1*xcord(i))  + rhs(i,j,k)


   cf1=5*2*pi 
   input(i,j,k)=           4*cos(cf1*ycord(j))    + input(i,j,k)
   rhs(i,j,k)=    -cf1*cf1*4*cos(cf1*ycord(j))    + rhs(i,j,k)

   if (nz2>1) then
      cf1=4*2*pi 
      input(i,j,k)=           5*cos(cf1*zcord(k))    + input(i,j,k)
      rhs(i,j,k)=    -cf1*cf1*5*cos(cf1*zcord(k))    + rhs(i,j,k)

      cf1=3*2*pi 
      input(i,j,k)=           2*sin(cf1*zcord(k))    + input(i,j,k)
      rhs(i,j,k)=    -cf1*cf1*2*sin(cf1*zcord(k))    + rhs(i,j,k)

   endif

   cf1=1*2*pi 
   input(i,j,k)= sin(cf1*xcord(i))*cos(cf1*ycord(j))    + input(i,j,k)
   rhs(i,j,k)=   -cf1*cf1*sin(cf1*xcord(i))*cos(cf1*ycord(j))   &
                 +sin(cf1*xcord(i))*(-cf1*cf1*cos(cf1*ycord(j))) &
                   + rhs(i,j,k)



enddo
enddo
enddo



rhs = alpha*input + beta*rhs
call poisson(rhs,work,alpha,beta)

#if 0
call fft3d(rhs,work)
call print_modes(rhs)
call ifft3d(rhs,work)
#endif



work=rhs-input
#if 0
do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
   if (abs(work(i,j,k))>1e-9) then
	print *,i,j,k,rhs(i,j,k),input(i,j,k)
   endif	
enddo
enddo
enddo
#endif
error=maxval(abs(work(nx1:nx2,ny1:ny2,nz1:nz2)));
write(message,'(a,f6.2,a,f6.2,a,e15.8)') 'Laplace solver alpha=',alpha,' beta=',beta,&
   '  error=',error
call print_message(message)





end subroutine




subroutine test_fft
use params
implicit none

real*8 work(nx,ny,nz)
real*8 output(nx,ny,nz)
real*8 px(nx,ny,nz)
real*8 pxx(nx,ny,nz)
real*8 input(nx,ny,nz)
real*8 inputx(nx,ny,nz)
real*8 inputxx(nx,ny,nz)
real*8 inputy(nx,ny,nz)
real*8 inputyy(nx,ny,nz)
real*8 inputz(nx,ny,nz)
real*8 inputzz(nx,ny,nz)
real*8 cf1,error
integer i,j,k
character*80 message


input = 0
inputx = 0
inputxx=0
inputy = 0
inputyy=0
inputz = 0
inputzz=0


do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
   input(i,j,k)=.333

   cf1=2*pi
   input(i,j,k)=           2*sin(cf1*xcord(i))  + input(i,j,k)
   inputx(i,j,k)=      cf1*2*cos(cf1*xcord(i))  + inputx(i,j,k)
   inputxx(i,j,k)=-cf1*cf1*2*sin(cf1*xcord(i))  + inputxx(i,j,k)

   cf1=3*2*pi
   input(i,j,k)=           3*cos(cf1*xcord(i))    + input(i,j,k)
   inputx(i,j,k)=     -cf1*3*sin(cf1*xcord(i))    + inputx(i,j,k)
   inputxx(i,j,k)=-cf1*cf1*3*cos(cf1*xcord(i))    + inputxx(i,j,k)


   cf1=5*2*pi 
   input(i,j,k)=           4*cos(cf1*ycord(j))    + input(i,j,k)
   inputy(i,j,k)=     -cf1*4*sin(cf1*ycord(j))    + inputy(i,j,k)
   inputyy(i,j,k)=-cf1*cf1*4*cos(cf1*ycord(j))    + inputyy(i,j,k)


   cf1=4*2*pi 
   input(i,j,k)=           5*cos(cf1*zcord(k))    + input(i,j,k)
   inputz(i,j,k)=     -cf1*5*sin(cf1*zcord(k))    + inputz(i,j,k)
   inputzz(i,j,k)=-cf1*cf1*5*cos(cf1*zcord(k))    + inputzz(i,j,k)

enddo
enddo
enddo





output=input
call fft3d(output,work)
call print_modes(output)
call ifft3d(output,work)

error=maxval(abs(input(nx1:nx2,ny1:ny2,nz1:nz2)-output(nx1:nx2,ny1:ny2,nz1:nz2)))
write(message,'(a,e15.10)') "x-direction Forward-Backward 3D FFT: error=",error
call print_message(message)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! d/dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call der(input,px,pxx,work,DX_AND_DXX,1)

error=maxval(abs(px(nx1:nx2,ny1:ny2,nz1:nz2)-inputx(nx1:nx2,ny1:ny2,nz1:nz2)))
write(message,'(a,e15.10)') "x-direction d/dx error=",error
call print_message(message)

error=maxval(abs(pxx(nx1:nx2,ny1:ny2,nz1:nz2)-inputxx(nx1:nx2,ny1:ny2,nz1:nz2)))
write(message,'(a,e15.10)') "x-direction d2/dxx error=",error
call print_message(message)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! d/dy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call der(input,px,pxx,work,DX_AND_DXX,2)

error=maxval(abs(px(nx1:nx2,ny1:ny2,nz1:nz2)-inputy(nx1:nx2,ny1:ny2,nz1:nz2)))
write(message,'(a,e15.10)') "y-direction d/dy error=",error
call print_message(message)

error=maxval(abs(pxx(nx1:nx2,ny1:ny2,nz1:nz2)-inputyy(nx1:nx2,ny1:ny2,nz1:nz2)))
write(message,'(a,e15.10)') "y-direction d2/dyy error=",error
call print_message(message)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! d/dz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call der(input,px,pxx,work,DX_AND_DXX,3)

error=maxval(abs(px(nx1:nx2,ny1:ny2,nz1:nz2)-inputz(nx1:nx2,ny1:ny2,nz1:nz2)))
write(message,'(a,e15.10)') "z-direction d/dz error=",error
call print_message(message)

error=maxval(abs(pxx(nx1:nx2,ny1:ny2,nz1:nz2)-inputzz(nx1:nx2,ny1:ny2,nz1:nz2)))
write(message,'(a,e15.10)') "z-direction d2/dzz error=",error
call print_message(message)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






end subroutine




subroutine print_modes(output)
use params
implicit none
real*8 :: output(nx,ny,nz)
character*80 message
integer i,j,k

do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
    if (abs(output(i,j,k)) > 1e-9) then	
       if (mod(i,2)==0) then
          write(message,'(a,i4,i4,i4,a,2f15.10)') '  sine mode=',(i-1)/2,(j-1)/2,(k-1)/2,&
               ' val=',output(i,j,k)
       else
          write(message,'(a,i4,i4,i4,a,2f15.10)') 'cosine mode=',(i-1)/2,(j-1)/2,(k-1)/2,&
               ' val=',output(i,j,k)
        endif
       call print_message(message)
    endif
enddo
enddo
enddo
end subroutine




