#include "macros.h"
subroutine test
use params

call test_transform
call test_fft
call test_poisson
call test_divfree

end subroutine



subroutine maxvalMPI(error)
use mpi
use params
implicit none
real*8 error,tmp
integer ierr
#ifdef USE_MPI
   tmp=error
   call MPI_allreduce(tmp,error,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
end subroutine





subroutine test_transform
use params
use transpose
implicit none

real*8 input(nx,ny,nz)
real*8 output(nx,ny,nz)
real*8 work(g_nz2,nx,ny)
real*8 mx
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k
character*80 message

input=0
output=0
do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
   input(i,j,k)=zcord(k)+xcord(i)+ycord(j)
enddo
enddo
enddo

output=0
call transpose_to_z(input,work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(work,output,n1,n1d,n2,n2d,n3,n3d)
mx=maxval(abs(input-output))
call maxvalMPI(mx)
write(message,*) 'transpose z round trip error:',mx
call print_message(message)

output=0
call transpose_to_y(input,work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_y(work,output,n1,n1d,n2,n2d,n3,n3d)
mx=maxval(abs(input-output))
call maxvalMPI(mx)
write(message,*) 'transpose y round trip error:',mx
call print_message(message)

output=0
call transpose_to_x(input,work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_x(work,output,n1,n1d,n2,n2d,n3,n3d)
mx=maxval(abs(input-output))
call maxvalMPI(mx)
write(message,*) 'transpose x round trip error:',mx
call print_message(message)


#if 0
if (my_z==1) then
   print *,'rank=',my_pe,'my coords  = ',my_x,my_y,my_z
   print *,'maxval round trip=',maxval(abs(input-output))
endif

if (my_z==1) then
   print *,'rank=',my_pe,'my coords  = ',my_x,my_y,my_z
   do k=1,g_nz2
      write(*,'(a,i5,3f10.4)') 'k,work ',k,work(k,1),work(k,1)-work(k-1,1)
   enddo
endif
#endif
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
do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2

input(i,j,k,1)= xcord(i)*(xcord(i)-1)*ycord(j)*(ycord(j)-1)  
!input(i,j,k,1)= sin(2*3*pi*xcord(i))*cos(2*4*pi*ycord(j))

input(i,j,k,2)= exp((xcord(i)-1))*exp(ycord(j))
!input(i,j,k,2)= sin(2*3*pi*xcord(i))*cos(2*4*pi*ycord(j))

if (nslabz>1)  then
   input(i,j,k,3)= xcord(i)*(xcord(i)-1)*ycord(j)*(ycord(j)-1)  
endif

enddo
enddo
enddo

#if 0
! remove that pesky highest cosine mode
! no longer needed, we tweaked laplacian so that it = div grad
do i=1,3
   call fft3d(input(1,1,1,i),work)
   call fft_filter_last(input(1,1,1,i))
   call ifft3d(input(1,1,1,i),work)
enddo
#endif



call divfree(input,p)

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
real*8 error,tmp
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

   if (g_nz>=6) then
      cf1=2*2*pi 
      input(i,j,k)=           5*cos(cf1*zcord(k))    + input(i,j,k)
      rhs(i,j,k)=    -cf1*cf1*5*cos(cf1*zcord(k))    + rhs(i,j,k)
   endif
   if (g_nz>=8) then
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


work=rhs-input
#if 0
do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
   if (abs(work(i,j,k))>1e-9) then
	write(*,'(3i5,2e15.7)') i,j,k,rhs(i,j,k),input(i,j,k)
   endif	
enddo
enddo
enddo
#endif
error=maxval(abs(work(nx1:nx2,ny1:ny2,nz1:nz2)));
call maxvalMPI(error)
write(message,'(a,f6.2,a,f6.2,a,e15.8)') 'Laplace solver alpha=',alpha,' beta=',beta,&
   '  error=',error
call print_message(message)


#if 0
call print_message("Laplace solver data:")
call print_message("analytic solution")
call fft3d(input,work)
call print_modes(input)
call ifft3d(input,work)

call print_message("computed solution")
call fft3d(rhs,work)
call print_modes(rhs)
call ifft3d(rhs,work)
#endif




end subroutine




subroutine test_fft
use params
use transpose
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
integer n1,n1d,n2,n2d,n3,n3d
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
enddo
enddo
enddo


do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
   cf1=2*pi
   input(i,j,k)=           2*sin(cf1*xcord(i))  + input(i,j,k)
   inputx(i,j,k)=      cf1*2*cos(cf1*xcord(i))  + inputx(i,j,k)
   inputxx(i,j,k)=-cf1*cf1*2*sin(cf1*xcord(i))  + inputxx(i,j,k)
enddo
enddo
enddo


do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
   cf1=3*2*pi
   input(i,j,k)=           3*cos(cf1*xcord(i))    + input(i,j,k)
   inputx(i,j,k)=     -cf1*3*sin(cf1*xcord(i))    + inputx(i,j,k)
   inputxx(i,j,k)=-cf1*cf1*3*cos(cf1*xcord(i))    + inputxx(i,j,k)
enddo
enddo
enddo



do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
   cf1=5*2*pi 
   input(i,j,k)=           4*cos(cf1*ycord(j))    + input(i,j,k)
   inputy(i,j,k)=     -cf1*4*sin(cf1*ycord(j))    + inputy(i,j,k)
   inputyy(i,j,k)=-cf1*cf1*4*cos(cf1*ycord(j))    + inputyy(i,j,k)

enddo
enddo
enddo


do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
   cf1=3*2*pi 
   input(i,j,k)=           5*cos(cf1*zcord(k))    + input(i,j,k)
   inputz(i,j,k)=     -cf1*5*sin(cf1*zcord(k))    + inputz(i,j,k)
   inputzz(i,j,k)=-cf1*cf1*5*cos(cf1*zcord(k))    + inputzz(i,j,k)
enddo
enddo
enddo





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! d/dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call der(input,px,pxx,work,DX_AND_DXX,1)

error=maxval(abs(px(nx1:nx2,ny1:ny2,nz1:nz2)-inputx(nx1:nx2,ny1:ny2,nz1:nz2)))
call maxvalMPI(error)
write(message,'(a,e15.10)') "x-direction d/dx error=",error
call print_message(message)


error=maxval(abs(pxx(nx1:nx2,ny1:ny2,nz1:nz2)-inputxx(nx1:nx2,ny1:ny2,nz1:nz2)))
call maxvalMPI(error)
write(message,'(a,e15.10)') "x-direction d2/dxx error=",error
call print_message(message)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! d/dy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call der(input,px,pxx,work,DX_AND_DXX,2)

error=maxval(abs(px(nx1:nx2,ny1:ny2,nz1:nz2)-inputy(nx1:nx2,ny1:ny2,nz1:nz2)))
call maxvalMPI(error)
write(message,'(a,e15.10)') "y-direction d/dy error=",error
call print_message(message)

error=maxval(abs(pxx(nx1:nx2,ny1:ny2,nz1:nz2)-inputyy(nx1:nx2,ny1:ny2,nz1:nz2)))
call maxvalMPI(error)
write(message,'(a,e15.10)') "y-direction d2/dyy error=",error
call print_message(message)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! d/dz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call der(input,px,pxx,work,DX_AND_DXX,3)

error=maxval(abs(px(nx1:nx2,ny1:ny2,nz1:nz2)-inputz(nx1:nx2,ny1:ny2,nz1:nz2)))
call maxvalMPI(error)
write(message,'(a,e15.10)') "z-direction d/dz error=",error
call print_message(message)

error=maxval(abs(pxx(nx1:nx2,ny1:ny2,nz1:nz2)-inputzz(nx1:nx2,ny1:ny2,nz1:nz2)))
call maxvalMPI(error)
write(message,'(a,e15.10)') "z-direction d2/dzz error=",error
call print_message(message)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3d fft
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
output=input
call fft3d(output,work)
call print_modes(output)
call ifft3d(output,work)

error=maxval(abs(input(nx1:nx2,ny1:ny2,nz1:nz2)-output(nx1:nx2,ny1:ny2,nz1:nz2)))
call maxvalMPI(error)
write(message,'(a,e15.10)') "x-direction Forward-Backward 3D FFT: error=",error
call print_message(message)



end subroutine




subroutine print_modes(output)
use params
implicit none
real*8 :: output(nx,ny,nz)
character*80 message
integer :: i,j,k,count=0

do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
    if (abs(output(i,j,k)) > 1e-9 .and. count<50) then	
       count=count+1
       write(message,'(i4,i4,i4,a,i4,i4,i4,a,2e15.7)') i,j,k,'mode=',imcord(i),jmcord(j),kmcord(k),&
            ' val=',output(i,j,k)
       call print_message(message)
       if (count==50) then
          call print_message("count>50 will not print any more modes")
       endif
    endif
enddo
enddo
enddo
end subroutine




