subroutine test

call test_fft
!call test_poisson
call test_divfree

end subroutine






subroutine test_divfree
use params
implicit none

real*8 work(nx,ny,nz)
real*8 input(nx,ny,nz,3)
real*8 p(nx,ny,nz)
real*8 d1(nx,ny,nz)
real*8 dummy

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

call filter(input,work)  ! remove that pesky highest cosine mode
call divfree(input)

! compute p = div(u)
i=1
call der(input(1,1,1,i),d1,dummy,work,1,i)
p = d1
i=2
call der(input(1,1,1,i),d1,dummy,work,1,i)
p = p + d1
i=3
call der(input(1,1,1,i),d1,dummy,work,1,i)
p = p + d1

print *,'maxval of divergence = ',maxval(p(nx1:nx2,ny1:ny2,nz1:nz2))
end subroutine







subroutine test_poisson
use params
implicit none

real*8 work(nx,ny,nz)
real*8 input(nx,ny,nz)
real*8 rhs(nx,ny,nz)
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


   cf1=4*2*pi 
   input(i,j,k)=           5*cos(cf1*zcord(k))    + input(i,j,k)
   rhs(i,j,k)=    -cf1*cf1*5*cos(cf1*zcord(k))    + rhs(i,j,k)

   cf1=3*2*pi 
   input(i,j,k)=           2*sin(cf1*zcord(k))    + input(i,j,k)
   rhs(i,j,k)=    -cf1*cf1*2*sin(cf1*zcord(k))    + rhs(i,j,k)

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
do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
   if (abs(work(i,j,k))>1e-9) then
	print *,i,j,k,rhs(i,j,k),input(i,j,k)
   endif	
enddo
enddo
enddo
error=maxval(work(nx1:nx2,ny1:ny2,nz1:nz2));
write(message,'(a,f6.2,a,f6.2,a,e15.10)') 'Laplace solver alpha=',alpha,' beta=',beta,&
   '  error=',error
call print_message(message)





end subroutine




subroutine test_fft
use params
implicit none

real*8 work(nx,ny,nz)
real*8 input(nx,ny,nz)
real*8 inputx(nx,ny,nz)
real*8 inputxx(nx,ny,nz)
real*8 inputy(nx,ny,nz)
real*8 inputyy(nx,ny,nz)
real*8 inputz(nx,ny,nz)
real*8 inputzz(nx,ny,nz)
real*8 cf1
integer i,j,k
integer n1,n1d,n2,n2d,n3,n3d


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


call print_message("direction1: ")
call fftest_direction1(input,inputx,inputxx,nx2,nx,ny2,ny,nz2,nz)
call print_message("")



n1=nx2
n1d=nx
n2=ny2
n2d=ny
n3=nz2
n3d=nz


call transpose12(input,work,0,n1,n1d,n2,n2d,n3,n3d)  ! x,y,z -> y,x,z
input=work
call transpose12(inputy,work,0,n1,n1d,n2,n2d,n3,n3d)  ! x,y,z -> y,x,z
inputy=work
call transpose12(inputyy,work,1,n1,n1d,n2,n2d,n3,n3d)  ! x,y,z -> y,x,z
inputyy=work
print *,'direction2: ',n1,n1d,n2,n2d,n3,n3d
call fftest_direction1(input,inputy,inputyy,n1,n1d,n2,n2d,n3,n3d)
call transpose12(input,work,1,n1,n1d,n2,n2d,n3,n3d)  ! x,y,z -> y,x,z
input=work
call print_message("");



call transpose13(input,work,0,n1,n1d,n2,n2d,n3,n3d)  
input=work
call transpose13(inputz,work,0,n1,n1d,n2,n2d,n3,n3d)  
inputz=work
call transpose13(inputzz,work,1,n1,n1d,n2,n2d,n3,n3d)  
inputzz=work
print *,'direction3',n1,n1d,n2,n2d,n3,n3d
call fftest_direction1(input,inputz,inputzz,n1,n1d,n2,n2d,n3,n3d)
call transpose13(input,work,1,n1,n1d,n2,n2d,n3,n3d)  
input=work
print *,'after direction3',n1,n1d,n2,n2d,n3,n3d
call print_message("");









end subroutine










subroutine fftest_direction1(input,inputx,inputxx,n1,n1d,n2,n2d,n3,n3d)
use fft_interface
implicit none
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k
real*8 input(n1d,n2d,n3d)
real*8 inputx(n1d,n2d,n3d)
real*8 inputxx(n1d,n2d,n3d)
real*8 output(n1d,n2d,n3d)
real*8 px(n1d,n2d,n3d)
real*8 pxx(n1d,n2d,n3d)
real*8 error
character*80 message

output=input
px=input

call fft(output,n1,n1d,n2,n2d,n3,n3d)
#if 1
do i=1,n1+2
do j=1,1  !n2
do k=1,1  !n3
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
#endif


call ifft(output,n1,n1d,n2,n2d,n3,n3d)
error=maxval(abs(input(1:n1,1:n2,1:n3)-output(1:n1,1:n2,1:n3)))
write(message,'(a,e15.10)') "x-direction Forward-Backward FFT: error=",error
call print_message(message)



call fft_derivatives(px,pxx,2,n1,n1d,n2,n2d,n3,n3d)
error=maxval(abs(px(1:n1,1:n2,1:n3)-inputx(1:n1,1:n2,1:n3)))
write(message,'(a,e15.10)') "x-direction d/dx error=",error
call print_message(message)


error=maxval(abs(pxx(1:n1,1:n2,1:n3)-inputxx(1:n1,1:n2,1:n3)))
write(message,'(a,e15.10)') "x-direction d2/dxx error=",error
call print_message(message)


end subroutine