#include "macros.h"

!
! note: most test subroutines are ifdef'd out because they
!       generate compilation errors (pgf90/linux) if the memeory
!       is 256^3 per cpu.  ifc does not complain, so maybe there is
!       a pgf90 compiler option to fix this?
!

subroutine test
use params



!call test_fft_fd
!call test_transform
!call test_fft
!call test_poisson
!call test_divfree
call test_poisson_dirichlet
!call test_poisson_ghost

end subroutine




subroutine test_poisson_dirichlet
use params
use ghost
implicit none
external helmholtz_dirichlet

real*8 b(nx,ny)
real*8 psi(nx,ny)
real*8 psi_exact(nx,ny)
real*8 lpsi(nx,ny)
real*8 work(nx,ny)
real*8 :: zero=0,one=1
real*8 :: tol=1e-12
integer i,j,k,n
b=0
psi=0
work=0
lpsi=0

g_bdy_x1=INFLOW0_ONESIDED
g_bdy_x2=INFLOW0_ONESIDED
g_bdy_y1=INFLOW0_ONESIDED
g_bdy_y2=INFLOW0_ONESIDED
offset_bdy=1
call init_grid()


! psi = x**2 + y**2    
! laplacian = 4
do j=by1,by2
do i=bx1,bx2
   !psi_exact(i,j)=xcord(i)**2 + ycord(j)**2
   !b(i,j)=4

   !psi_exact(i,j)=xcord(i)**3 + ycord(j)**3
   !b(i,j)=3*2*xcord(i) + 3*2*ycord(j)

   psi_exact(i,j)=100 + sin(pi*xcord(i))*sin(pi*ycord(j))
   b(i,j)=-2*pi*pi*sin(pi*xcord(i))*sin(pi*ycord(j))

   !psi_exact(i,j)=cos(pi*xcord(i))*cos(pi*ycord(j))
   !b(i,j)=-2*pi*pi*cos(pi*xcord(i))*cos(pi*ycord(j))
enddo
enddo



psi=psi_exact
! muck around with psi on the interior, to give CG something to iterate on:
do j=inty1,inty2
do i=intx1,intx2
   psi(i,j)=0
enddo
enddo


! copy b.c. from psi into 'b', and apply compact correction to b:
call ghost_update_x(b,1)
call ghost_update_y(b,1)
call helmholtz_dirichlet_setup(b,psi,work,1)


!call jacobi(psi,b,zero,one,tol,work,helmholtz_dirichlet,.false.)
!call cgsolver(psi,b,zero,one,tol,work,helmholtz_dirichlet,.false.)
psi=b; call helmholtz_dirichlet_inv(psi,work,zero,one)


print *,'helmholtz_dirichlet solver error: ',&
maxval(abs(psi(  bx1:bx2,by1:by2)-psi_exact(bx1:bx2,by1:by2)     ))
do j=by1,by2
do i=bx1,bx2
!   print *,i,j,psi(i,j),psi_exact(i,j)
enddo
enddo


#if 0
 jacobi results  COMPACT WORKING!
 compact:      ITER         ERROR
      16x16     513         6.18e-5
      32x32     830         4.12e-6       15x    
      64x64    1795         2.47e-7       17x
 
 2nd:
      16x16   1171           1.31e-2 
      32x32   3447           3.23e-3      4x
      64x64   8059           8.05e-4      4x


#endif
stop
end subroutine



#if 0


subroutine test_poisson_ghost
use params
use ghost
implicit none
external helmholtz_periodic_ghost

real*8 b(nx,ny)
real*8 psi(nx,ny)
real*8 psi_exact(nx,ny)
real*8 lpsi(nx,ny)
real*8 work(nx,ny)
real*8 :: zero=0,one=1
real*8 :: tol=1e-8
integer i,j,k,n
b=0
psi=0
work=0
lpsi=0

! psi = x**2 + y**2    
! laplacian = 4

! set b.c. periodoic in X
! set b.c. reflect-odd in Y

g_bdy_x1=REFLECT
g_bdy_x2=REFLECT
g_bdy_y1=REFLECT_ODD
g_bdy_y2=REFLECT_ODD
call init_grid()

do j=ny1,ny2
do i=nx1,nx2
!   psi_exact(i,j)=xcord(i)/delx
!   b(i,j)=0

   psi_exact(i,j)=cos(pi*xcord(i))*sin(pi*ycord(j))
   b(i,j)=-2*pi*pi*cos(pi*xcord(i))*sin(pi*ycord(j))
enddo
enddo



psi=psi_exact
! muck around with psi to give CG something to iterate on:
psi=0


! copy b.c. from psi into 'b', and apply compact correction to b:
call ghost_update_x(b,1)
call ghost_update_y(b,1)
call helmholtz_dirichlet_setup(b,psi,work,0)

call cgsolver(psi,b,zero,one,tol,work,helmholtz_periodic_ghost,.false.)


print *,'helmholtz_periodic_ghost solver error: ',&
maxval(abs(psi(  nx1:nx2,ny1:ny2)-psi_exact(nx1:nx2,ny1:ny2)     ))
stop

end subroutine










subroutine test_fft_fd
use params
use transpose
implicit none

real*8 work(nx,ny,nz)
real*8 work2(nx,ny,nz)
real*8 input(nx,ny,nz)
real*8 div(nx,ny,nz)
real*8 grad(nx,ny,nz,3)
real*8 dummy(1)
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k,i0,i1,i2,i3,im,j1,k1
character(len=80) message

input = 0

do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
   input(i,j,k)=0
   if (imcord(i)==-2 .and. &
       jmcord(j)== 1 .and. &
       kmcord(k)==-3 ) then
      input(i,j,k)=1
   endif
enddo
enddo
enddo
call ifft3d(input,work)
! compute grad
do i=1,3
   call der(input,grad(1,1,1,i),dummy,work,DX_ONLY,i)
enddo

! compute div(grad)
div=0
do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
   i1=i+1
   if (i1>nx2) i1=i1-nx2
   j1=j+1
   if (j1>ny2) j1=j1-ny2
   k1=k+1
   if (k1>nz2) k1=k1-nz2
   div(i,j,k)=( &
     (grad(i1,j ,k, 1)-grad(i,j, k, 1)) + &
     (grad(i1,j1,k, 1)-grad(i,j1,k, 1)) + &
     (grad(i1,j ,k1,1)-grad(i,j, k1,1)) + &
     (grad(i1,j1,k1,1)-grad(i,j1,k1,1))  ) / (4*delx)

   div(i,j,k)=div(i,j,k) + ( &
     (grad(i ,j1,k ,2)-grad(i ,j,k ,2)) + &
     (grad(i1,j1,k ,2)-grad(i1,j,k ,2)) + &
     (grad(i ,j1,k1,2)-grad(i ,j,k1,2)) + &
     (grad(i1,j1,k1,2)-grad(i1,j,k1,2))  ) / (4*dely)

   div(i,j,k)=div(i,j,k) + ( &
     (grad(i ,j ,k1,3)-grad(i ,j ,k,3)) + &
     (grad(i1,j ,k1,3)-grad(i1,j ,k,3)) + &
     (grad(i ,j1,k1,3)-grad(i ,j1,k,3)) + &
     (grad(i1,j1,k1,3)-grad(i1,j1,k,3))  ) / (4*delz)

enddo
enddo
enddo
!call divergence(div,grad,work,work2)

call fft3d(div,work)
call print_modes(div)

do i=1,3
print *,'grad modes'
call fft3d(grad(1,1,1,i),work)
call print_modes(grad(1,1,1,i))
enddo

stop
end subroutine












subroutine maxvalMPI(error)
use mpi
use params
implicit none
real*8 error,tmp
integer ierr
#ifdef USE_MPI
   tmp=error
   call mpi_allreduce(tmp,error,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
#endif
end subroutine



subroutine test_transform
use params
use transpose
implicit none

real*8 input(nx,ny,nz)
real*8 output(nx,ny,nz)
real*8 work(nx,ny,nz)
real*8 mx,tx1,tx2,tmax,tave
integer n1,n1d,n2,n2d,n3,n3d
integer i,j,k,n,ierr
character(len=80) message

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
call transpose_to_x(input,work,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_x(work,output,n1,n1d,n2,n2d,n3,n3d)
mx=maxval(abs(input-output))
call maxvalMPI(mx)
write(message,*) 'transpose x round trip error:',mx
call print_message(message)


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






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  benchmark the transforms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
n=1

call wallclock(tx1)
do i=1,n
call transpose_to_z(input,work,n1,n1d,n2,n2d,n3,n3d)
call transpose_to_z(work,input,n1,n1d,n2,n2d,n3,n3d)
enddo
call wallclock(tx2)

tmax=tx2-tx1
tave=tx2-tx1
#ifdef USE_MPI
tx2=tmax
call mpi_allreduce(tx2,tmax,1,MPI_REAL8,MPI_MAX,comm_3d,ierr)
tx2=tave
call mpi_allreduce(tx2,tave,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
tave=tave/ncpus
#endif

write(message,'(a,i3,a,2f10.5)') 'wall clock transpose_to_z n=',n,' time=',&
     tave/(2*n),tmax/(2*n)
call print_message(message)


call wallclock(tx1)
do i=1,n
call transpose_from_z(work,output,n1,n1d,n2,n2d,n3,n3d)
call transpose_from_z(output,work,n1,n1d,n2,n2d,n3,n3d)
enddo
call wallclock(tx2)



tmax=tx2-tx1
tave=tx2-tx1
#ifdef USE_MPI
tx2=tmax
call mpi_allreduce(tx2,tmax,1,MPI_REAL8,MPI_MAX,comm_3d,ierr)
tx2=tave
call mpi_allreduce(tx2,tave,1,MPI_REAL8,MPI_SUM,comm_3d,ierr)
tave=tave/ncpus
#endif

write(message,'(a,i3,a,2f10.5)') 'wall clock transpose_from_z n=',n,' time=',&
     tave/(2*n),tmax/(2*n)
call print_message(message)



end subroutine







subroutine test_divfree
use params
use transpose
use fft_interface
implicit none

real*8 work(nx,ny,nz)
real*8 input(nx,ny,nz,3)
real*8 input2(nx,ny,nz,3)
real*8 p(nx,ny,nz)
real*8 div(nx,ny,nz),div2(nx,ny,nz)
real*8 d1(nx,ny,nz)
real*8 d2(nx,ny,nz)
real*8 dummy
character(len=80) message

integer i,j,k,dim,n1,n2,n3,n1d,n2d,n3d

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



call divfree_gridspace(input,p,d1,work)

! compute p = div(u)
i=3
call der(input(1,1,1,i),d1,dummy,work,DX_ONLY,i)
p = d1
i=2
call der(input(1,1,1,i),d1,dummy,work,DX_ONLY,i)
p = p + d1
i=1
call der(input(1,1,1,i),d1,dummy,work,DX_ONLY,i)
p = p + d1

write(message,'(a,e10.5)') 'maxval of divergence = ',&
   maxval(abs(p(nx1:nx2,ny1:ny2,nz1:nz2)))
call print_message(message)
#if 0
print *,'modes of divergence:'
call fft3d(p,work)
call print_modes(p)
call ifft3d(p,work)
#endif
end subroutine







subroutine test_poisson
use params
implicit none

real*8 work(nx,ny,nz)
real*8 d1(nx,ny,nz)
real*8 d2(nx,ny,nz)
real*8 rhs(nx,ny,nz)
real*8 input(nx,ny,nz)
integer i,j,k
real*8 cf1,alpha,beta
real*8 error,tmp
character(len=80) message

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
call helmholtz_periodic_inv(rhs,work,alpha,beta)

work=rhs-input
error=maxval(abs(work(nx1:nx2,ny1:ny2,nz1:nz2)));
call maxvalMPI(error)
call print_message("analytic Helmholtz problem:")
write(message,'(a,f6.2,a,f6.2,a,e15.8)') 'Laplace solver alpha=',alpha,' beta=',beta,&
   '  error=',error
call print_message(message)


! descrete problem


! compute descrete laplacian:
call der(input,d1,d2,work,DX_AND_DXX,1)
rhs=d2
call der(input,d1,d2,work,DX_AND_DXX,2)
rhs=rhs+d2
call der(input,d1,d2,work,DX_AND_DXX,3)
rhs=rhs+d2

rhs=alpha*input + beta*rhs
call helmholtz_periodic_inv(rhs,work,alpha,beta)
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
call print_message("discrete Helmholtz problem:")
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
real*8 cf1,error,xm,ym,zm
integer i,j,k,i0,i1,i2,i3,im
character(len=80) message


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


do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
   cf1=3*2*pi 
   input(i,j,k)=           6*sin(cf1*zcord(k))*sin(cf1*xcord(i))    + input(i,j,k)
   inputz(i,j,k)=      cf1*6*cos(cf1*zcord(k))*sin(cf1*xcord(i))    + inputz(i,j,k)
   inputzz(i,j,k)=-cf1*cf1*6*sin(cf1*zcord(k))*sin(cf1*xcord(i))    + inputzz(i,j,k)
   inputx(i,j,k)=      cf1*6*sin(cf1*zcord(k))*cos(cf1*xcord(i))    + inputx(i,j,k)
   inputxx(i,j,k)=-cf1*cf1*6*sin(cf1*zcord(k))*sin(cf1*xcord(i))    + inputxx(i,j,k)

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
real*8 :: output(nx,ny,nz),wn
character(len=80) message
integer :: i,j,k,count

count=0
do i=nx1,nx2
do j=ny1,ny2
do k=nz1,nz2
    wn=sqrt(real(imcord(i)**2+jmcord(j)**2+kmcord(k)**2))
    if (abs(output(i,j,k)) > 1e-10 .and. count<100 ) then	
       count=count+1
       write(*,'(i4,i4,i4,a,i4,i4,i4,f6.2,a,2e15.7)') i,j,k,&
           ' mode=',imcord(i),jmcord(j),kmcord(k),wn,&
            ' val=',output(i,j,k)
       if (count==100) then
          call print_message("count>100 will not print any more modes")
       endif
    endif
enddo
enddo
enddo
end subroutine




#endif