subroutine init_data(Q)
use params
use fft_interface
implicit none
real*8 :: Q(nx,ny,nz,n_var)

! local variables
integer i,j,k
real*8 eps

! uniform flow to the right
Q=0
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   eps=.05*sin(2*pi*xcord(i))
   if (ycord(J)<.33+eps .or. ycord(j)>.66+eps) then
      Q(i,j,k,1)=1
   else
      Q(i,j,k,1)=-1
   endif
enddo
enddo
enddo

call bc_preloop


end subroutine









