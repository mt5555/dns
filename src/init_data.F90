subroutine init_data(Q)
use params
use fft_interface
implicit none
real*8 :: Q(nx,ny,nz,n_var)
integer i,j,k

! uniform flow to the right
Q=0
do j=ny1,ny2
   if (ycord(J)<.33 .or. ycord(j)>.66) then
      Q(:,j,:,1)=1
   else
      Q(:,j,:,1)=-1
   endif
enddo

call bc_preloop


end subroutine









