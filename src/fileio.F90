subroutine out(itime,time,Q)
use params
implicit none
real*8 :: Q(nx,ny,nz,n_var)
real*8 :: time
integer :: itime


! local variables
integer i,j,k,n
real*8 xnx,xny,xnz,xnv
character*80 message

write(message,'(f9.4)') 1000.0000 + time
message = "data" // message(2:9) // ".out"
open(unit=10,file=message,form='binary')

write(10) time
xnv=n_var
xnx=nslabx
xny=nslaby
xnz=nslabz
write(10) xnx,xny,xnz,xnv
do n=1,n_var
do k=nz1,nz2
do j=ny1,ny2
do i=nx1,nx2
   write(10) Q(i,j,k,n)	
enddo
enddo
enddo
enddo
close(10)

end subroutine

