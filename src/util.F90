subroutine abort(message)
implicit none
character*(*) message
write(*,*) message
halt
end subroutine


subroutine transpose12(p,pt,nx,nxd,ny,nyd,nz,nzd)
implicit none
real*8 p(nxd,nyd,nzd)
real*8 pt(nyd,nxd,nzd)
integer nx,ny,nz,nxd,nyd,nzd

do k=1,nz
do j=1,ny
do i=1,nx
   pt(j,i,k)=p(i,j,k)
enddo
enddo
enddo

i=ny
ny=nx
nx=i

i=nyd
nyd=nxd
nxd=i
end



subroutine transpose13(p,pt,nx,nxd,ny,nyd,nz,nzd)
implicit none
real*8 p(nxd,nyd,nzd)
real*8 pt(nzd,nyd,nxd)
integer nx,ny,nz,nxd,nyd,nzd

integer i,j,k

do i=1,nx
do j=1,ny
do k=1,nz
   pt(k,j,i)=p(i,j,k)
enddo
enddo
enddo

i=nz
nz=nx
nx=i

i=nzd
nzd=nxd
nxd=i

end







