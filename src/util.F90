
fft_derivatives(Q,Qx,Qxx,numder,n1,n1d,n2,n2d,n3,n3d)
!
!  if numder=1   compute Qx
!  if numder=2   compute Qx and Qxx
!
implicit none

call fft(Q,Qhat,n1,n1d,n2,n2d,n3,n3d)

if (numder>=1) then
   Qhat *= m
   call ifft(Qhat,Qx,n1,n1d,n2,n2d,n3,n3d)
endif
if (numder>=2) then
   Qxx *= m
   call ifft(Qhat,Qxx,n1,n1d,n2,n2d,n3,n3d)
endif
end




subroutine transpose12(p,pt,nx,nxd,ny,nyd,nz,nzd)
implicit none
real*8 p(nxd,nyd,nzd)
real*8 pt(nyd,nxd,nzd)

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


