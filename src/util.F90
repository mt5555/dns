subroutine abort(message)
implicit none
character*(*) message
character*15 :: pre="ASSERT FAILURE "
call print_message(pre // message)
stop
end subroutine


subroutine print_message(message)
implicit none
character*(*) message
write(*,'(a)') message
end subroutine


subroutine transpose12(p,pt,swap_index,nx,nxd,ny,nyd,nz,nzd)
implicit none
real*8 p(nxd,nyd,nzd)
real*8 pt(nyd,nxd,nzd)
integer nx,ny,nz,nxd,nyd,nzd,swap_index

integer i,j,k

do k=1,nz
do j=1,ny
do i=1,nx
   pt(j,i,k)=p(i,j,k)
enddo
enddo
enddo

if (swap_index<>0) then
   i=ny
   ny=nx
   nx=i

   i=nyd
   nyd=nxd
   nxd=i
endif

end



subroutine transpose13(p,pt,swap_index,nx,nxd,ny,nyd,nz,nzd)
implicit none
real*8 p(nxd,nyd,nzd)
real*8 pt(nzd,nyd,nxd)
integer nx,ny,nz,nxd,nyd,nzd,swap_index

integer i,j,k

do i=1,nx
do j=1,ny
do k=1,nz
   pt(k,j,i)=p(i,j,k)
enddo
enddo
enddo

if (swap_index<>0) then
   i=nz
   nz=nx
   nx=i

   i=nzd
   nzd=nxd
   nxd=i
endif

end







