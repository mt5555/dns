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




#if 0
subroutine transpose_to_y(p,pt,n1,n1d,n2,n2d,n3,n3d)
real*8 p(nx,ny,nz)
real*8 pt(g_ny,nz_2d,nx)

!nz_2d = (nz2-nz1+1)/ncpu_y
nslaby = (ny2-ny1+1)


do iproc=0,ncpu_y-1
   if (iproc==myproc_y) then
      do i=nx1,nx2
         do j=ny1,ny2
         do k=1,nslab
            pt(j+iproc*nslaby,k,i)=p(i,j,k+iproc*nz_2d)
         enddo
         enddo
      enddo
   else
      l=0
      do i=nx1,nx2
         do j=ny1,ny2
         do k=1,nslab
            l=l+1
            sendbuf(l)=p(i,j,k+iproc*nz_2d)
         enddo
         enddo
      enddo

!     send buffer to (myproc_x,iproc,mproc_z)
!     rec  buffer from (myproc_x,iproc,mproc_z)

      l=0
      do i=nx1,nx2
         do j=ny1,ny2
         do k=1,nslab
            l=l+1
            pt(j+iproc*nslaby,k,i)=recbuf(l)
         enddo
         enddo
      enddo
   endif
enddo

end subroutine
#endif

