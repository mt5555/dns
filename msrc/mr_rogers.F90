program 3D_mr_rogers
use mod-mpi
use mod-debug




!call mpi_cart_create(old_comm,ndims,dims,periodic,reorganize,  &
!                     comm_3d,ierr4)

call mpi_cart_coords(comm_3d,me_global,ndims,p_cart_coords,ierr3)

coords=p_cart_coords
do ix=1,dims(1)
coords(1)=ix-1

call mpi_cart_rank(comm_3d,coords,you_cart,ierr4)
you_cart_ranks_x(ix)=you_cart

end do


do iy=1,dims(2)
coords(2)=iy-1

call mpi_cart_rank(comm_3d,coords,you_cart,ierr4)
you_cart_ranks_y(iy)=you_cart

end do

do iz=1,dims(3)
coords(3)=iz-1

call mpi_cart_rank(comm_3d,coords,you_cart,ierr4)
you_cart_ranks_y(iz)=you_cart

end do



if(debug_mpi.ne.0) then 
ierr=ierr4+iedd3+ierr2
if(ierr>)) write(*,*) "PDNS ERROR in 3D_setup",ierr4,iedd3,ierr2
end if

write(*,*) me_global,"Cart Coords", p_cart_coords

end program 3D_mr_rogers

