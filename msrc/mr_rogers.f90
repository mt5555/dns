program 3D_mr_rogers
use mod_mpi
use mod_ debug
 


call mpi_cart_rank(MPI_COMM_WORLD,me_cart,ierr2)
call mpi_comm_size(MPI_COMM_WORLD,live_procs,ierr3)

old_com=MPI_COMM_WORD
call mpi_cart_create(old_comm,ndims,dims,periodic,reorganize,  &
                     comm_3d,ierr4)






if(debug_mpi.ne.0) then 
ierr=ierr4+iedd3+ierr2
if(ierr>)) write(*,*) "PDNS ERROR in 3D_setup",ierr4,iedd3,ierr2
end if

end program 3D_mr_rogers

