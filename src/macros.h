#ifdef MPI_UNDERSCORE
#define mpi_reduce mpi_reduce_
#define mpi_allreduce mpi_allreduce_
#define mpi_bcast mpi_bcast_
#define mpi_abort mpi_abort_
#define mpi_barrier mpi_barrier_
#define mpi_cart_rank mpi_cart_rank_
#define mpi_irecv mpi_irecv_
#define mpi_isend  mpi_isend_
#define mpi_cart_create mpi_cart_create_
#define mpi_cart_get mpi_cart_get_
#define mpi_comm_free mpi_comm_free_
#define mpi_comm_rank mpi_comm_rank_
#define mpi_comm_size mpi_comm_size_
#define mpi_comm_split mpi_comm_split_
#define mpi_finalize mpi_finalize_
#define mpi_init mpi_init_
#define mpi_irecv mpi_irecv_
#define mpi_isent mpi_isend_
#define mpi_waitall mpi_waitall_
#define mpi_wtime mpi_wtime_
#endif


#ifdef NDEBUG
#define ASSERT(str,cond) 
#else
#define ASSERT(str,cond) if (.not. cond) call abort(str)
#endif


#define DX_AND_DXX 2
#define DX_ONLY 1

#define my_x    mpicoords(1)
#define my_y    mpicoords(2)
#define my_z    mpicoords(3)

! variable used to hold C style pointers.  should be at least 64 bits
#define CPOINTER real*8

! equations
#define NS_UVW  0
#define SHALLOW 1
#define NS_PSIVOR 2
#define CNS       3

! numerical methods
#define FOURIER        0
#define FOURTH_ORDER    1

! types of boundary conditions
! these boundaries are not 'real', and are just handled by the
! parallel ghost cell update:
#define PERIODIC 0
#define REFLECT  1
#define REFLECT_ODD        2
#define INTERNAL           3

! these are the real boundaries and must be
! treated with code outside of the ghost cell update:
! to check if a boundary is a 'real' boundary, use for example:  
!  REALBOUNDARY(bdy_x1)
#define REALBOUNDARY(bdy)    (bdy>=100)
#define INFLOW0            100
#define INFLOW0_ONESIDED   101






