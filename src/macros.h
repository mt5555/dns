#define u_index 1
#define v_index 2
#define w_index 3
#define rho_index 4
#define e_index 5

#ifdef NDEBUG
#else
#define ASSERT(str,cond) if (.not. cond) call abort(str)
#endif



#define DX_AND_DXX 2
#define DX_ONLY 1
#define DXONLY 1

#define ncpu_x  mpidims(1)
#define ncpu_y  mpidims(2)
#define ncpu_z  mpidims(3)
#define my_x    mpicoords(1)
#define my_y    mpicoords(2)
#define my_z    mpicoords(3)
