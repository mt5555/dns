#ifdef NDEBUG
#else
#define ASSERT(str,cond) if (.not. cond) call abort(str)
#endif



#define DX_AND_DXX 2
#define DX_ONLY 1

#define my_x    mpicoords(1)
#define my_y    mpicoords(2)
#define my_z    mpicoords(3)
