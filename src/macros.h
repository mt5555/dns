#ifdef NDEBUG
#else
#define ASSERT(str,cond) if (.not. cond) call abort(str)
#endif

! this is the default:
#define TRANSPOSE_X_SPLIT_Z
! undef it to use TRANSPOSE_X_SPLIT_Y, which is usefull if nslabz=1



#define DX_AND_DXX 2
#define DX_ONLY 1

#define my_x    mpicoords(1)
#define my_y    mpicoords(2)
#define my_z    mpicoords(3)

! variable used to hold C style pointers.  should be at least 64 bits
#define CPOINTER real*8
