#define u_index 1
#define v_index 2
#define w_index 3
#define rho_index 4
#define e_index 5

#ifdef NDEBUG
#else
#define ASSERT(str,cond) if (.not. cond) call abort(str)
#endif



#define DX_ONLY     1
#define DX_AND_DXX  2

