#ifdef USE_STKFFT
#include "fftstk_interface.F90"
#elif (defined USE_SGIFFT)
#include "fftsgi_interface.F90"
#elif (defined USE_CPQFFT)
#include "fftcpq_interface.F90"
#elif (defined USE_FFT99)
#include "fft99_interface.F90"
#else
#include "fft99_interface.F90"
#endif




