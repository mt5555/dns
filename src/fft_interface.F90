#ifdef USE_STKFFT
#include "fftstk_interface.F90"
#elif (defined USE_FFT99 || defined USE_SGIFFT)
#include "fft99_interface.F90"
#elif (defined USE_CPQFFT)
#include "fftcpq_interface.F90"
#elif 
#include "fft99_interface.F90"
#endif




