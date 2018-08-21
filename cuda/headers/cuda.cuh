#ifndef CUDA_CUH_INCLUDED
#define CUDA_CUH_INCLUDED

#include <cstdio>
#include <stdexcept>

#include "disptools.h"


#define SQUARE_FLOATING(a) ((a) * (a))


#if FLOAT_SIZE == 32
    #define FLOATING2 float2
    #define FLOATING3 float3
    #define MAX_FLOATING fmaxf
    #define ABS_FLOATING fabsf
#elif FLOAT_SIZE == 64
    #define FLOATING2 double2
    #define FLOATING3 double3
    #define MAX_FLOATING fmax
    #define ABS_FLOATING fabs
#else
    #error "Invalid FLOAT_SIZE. Accepted values are 32 and 64."
#endif
 

#ifdef DISPTOOLS_CUDA_ERROR_CHECK
    #define cuda_safe_call(e) __cuda_safe_call(e, #e, __FILE__, __LINE__)
    #define cuda_check_error() __cuda_check_error(__FILE__, __LINE__)
    #define cuda_return_error(e) \
        if (cudaSuccess != e) { \
            return e; \
        }
#else
    #define cuda_safe_call(e) e
    #define cuda_check_error()
    #define cuda_return_error(e) e
#endif
 

inline void __cuda_safe_call(cudaError e, const char *call, const char *file, const int line)
{
    if (cudaSuccess != e) {
        char msg[1024]; 
        sprintf(msg, "cuda_safe_call(%s) failed at %s:%i : %s\n",
                call, file, line, cudaGetErrorString(e));
        DISPTOOLS_SET_ERROR(true, msg);
    }
 
    return;
}
 

inline void __cuda_check_error(const char *file, const int line)
{
    cudaError e = cudaGetLastError();
    if (cudaSuccess != e) {
        char msg[1024];
        sprintf(msg, "cuda_check_error() failed at %s:%i : %s\n",
                file, line, cudaGetErrorString(e));
        DISPTOOLS_SET_ERROR(true, msg);
    }
 
#ifdef DISPTOOLS_CUDA_ERROR_CHECK_SYNC
    // More careful checking, with an extra performance cost
    e = cudaDeviceSynchronize();
    if(cudaSuccess != e)
    {
        char msg[1024];
        sprintf(msg, "cuda_check_error() with sync failed at %s:%i : %s\n",
                file, line, cudaGetErrorString(e));
        DISPTOOLS_SET_ERROR(true, msg);
    }
#endif
 
    return;
}


#endif /* ifndef CUDA_CUH_INCLUDED */
