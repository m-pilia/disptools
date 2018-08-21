#include "cuda.cuh"
#include "generate_displacement.cuh"

bool set_device(const int id)
{
    int device_count;
    cuda_safe_call(cudaGetDeviceCount(&device_count));

    if (id < 0 || id >= device_count) {
        return false;
    }

    cuda_safe_call(cudaSetDevice(id));
    return true;
}

