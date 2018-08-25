#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cuda.cuh"
#include "disptools.cuh"
#include "disptools.h"

Image new_gpu_image(
        const size_t nd,
        const size_t nx,
        const size_t ny,
        const size_t nz,
        const FLOATING dx,
        const FLOATING dy,
        const FLOATING dz
        )
{
    Image img = {nd, nx, ny, nz, dx, dy, dz, NULL};
    cuda_safe_call(cudaMalloc((void**) &img.data, nd * nx * ny * nz * sizeof (FLOATING)));
    if (disptools_has_error()) {
        return img;
    }
    cuda_safe_call(cudaMemset(img.data, 0, nd * nx * ny * nz * sizeof (FLOATING)));
    return img;
}

Mask new_gpu_mask(
        const size_t nx,
        const size_t ny,
        const size_t nz
        )
{
    Mask mask = {nx, ny, nz, NULL};
    cuda_safe_call(cudaMalloc((void**) &mask.data, nx * ny * nz * sizeof (bool)));
    if (disptools_has_error()) {
        return mask;
    }
    cuda_safe_call(cudaMemset(mask.data, 0, nx * ny * nz * sizeof (bool)));
    return mask;
}

void delete_gpu_image(Image *img)
{
    if (img->data) {
        cuda_safe_call(cudaFree(img->data));
        img->data = NULL;
    }
}

void delete_gpu_mask(Mask *mask)
{
    if (mask->data) {
        cuda_safe_call(cudaFree(mask->data));
        mask->data = NULL;
    }
}

