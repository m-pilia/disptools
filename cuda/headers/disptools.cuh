#ifndef FIELD_CUH_INCLUDED
#define FIELD_CUH_INCLUDED value

#include "disptools.h"

Image new_gpu_image(
        const size_t nd,
        const size_t nx,
        const size_t ny,
        const size_t nz,
        const FLOATING dx,
        const FLOATING dy,
        const FLOATING dz
        );

Mask new_gpu_mask(
        const size_t nx,
        const size_t ny,
        const size_t nz
        );

void delete_gpu_image(Image *img);

void delete_gpu_mask(Mask *mask);

#endif /* ifndef  */
