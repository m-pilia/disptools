#include "jacobian.cuh"

// Central partial differences with second order accuracy 
#define dfx_dx_2(f, x, y, z, idx) \
         ((_((f), (x)+1, (y),  (z),   X) - _((f), (x)-1, (y),  (z),   X)) * (idx) * .5)
#define dfx_dy_2(f, x, y, z, idy) \
         ((_((f), (x),  (y)+1, (z),   X) - _((f), (x),  (y)-1, (z),   X)) * (idy) * .5)
#define dfx_dz_2(f, x, y, z, idz) \
         ((_((f), (x),  (y),   (z)+1, X) - _((f), (x),  (y),   (z)-1, X)) * (idz) * .5)
#define dfy_dx_2(f, x, y, z, idx) \
         ((_((f), (x)+1, (y),  (z),   Y) - _((f), (x)-1, (y),  (z),   Y)) * (idx) * .5)
#define dfy_dy_2(f, x, y, z, idy) \
         ((_((f), (x),  (y)+1, (z),   Y) - _((f), (x),  (y)-1, (z),   Y)) * (idy) * .5)
#define dfy_dz_2(f, x, y, z, idz) \
         ((_((f), (x),  (y),   (z)+1, Y) - _((f), (x),  (y),   (z)-1, Y)) * (idz) * .5)
#define dfz_dx_2(f, x, y, z, idx) \
         ((_((f), (x)+1, (y),  (z),   Z) - _((f), (x)-1, (y),  (z),   Z)) * (idx) * .5)
#define dfz_dy_2(f, x, y, z, idy) \
         ((_((f), (x),  (y)+1, (z),   Z) - _((f), (x),  (y)-1, (z),   Z)) * (idy) * .5)
#define dfz_dz_2(f, x, y, z, idz) \
         ((_((f), (x),  (y),   (z)+1, Z) - _((f), (x),  (y),   (z)-1, Z)) * (idz) * .5)

// Approximate the Jacobian
// Add 1.0 along the diagonal to add the identity component 
// of the transform to the displacement field.

// Second order accuracy
#define Jacobian_2(f, x, y, z, idx, idy, idz) \
   det3j(dfx_dx_2((f),(x),(y),(z),(idx)), dfx_dy_2((f),(x),(y),(z),(idy)), dfx_dz_2((f),(x),(y),(z),(idz)),  \
         dfy_dx_2((f),(x),(y),(z),(idx)), dfy_dy_2((f),(x),(y),(z),(idy)), dfy_dz_2((f),(x),(y),(z),(idz)),  \
         dfz_dx_2((f),(x),(y),(z),(idx)), dfz_dy_2((f),(x),(y),(z),(idy)), dfz_dz_2((f),(x),(y),(z),(idz)));


/*!
 * \brief Compute the Jacobian determinant associated to a 
 *        tridimensional displacement field.
 *
 * The Jacobian is a scalar map defined on a grid with the same
 * size of the image. It is computed by approximating partial
 * derivatives with first order central differences.
 *
 * Since the displacement field denotes the difference between the transform
 * that warps the image and the identity, the Jacobian matrix of the identity is
 * added to the actual Jacobian matrix of the displacement field, in order to
 * obtain the Jacobian of the whole transformation.
 */
__global__ void jacobian(
        const Image f,
        const FLOATING3 is,
        Image J
        )
{
    // +1 to not compute along the boundary
    const size_t x = blockIdx.x * blockDim.x + threadIdx.x + 1;
    const size_t y = blockIdx.y * blockDim.y + threadIdx.y + 1;
    const size_t z = blockIdx.z * blockDim.z + threadIdx.z + 1;

    if (x >= f.nx - 1 || y >= f.ny - 1 || z >= f.nz - 1) {
        return;
    }

    __(J, x, y, z) = Jacobian_2(f, x, y, z, is.x, is.y, is.z);
}


#undef dfx_dx_2
#undef dfx_dy_2
#undef dfx_dz_2
#undef dfy_dx_2
#undef dfy_dy_2
#undef dfy_dz_2
#undef dfz_dx_2
#undef dfz_dy_2
#undef dfz_dz_2

#undef Jacobian_2

