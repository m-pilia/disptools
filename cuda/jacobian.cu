#include "jacobian.cuh"
#include "jacobian_macros.h"

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

    #if ORDER_PD == 2
        __(J, x, y, z) = Jacobian_2(f, x, y, z, is.x, is.y, is.z);
    #elif ORDER_PD == 4
        if (x == 1 || x == f.nx - 2 ||
            y == 1 || y == f.ny - 2 ||
            z == 1 || z == f.nz - 2)
        {
            __(J, x, y, z) = Jacobian_2(f, x, y, z, is.x, is.y, is.z);

        }
        else
        {
            __(J, x, y, z) = Jacobian_4(f, x, y, z, is.x, is.y, is.z);
        }
    #else
        #error "Unsupported order for partial derivatives"
    #endif
}

