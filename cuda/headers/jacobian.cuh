#ifndef JACOBIAN_CUH_INCLUDED
#define JACOBIAN_CUH_INCLUDED

#include "disptools.h"
#include "cuda.cuh"

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
        );

#endif /* ifndef JACOBIAN_CUH_INCLUDED */
