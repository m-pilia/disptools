#ifndef ERROR_CUH_INCLUDED
#define ERROR_CUH_INCLUDED

#include "disptools.h"
#include "cuda.cuh"


/*!
 * \brief Compute a tolerance map over the volume.
 */
__global__ void create_tolerance_map(
        const Mask mask,
        const FLOATING tolerance,
        Image tolerance_map
        );

 
/*!
 * \brief Compute the error on the Jacobian of the moving field.
 */
FLOATING2 compute_error(
        dim3 block_size,          /*!< CUDA block size */
        const Image J,            /*!< Target Jacobian */
        const Image J_field,      /*!< Current Jacobian */
        const Mask mask,          /*!< Body mask */
        const FLOATING tolerance, /*!< Jacobian tolerance on background */
        Image voxel_error         /*!< Error on the Jacobian */
        );

#endif /* ifndef ERROR_CUH_INCLUDED */
