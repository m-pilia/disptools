#include "cuda.cuh"
#include "error.cuh"
#include "reduction.cuh"


struct ErrorInit
{
    __device__ __forceinline__
    FLOATING2 operator()(const FLOATING &a) const {
        return {a * a, ABS_FLOATING(a)};
    }
};

struct ErrorFold
{
    __device__ __forceinline__
    FLOATING2 operator()(const FLOATING2 &a, const FLOATING2 &b) const {
        return {a.x + b.x, MAX_FLOATING(a.y, b.y)};
    }
};


/*!
 * \brief Compute a tolerance map over the volume.
 *
 * The tolerance is zero within the body volume and along the boundary,
 * while it is set to `tolerance' on the remaining background space.
 * Zero tolerance on the boundary prevents the vector field from 
 * flowing outside the image space in the contours of the regions
 * where the body touches the boundary.
 */
__global__ void create_tolerance_map(
        const Mask mask,
        const FLOATING tolerance,
        Image tolerance_map
        )
{
    const size_t x = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t y = blockIdx.y * blockDim.y + threadIdx.y;
    const size_t z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x >= mask.nx || y >= mask.ny || z >= mask.nz) {
        return;
    }
    
    if (__(mask, x, y, z) || 
            x == 0         || y == 0         || z == 0        || 
            x == mask.nx-1 || y == mask.ny-1 || z == mask.nz-1) {
        __(tolerance_map, x, y, z) = 0.0;
    }
    else {
        __(tolerance_map, x, y, z) = tolerance;
    }
}

 
/*!
 * \brief Kernel to compute the error over the Jacobian map.
 */
__global__ void compute_voxel_error(
        const Image J,            /*!< Target Jacobian */
        const Image J_field,      /*!< Current Jacobian */
        const Mask mask,          /*!< Body mask */
        const FLOATING tolerance, /*!< Jacobian tolerance on background */
        const Image voxel_error   /*!< Error on the Jacobian */
        )
{
    // Do not iterate over the voxels on the boundary
    // because in this implementation the Jacobian is not
    // computed along the boundary
    const size_t x = blockIdx.x * blockDim.x + threadIdx.x + 1;
    const size_t y = blockIdx.y * blockDim.y + threadIdx.y + 1;
    const size_t z = blockIdx.z * blockDim.z + threadIdx.z + 1;

    if (x >= J.nx - 1 || y >= J.ny - 1 || z >= J.nz - 1) {
        return;
    }

    // Compute the error on the voxel
    FLOATING error = __(J_field, x, y, z) - __(J, x, y, z);

    // Full error inside the body mask
    if (__(mask, x, y, z)) {
        __(voxel_error, x, y, z) = error;
    }

    // Tolerate some error in the background
    else {
        const FLOATING abs_error = fabs(error);
        __(voxel_error, x, y, z) =
            abs_error < tolerance ?
            0.0 :
            copysign(abs_error - tolerance, error);
    }
}


/*!
 * \brief Compute the error on the Jacobian of the moving field.
 */
FLOATING2 compute_error(
        dim3 img_block_size,      /*!< CUDA block size */
        const Image J,            /*!< Target Jacobian */
        const Image J_field,      /*!< Current Jacobian */
        const Mask mask,          /*!< Body mask */
        const FLOATING tolerance, /*!< Jacobian tolerance on background */
        Image voxel_error         /*!< Error on the Jacobian */
        )
{
    dim3 img_grid_size = {
        (unsigned int) (J.nx + img_block_size.x - 1) / img_block_size.x,
        (unsigned int) (J.ny + img_block_size.y - 1) / img_block_size.y,
        (unsigned int) (J.nz + img_block_size.z - 1) / img_block_size.z,
    };

    compute_voxel_error<<<img_grid_size, img_block_size>>>(J, J_field, mask, tolerance, voxel_error);

    const int voxel_count = voxel_error.nx * voxel_error.ny * voxel_error.nz;
    const int block_size = 256;
    FLOATING2 result;
    struct ErrorFold op;
    struct ErrorInit init;
    const FLOATING2 null = {0.0, 0.0};

    cuda_safe_call((reduce_array<FLOATING, FLOATING2, ErrorFold, ErrorInit>
                                (voxel_error.data, &result, voxel_count, block_size, op, init, null)));

    if (disptools_error.error) {
        return {-1.0, -1.0};
    }

    return result;
}

