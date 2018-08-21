#include "disptools.cuh"
#include "jacobian.cuh"
#include "error.cuh"
#include "generate_displacement.cuh"

/*!
 * \brief Update the displacement field according to the error
 *        on the Jacobian map.
 *
 * Each component of the displacement field is updated independently.
 * For each component, check the error on the two neighbours. If the
 * error is positive (negative) on the preceding (following) voxel
 * along the direction, then the value of the field must be decreased,
 * vice versa for errors of the opposite sign.
 */
__global__ void greedy_step(
        const Image old_field,   /*!< Current displacement field */
        const Image new_field,   /*!< New displacement field */
        const Image voxel_error, /*!< Error on the Jacobian */
        const FLOATING3 stepl    /*!< Step length parameter */
        )
{
    // +1 to not compute over the boundary voxels
    const size_t x = blockIdx.x * blockDim.x + threadIdx.x + 1;
    const size_t y = blockIdx.y * blockDim.y + threadIdx.y + 1;
    const size_t z = blockIdx.z * blockDim.z + threadIdx.z + 1;

    // Do not iterate over the voxels on the boundary
    if (x >= old_field.nx - 1 ||
        y >= old_field.ny - 1 ||
        z >= old_field.nz - 1) {
        return;
    }

    // Update the displacement vectors in each voxel according to the error
    const FLOATING delta_x = stepl.x * (__(voxel_error, x+1, y,   z  ) -
                                        __(voxel_error, x-1, y,   z  ));
    const FLOATING delta_y = stepl.y * (__(voxel_error, x,   y+1, z  ) -
                                        __(voxel_error, x,   y-1, z  ));
    const FLOATING delta_z = stepl.z * (__(voxel_error, x,   y,   z+1) -
                                        __(voxel_error, x,   y,   z-1));

    _(new_field, x, y, z, X) = _(old_field, x, y, z, X) + delta_x;
    _(new_field, x, y, z, Y) = _(old_field, x, y, z, Y) + delta_y;
    _(new_field, x, y, z, Z) = _(old_field, x, y, z, Z) + delta_z;
}

/*!
 * \brief Find a displacement field that realises the given Jacobian.
 *
 * Employ a greedy search, starting from an initial guess of the
 * displacement field (passed in the `field' argument). At each
 * iteration, compute the Jacobian of the current displacement field,
 * then correct the components of the field on each voxel according to
 * the error on the neighbours. Each component is corrected
 * independently, according to the two neighbours along that direction.
 *
 * Use two couples of buffers, to store a copy of the displacement field
 * and its Jacobian at the current iteration, before and after the
 * correction. If the correction improves the result, then switch the
 * buffers and proceed with the next iteration, otherwise keep the
 * current displacement field.
 */
void generate_displacement_greedy_cuda(
        const size_t nx,          /*!< Width of the image  */
        const size_t ny,          /*!< Length of the image */
        const size_t nz,          /*!< Depth of the image  */
        const FLOATING dx,        /*!< x spacing */
        const FLOATING dy,        /*!< y spacing */
        const FLOATING dz,        /*!< z spacing */
        const FLOATING *J,        /*!< Target Jacobian */
        const bool *mask,         /*!< Body mask */
        const FLOATING epsilon,   /*!< Tolerance on the Jacobian per voxel */
        const FLOATING tolerance, /*!< Jacobian tolerance on background */
        FLOATING eta,             /*!< Initial step length for the optimisation */
        const FLOATING eta_max,   /*!< Maximum step length allowed */
        const FLOATING alpha,     /*!< Step length increase coefficient */
        const FLOATING beta,      /*!< Step length decrease coefficient */
        const FLOATING gamma,     /*!< Armijo-Goldstein parameter */
        const FLOATING delta,     /*!< Jacobian regularisation threshold */
        const FLOATING zeta,      /*!< Jacobian regularisation weight */
        const FLOATING theta,     /*!< Termination condition based on improvement */
        const FLOATING iota,      /*!< Termination condition based on eta */
        const bool strict,        /*!< Always improve maximum voxel error */
        const size_t it_max,      /*!< Maximum number of iterations */
        FLOATING *field           /*!< Resulting displacement field */
        )
{
    ASSERT_PARAMETERS;
    disptools_error.error = false;

    dim3 block_size = {8, 8, 8};
    dim3 grid_size = {
        (unsigned int) ceilf((float) nx / block_size.x),
        (unsigned int) ceilf((float) ny / block_size.y),
        (unsigned int) ceilf((float) nz / block_size.z),
    };

    FLOATING3 inverse_spacing = {1 / dx, 1 / dy, 1 / dz};

    // Image size
    const size_t voxel_number = nx * ny * nz;
    const size_t image_size = voxel_number * sizeof (FLOATING);
    const size_t mask_size = voxel_number * sizeof (bool);

    // Use two buffers that are swapped
    unsigned old_buffer = 0, new_buffer = 1;

    // Step length
    const FLOATING3 stepl = {.5f * eta * dx, .5f * eta * dy, .5f * eta * dz};

    // Wrap arrays in data strucutres
    Image J_ = new_gpu_image(1, nx, ny, nz, dx, dy, dz);
    cuda_safe_call(cudaMemcpy(J_.data, J, image_size, cudaMemcpyHostToDevice));
    Mask mask_ = new_gpu_mask(nx, ny, nz);
    cuda_safe_call(cudaMemcpy(mask_.data, mask, mask_size, cudaMemcpyHostToDevice));

    // Allocate memory for the Jacobian map of the moving field
    // Use two buffers
    Image J_field_[2] = {
        new_gpu_image(1, nx, ny, nz, dx, dy, dz),
        new_gpu_image(1, nx, ny, nz, dx, dy, dz),
    };

    // Allocate memory for the moving field
    // Use two buffers
    Image field_[2] = {
        new_gpu_image(3, nx, ny, nz, dx, dy, dz),
        new_gpu_image(3, nx, ny, nz, dx, dy, dz),
    };

    // Allocate memory for the voxel error term
    Image voxel_error = new_gpu_image(1, nx, ny, nz, dx, dy, dz);

    FLOATING2 last_error = {FLOATING_MIN, 0.0};
    FLOATING2 error = {FLOATING_MAX, 0.0};

    if (disptools_error.error) {
        goto cleanup;
    }

    // Copy initial guess in the buffer
    cuda_safe_call(cudaMemcpy(field_[old_buffer].data,
                              field,
                              3 * image_size,
                              cudaMemcpyHostToDevice));

    // Compute the error of the initial guess
    jacobian<<<grid_size, block_size>>>(field_[old_buffer],
                                        inverse_spacing,
                                        J_field_[old_buffer]);

    last_error = compute_error(block_size,
                               J_,
                               J_field_[old_buffer],
                               mask_,
                               tolerance,
                               voxel_error
                               );

    cuda_check_error();
    if (disptools_error.error) {
        goto cleanup;
    }

    size_t it;
    for (it = 1; it <= it_max; ++it) {

        // Update the moving displacement field
        greedy_step<<<grid_size, block_size>>>(field_[old_buffer],
                                               field_[new_buffer],
                                               voxel_error,
                                               stepl
                                               );

        // Compute the Jacobian map of the moving displacement field
        jacobian<<<grid_size, block_size>>>(field_[new_buffer],
                                            inverse_spacing,
                                            J_field_[new_buffer]);

        // Compute the error of the moving field
        error = compute_error(block_size,
                              J_,
                              J_field_[new_buffer],
                              mask_,
                              tolerance,
                              voxel_error
                              );

        cuda_check_error();
        if (disptools_error.error) {
            goto cleanup;
        }

        // Verbose feedback
        verbose_printf(true,
                       "Iteration %5ld:  "
                       "total error %6e  "
                       "max voxel error %6e  "
                       "eta %6e\n",
                       it, error.x, error.y, eta);

        // Stopping conditions

        if (!isnormal(error.x)) {
            verbose_printf(true, "Terminating: error exploded.\n");
            break;
        }

        if (error.x >= last_error.x) {
            // Try to reduce the step size
            eta *= beta;

            // Terminate if eta is too small
            if (eta < iota) {
                verbose_printf(true, "Error not decreasing, terminating.\n");
                break;
            }

            // Otherwise, repeat the last iteration with the new eta
            --it;
            verbose_printf(true, "Error not decreasing, "
                                 "reducing step size to %.4e\n", eta);
            continue;
        }

        if (1.0 - error.x / last_error.x < theta) {
            verbose_printf(true, "Error not decreasing, terminating.\n");
            break;
        }

        if (!isnormal(error.y)) {
            verbose_printf(true, "Terminating: voxel error exploded.\n");
            break;
        }

        if (error.y < epsilon) {
            verbose_printf(true, "Terminating: reached desired tolerance.\n");
            break;
        }

        // Save error and swap the buffers
        last_error = error;
        XOR_SWAP(old_buffer, new_buffer);
    }

    verbose_printf(it == it_max, "Terminating: reached maximum number of iterations.\n");

    // Copy result for the caller
    cuda_safe_call(cudaMemcpy(field,
                              field_[old_buffer].data,
                              3 * image_size,
                              cudaMemcpyDeviceToHost));

cleanup:
    // Release buffers
    delete_gpu_image(&J_);
    delete_gpu_mask(&mask_);
    delete_gpu_image(&field_[0]);
    delete_gpu_image(&field_[1]);
    delete_gpu_image(&J_field_[0]);
    delete_gpu_image(&J_field_[1]);
    delete_gpu_image(&voxel_error);
}

