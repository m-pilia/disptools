#include "disptools.cuh"
#include "jacobian.cuh"
#include "error.cuh"
#include "generate_displacement.cuh"
#include "reduction.cuh"

#define dfx_dx(f,x,y,z,idx) ((_((f), (x)+1,(y),  (z),   X) - _((f), (x)-1,(y),  (z),   X)) * (idx) * .5 + 1.0)
#define dfx_dy(f,x,y,z,idy) ((_((f), (x),  (y)+1,(z),   X) - _((f), (x),  (y)-1,(z),   X)) * (idy) * .5)
#define dfx_dz(f,x,y,z,idz) ((_((f), (x),  (y),  (z)+1, X) - _((f), (x),  (y),  (z)-1, X)) * (idz) * .5)
#define dfy_dx(f,x,y,z,idx) ((_((f), (x)+1,(y),  (z),   Y) - _((f), (x)-1,(y),  (z),   Y)) * (idx) * .5)
#define dfy_dy(f,x,y,z,idy) ((_((f), (x),  (y)+1,(z),   Y) - _((f), (x),  (y)-1,(z),   Y)) * (idy) * .5 + 1.0)
#define dfy_dz(f,x,y,z,idz) ((_((f), (x),  (y),  (z)+1, Y) - _((f), (x),  (y),  (z)-1, Y)) * (idz) * .5)
#define dfz_dx(f,x,y,z,idx) ((_((f), (x)+1,(y),  (z),   Z) - _((f), (x)-1,(y),  (z),   Z)) * (idx) * .5)
#define dfz_dy(f,x,y,z,idy) ((_((f), (x),  (y)+1,(z),   Z) - _((f), (x),  (y)-1,(z),   Z)) * (idy) * .5)
#define dfz_dz(f,x,y,z,idz) ((_((f), (x),  (y),  (z)+1, Z) - _((f), (x),  (y),  (z)-1, Z)) * (idz) * .5 + 1.0)


struct SquaredNorm
{
    __device__ __forceinline__
    FLOATING operator()(const FLOATING3 &a) const {
        return a.x * a.x + a.y * a.y + a.z * a.z;
    }
};


/*!
 * \brief Compute the gradient of the loss function.
 */
__global__ void gradient(
        const Image f,
        const Image g,
        const Image J,
        const FLOATING3 is,
        const Image voxel_error,
        const FLOATING delta,
        const FLOATING zeta
        )
{
    // +1 to not compute over the boundary voxels
    const size_t x = blockIdx.x * blockDim.x + threadIdx.x + 1;
    const size_t y = blockIdx.y * blockDim.y + threadIdx.y + 1;
    const size_t z = blockIdx.z * blockDim.z + threadIdx.z + 1;

    // Do not iterate over the voxels on the boundary
    if (x >= f.nx - 1 || y >= f.ny - 1 || z >= f.nz - 1) {
        return;
    }

    // Actual cost
    FLOATING error_xb = __(voxel_error, x-1, y,   z  );
    FLOATING error_xf = __(voxel_error, x+1, y,   z  );
    FLOATING error_yb = __(voxel_error, x,   y-1, z  );
    FLOATING error_yf = __(voxel_error, x,   y+1, z  );
    FLOATING error_zb = __(voxel_error, x,   y,   z-1);
    FLOATING error_zf = __(voxel_error, x,   y,   z+1);

    // Regularisation terms
    if (__(J, x-1, y, z) < delta) {
        error_xb += zeta * (__(J, x-1, y, z) - delta);
    }
    if (__(J, x+1, y, z) < delta) {
        error_xf += zeta * (__(J, x+1, y, z) - delta);
    }
    if (__(J, x, y-1, z) < delta) {
        error_yb += zeta * (__(J, x, y-1, z) - delta);
    }
    if (__(J, x, y+1, z) < delta) {
        error_yf += zeta * (__(J, x, y+1, z) - delta);
    }
    if (__(J, x, y, z-1) < delta) {
        error_zb += zeta * (__(J, x, y, z-1) - delta);
    }
    if (__(J, x, y, z+1) < delta) {
        error_zf += zeta * (__(J, x, y, z+1) - delta);
    }

    // dc/dx
    _(g, x, y, z, X) =
        (+is.x * (+error_xb * (dfy_dy(f, x-1, y, z, is.y) * dfz_dz(f, x-1, y, z, is.z) -
                               dfy_dz(f, x-1, y, z, is.z) * dfz_dy(f, x-1, y, z, is.y))
                  -error_xf * (dfy_dy(f, x+1, y, z, is.y) * dfz_dz(f, x+1, y, z, is.z) -
                               dfy_dz(f, x+1, y, z, is.z) * dfz_dy(f, x+1, y, z, is.y)))
         -is.y * (+error_yb * (dfy_dx(f, x, y-1, z, is.x) * dfz_dz(f, x, y-1, z, is.z) -
                               dfy_dz(f, x, y-1, z, is.z) * dfz_dx(f, x, y-1, z, is.x))
                  -error_yf * (dfy_dx(f, x, y+1, z, is.x) * dfz_dz(f, x, y+1, z, is.z) -
                               dfy_dz(f, x, y+1, z, is.z) * dfz_dx(f, x, y+1, z, is.x)))
         +is.z * (+error_zb * (dfy_dx(f, x, y, z-1, is.x) * dfz_dy(f, x, y, z-1, is.y) -
                               dfy_dy(f, x, y, z-1, is.y) * dfz_dx(f, x, y, z-1, is.x))
                  -error_zf * (dfy_dx(f, x, y, z+1, is.x) * dfz_dy(f, x, y, z+1, is.y) -
                               dfy_dy(f, x, y, z+1, is.y) * dfz_dx(f, x, y, z+1, is.x))));

    // dc/dy
    _(g, x, y, z, Y) =
        (-is.x * (+error_xb * (dfx_dy(f, x-1, y, z, is.y) * dfz_dz(f, x-1, y, z, is.z) -
                               dfx_dz(f, x-1, y, z, is.z) * dfz_dy(f, x-1, y, z, is.y))
                  -error_xf * (dfx_dy(f, x+1, y, z, is.y) * dfz_dz(f, x+1, y, z, is.z) -
                               dfx_dz(f, x+1, y, z, is.z) * dfz_dy(f, x+1, y, z, is.y)))
         +is.y * (+error_yb * (dfx_dx(f, x, y-1, z, is.x) * dfz_dz(f, x, y-1, z, is.z) -
                               dfx_dz(f, x, y-1, z, is.z) * dfz_dx(f, x, y-1, z, is.x))
                  -error_yf * (dfx_dx(f, x, y+1, z, is.x) * dfz_dz(f, x, y+1, z, is.z) -
                               dfx_dz(f, x, y+1, z, is.z) * dfz_dx(f, x, y+1, z, is.x)))
         -is.z * (+error_zb * (dfx_dx(f, x, y, z-1, is.x) * dfz_dy(f, x, y, z-1, is.y) -
                               dfx_dy(f, x, y, z-1, is.y) * dfz_dx(f, x, y, z-1, is.x))
                  -error_zf * (dfx_dx(f, x, y, z+1, is.x) * dfz_dy(f, x, y, z+1, is.y) -
                               dfx_dy(f, x, y, z+1, is.y) * dfz_dx(f, x, y, z+1, is.x))));

    // dc/dz
    _(g, x, y, z, Z) =
        (+is.x * (+error_xb * (dfx_dy(f, x-1, y, z, is.y) * dfy_dz(f, x-1, y, z, is.z) -
                               dfx_dz(f, x-1, y, z, is.z) * dfy_dy(f, x-1, y, z, is.y))
                  -error_xf * (dfx_dy(f, x+1, y, z, is.y) * dfy_dz(f, x+1, y, z, is.z) -
                               dfx_dz(f, x+1, y, z, is.z) * dfy_dy(f, x+1, y, z, is.y)))
         -is.y * (+error_yb * (dfx_dx(f, x, y-1, z, is.x) * dfy_dz(f, x, y-1, z, is.z) -
                               dfx_dz(f, x, y-1, z, is.z) * dfy_dx(f, x, y-1, z, is.x))
                  -error_yf * (dfx_dx(f, x, y+1, z, is.x) * dfy_dz(f, x, y+1, z, is.z) -
                               dfx_dz(f, x, y+1, z, is.z) * dfy_dx(f, x, y+1, z, is.x)))
         +is.z * (+error_zb * (dfx_dx(f, x, y, z-1, is.x) * dfy_dy(f, x, y, z-1, is.y) -
                               dfx_dy(f, x, y, z-1, is.y) * dfy_dx(f, x, y, z-1, is.x))
                  -error_zf * (dfx_dx(f, x, y, z+1, is.x) * dfy_dy(f, x, y, z+1, is.y) -
                               dfx_dy(f, x, y, z+1, is.y) * dfy_dx(f, x, y, z+1, is.x))));
}

/*!
 * \brief Move `old_field` along direction `g` with step size `eta`.
 */
__global__ void move_field(
        const Image old_field,
        const Image new_field,
        const Image g,
        const FLOATING eta
        )
{
    // +1 to not compute over the boundary voxels
    const size_t x = blockIdx.x * blockDim.x + threadIdx.x + 1;
    const size_t y = blockIdx.y * blockDim.y + threadIdx.y + 1;
    const size_t z = blockIdx.z * blockDim.z + threadIdx.z + 1;

    // Do not iterate over the voxels on the boundary
    if (x >= old_field.nx - 1 || y >= old_field.ny - 1 || z >= old_field.nz - 1) {
        return;
    }

    _(new_field, x, y, z, X) = _(old_field, x, y, z, X) - eta * _(g, x, y, z, X);
    _(new_field, x, y, z, Y) = _(old_field, x, y, z, Y) - eta * _(g, x, y, z, Y);
    _(new_field, x, y, z, Z) = _(old_field, x, y, z, Z) - eta * _(g, x, y, z, Z);
}


/*!
 * \brief Find a displacement field that realises the given Jacobian.
 *
 * Employ a greedy search, starting from an initial guess of the
 * displacement field (passed in the `field' argument). At each
 * iteration, compute the Jacobian of the current displacement field,
 * then correct the components of the field on each voxel.
 *
 * Use two couples of buffers, to store a copy of the displacement field
 * and its Jacobian at the current iteration, before and after the
 * correction. If the correction improves the result, then switch the
 * buffers and proceed with the next iteration, otherwise keep the
 * current displacement field.
 */
void generate_displacement_gradient_cuda(
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
    FLOATING eta,             /*!< Step length for the optimisation */
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

    const Sum<FLOATING> op;
    const SquaredNorm init;

    // Image size
    const size_t voxel_number = nx * ny * nz;
    const size_t image_size = voxel_number * sizeof (FLOATING);
    const size_t mask_size = voxel_number * sizeof (bool);

    // Use two buffers that are swapped
    unsigned old_buffer = 0, new_buffer = 1;

    // Wrap arrays in data strucutres
    Image J_ = new_gpu_image(1, nx, ny, nz, dx, dy, dz);
    cuda_safe_call(cudaMemcpy(J_.data, J, image_size, cudaMemcpyHostToDevice));
    Mask mask_ = new_gpu_mask(nx, ny, nz);
    cuda_safe_call(cudaMemcpy(mask_.data, mask, mask_size, cudaMemcpyHostToDevice));

    // Allocate memory for the Jacobian map of the moving field
    // Use two buffers
    Image J_field_[2] = {
        new_gpu_image(3, nx, ny, nz, dx, dy, dz),
        new_gpu_image(3, nx, ny, nz, dx, dy, dz),
    };

    // Allocate memory for the moving field
    // Use two buffers
    Image field_[2] = {
        new_gpu_image(3, nx, ny, nz, dx, dy, dz),
        new_gpu_image(3, nx, ny, nz, dx, dy, dz),
    };

    // Allocate memory for the voxel error term
    Image voxel_error = new_gpu_image(1, nx, ny, nz, dx, dy, dz);

    // Allocate memory for the gradient
    Image g = new_gpu_image(3, nx, ny, nz, dx, dy, dz);

    FLOATING2 last_error = {FLOATING_MIN, 0.0};
    FLOATING2 error = {FLOATING_MAX, 0.0};
    FLOATING g_norm_2 = 0.0;

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

    // Verbose feedback
    verbose_printf(true,
                   "Iteration %5ld:  "
                   "total error %6e  "
                   "max voxel error %6e  "
                   "eta %6e\n",
                   0l, last_error.x, last_error.y, eta);

    // Compute gradient
    gradient<<<grid_size, block_size>>>(field_[old_buffer],
                                        g,
                                        J_field_[old_buffer],
                                        inverse_spacing,
                                        voxel_error,
                                        delta,
                                        zeta
                                        );

    cuda_safe_call((reduce_array<FLOATING3, FLOATING, Sum<FLOATING>, SquaredNorm>
                                ((FLOATING3*) g.data,
                                 &g_norm_2,
                                 voxel_number,
                                 256,
                                 op,
                                 init,
                                 0.0
                                 )));

    cuda_check_error();
    if (disptools_error.error) {
        goto cleanup;
    }

    // Find an high initial eta
    do {
        eta *= alpha;

        // Update the moving displacement field
        move_field<<<grid_size, block_size>>>(field_[old_buffer],
                                              field_[new_buffer],
                                              g,
                                              eta
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

        // Armijo-Goldstein condition
    } while (eta < eta_max && error.x - last_error.x > -gamma * eta * g_norm_2);

    // One alpha in excess from the last iteration, one to compensate
    // for the increment before the first iteration
    eta /= alpha * alpha;

    // Recompute the initial error
    last_error = compute_error(block_size,
                               J_,
                               J_field_[old_buffer],
                               mask_,
                               tolerance,
                               voxel_error
                               );

    size_t it;
    for (it = 1; it <= it_max; ++it) {

        // Compute gradient
        gradient<<<grid_size, block_size>>>(field_[old_buffer],
                                            g,
                                            J_field_[old_buffer],
                                            inverse_spacing,
                                            voxel_error,
                                            delta,
                                            zeta
                                            );

        cuda_safe_call((reduce_array<FLOATING3, FLOATING, Sum<FLOATING>, SquaredNorm>
                                    ((FLOATING3*)g.data,
                                     &g_norm_2,
                                     voxel_number,
                                     256,
                                     op,
                                     init,
                                     0
                                     )));

        // Backtracking line search
        eta *= alpha;
        eta = eta > eta_max ? eta_max : eta;
        while (true) {

            // Update the moving displacement field
            move_field<<<grid_size, block_size>>>(field_[old_buffer],
                                                  field_[new_buffer],
                                                  g,
                                                  eta
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

            // Armijo-Goldstein condition
            const bool eta_good = eta >= iota;
            const bool ag_condition = error.x - last_error.x > -gamma * eta * g_norm_2;
            const bool strict_condition = strict && error.y > last_error.y;
            if (eta_good && (ag_condition || strict_condition)) {
                eta *= beta;
            }
            else {
                break;
            }
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

        if (eta < iota) {
            verbose_printf(true, "Error not decreasing, terminating.\n");
            break;
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
    delete_gpu_image(&g);
    delete_gpu_image(&voxel_error);
}
