#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <stdbool.h>

#include "headers/field.h"
#include "headers/displacement_field_gradient.h"
#include "headers/jacobian.h"
#include "headers/error.h"

#define dfx_dx(f,x,y,z,idx) ((_((f), (x)+1,(y),  (z),   X) - _((f), (x)-1,(y),  (z),   X)) * (idx) * .5 + 1.0)
#define dfx_dy(f,x,y,z,idy) ((_((f), (x),  (y)+1,(z),   X) - _((f), (x),  (y)-1,(z),   X)) * (idy) * .5)
#define dfx_dz(f,x,y,z,idz) ((_((f), (x),  (y),  (z)+1, X) - _((f), (x),  (y),  (z)-1, X)) * (idz) * .5)
#define dfy_dx(f,x,y,z,idx) ((_((f), (x)+1,(y),  (z),   Y) - _((f), (x)-1,(y),  (z),   Y)) * (idx) * .5)
#define dfy_dy(f,x,y,z,idy) ((_((f), (x),  (y)+1,(z),   Y) - _((f), (x),  (y)-1,(z),   Y)) * (idy) * .5 + 1.0)
#define dfy_dz(f,x,y,z,idz) ((_((f), (x),  (y),  (z)+1, Y) - _((f), (x),  (y),  (z)-1, Y)) * (idz) * .5)
#define dfz_dx(f,x,y,z,idx) ((_((f), (x)+1,(y),  (z),   Z) - _((f), (x)-1,(y),  (z),   Z)) * (idx) * .5)
#define dfz_dy(f,x,y,z,idy) ((_((f), (x),  (y)+1,(z),   Z) - _((f), (x),  (y)-1,(z),   Z)) * (idy) * .5)
#define dfz_dz(f,x,y,z,idz) ((_((f), (x),  (y),  (z)+1, Z) - _((f), (x),  (y),  (z)-1, Z)) * (idz) * .5 + 1.0)

/*!
 * \brief Compute the gradient of the loss function.
 */
static inline void gradient(
        const Image f,
        const Image g,
        const Image J,
        FLOATING *g_norm_2,
        const Image voxel_error,
        const FLOATING delta,
        const FLOATING zeta
        )
{
    // Do not iterate over the voxels on the boundary
    const size_t x_max = f.nx - 1;
    const size_t y_max = f.ny - 1;
    const size_t z_max = f.nz - 1;

    // Precompute the step for finite differences
    const FLOATING idx = 1.0 / f.dx;
    const FLOATING idy = 1.0 / f.dy;
    const FLOATING idz = 1.0 / f.dz;

    // Local variable for the squared norm
    FLOATING squared_norm = 0.0;

    #pragma omp parallel for reduction(+: squared_norm) collapse(3) schedule(static)
    for (size_t z = 1; z < z_max; ++z) {
        for (size_t y = 1; y < y_max; ++y) {
            for (size_t x = 1; x < x_max; ++x) {

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
                    (+idx * (+error_xb * (dfy_dy(f, x-1, y, z, idy) * dfz_dz(f, x-1, y, z, idz) -
                                          dfy_dz(f, x-1, y, z, idz) * dfz_dy(f, x-1, y, z, idy))
                             -error_xf * (dfy_dy(f, x+1, y, z, idy) * dfz_dz(f, x+1, y, z, idz) -
                                          dfy_dz(f, x+1, y, z, idz) * dfz_dy(f, x+1, y, z, idy)))
                     -idy * (+error_yb * (dfy_dx(f, x, y-1, z, idx) * dfz_dz(f, x, y-1, z, idz) -
                                          dfy_dz(f, x, y-1, z, idz) * dfz_dx(f, x, y-1, z, idx))
                             -error_yf * (dfy_dx(f, x, y+1, z, idx) * dfz_dz(f, x, y+1, z, idz) -
                                          dfy_dz(f, x, y+1, z, idz) * dfz_dx(f, x, y+1, z, idx)))
                     +idz * (+error_zb * (dfy_dx(f, x, y, z-1, idx) * dfz_dy(f, x, y, z-1, idy) -
                                          dfy_dy(f, x, y, z-1, idy) * dfz_dx(f, x, y, z-1, idx))
                             -error_zf * (dfy_dx(f, x, y, z+1, idx) * dfz_dy(f, x, y, z+1, idy) -
                                          dfy_dy(f, x, y, z+1, idy) * dfz_dx(f, x, y, z+1, idx))));

                // dc/dy
                _(g, x, y, z, Y) =
                    (-idx * (+error_xb * (dfx_dy(f, x-1, y, z, idy) * dfz_dz(f, x-1, y, z, idz) -
                                          dfx_dz(f, x-1, y, z, idz) * dfz_dy(f, x-1, y, z, idy))
                             -error_xf * (dfx_dy(f, x+1, y, z, idy) * dfz_dz(f, x+1, y, z, idz) -
                                          dfx_dz(f, x+1, y, z, idz) * dfz_dy(f, x+1, y, z, idy)))
                     +idy * (+error_yb * (dfx_dx(f, x, y-1, z, idx) * dfz_dz(f, x, y-1, z, idz) -
                                          dfx_dz(f, x, y-1, z, idz) * dfz_dx(f, x, y-1, z, idx))
                             -error_yf * (dfx_dx(f, x, y+1, z, idx) * dfz_dz(f, x, y+1, z, idz) -
                                          dfx_dz(f, x, y+1, z, idz) * dfz_dx(f, x, y+1, z, idx)))
                     -idz * (+error_zb * (dfx_dx(f, x, y, z-1, idx) * dfz_dy(f, x, y, z-1, idy) -
                                          dfx_dy(f, x, y, z-1, idy) * dfz_dx(f, x, y, z-1, idx))
                             -error_zf * (dfx_dx(f, x, y, z+1, idx) * dfz_dy(f, x, y, z+1, idy) -
                                          dfx_dy(f, x, y, z+1, idy) * dfz_dx(f, x, y, z+1, idx))));

                // dc/dz
                _(g, x, y, z, Z) =
                    (+idx * (+error_xb * (dfx_dy(f, x-1, y, z, idy) * dfy_dz(f, x-1, y, z, idz) -
                                          dfx_dz(f, x-1, y, z, idz) * dfy_dy(f, x-1, y, z, idy))
                             -error_xf * (dfx_dy(f, x+1, y, z, idy) * dfy_dz(f, x+1, y, z, idz) -
                                          dfx_dz(f, x+1, y, z, idz) * dfy_dy(f, x+1, y, z, idy)))
                     -idy * (+error_yb * (dfx_dx(f, x, y-1, z, idx) * dfy_dz(f, x, y-1, z, idz) -
                                          dfx_dz(f, x, y-1, z, idz) * dfy_dx(f, x, y-1, z, idx))
                             -error_yf * (dfx_dx(f, x, y+1, z, idx) * dfy_dz(f, x, y+1, z, idz) -
                                          dfx_dz(f, x, y+1, z, idz) * dfy_dx(f, x, y+1, z, idx)))
                     +idz * (+error_zb * (dfx_dx(f, x, y, z-1, idx) * dfy_dy(f, x, y, z-1, idy) -
                                          dfx_dy(f, x, y, z-1, idy) * dfy_dx(f, x, y, z-1, idx))
                             -error_zf * (dfx_dx(f, x, y, z+1, idx) * dfy_dy(f, x, y, z+1, idy) -
                                          dfx_dy(f, x, y, z+1, idy) * dfy_dx(f, x, y, z+1, idx))));

                // (dc/dx)^2 + (dc/dy)^2 + (dc/dz)^2
                squared_norm += _(g, x, y, z, X) * _(g, x, y, z, X) +
                                _(g, x, y, z, Y) * _(g, x, y, z, Y) +
                                _(g, x, y, z, Z) * _(g, x, y, z, Z);
            }
        }
    }

    // Return squared norm 
    *g_norm_2 = squared_norm;
}

/*!
 * \brief Move `old_field` along direction `g` with step size `eta`.
 */
static inline void move_field(
        const Image old_field,
        const Image new_field,
        const Image g,
        const FLOATING eta
        )
{
    #pragma omp parallel for collapse(3) schedule(static)
    for (size_t z = 0; z < new_field.nz; ++z) {
        for (size_t y = 0; y < new_field.ny; ++y) {
            for (size_t x = 0; x < new_field.nx; ++x) {
                _(new_field, x, y, z, X) = _(old_field, x, y, z, X) - eta * _(g, x, y, z, X);
                _(new_field, x, y, z, Y) = _(old_field, x, y, z, Y) - eta * _(g, x, y, z, Y);
                _(new_field, x, y, z, Z) = _(old_field, x, y, z, Z) - eta * _(g, x, y, z, Z);
            }
        }
    }
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
void generate_displacement_gradient(
        const size_t nx,               /*!< Width of the image  */
        const size_t ny,               /*!< Length of the image */
        const size_t nz,               /*!< Depth of the image  */
        const FLOATING dx,             /*!< x spacing */
        const FLOATING dy,             /*!< y spacing */
        const FLOATING dz,             /*!< z spacing */
        const FLOATING J[nz][ny][nx],  /*!< Target Jacobian */
        const bool mask[nz][ny][nx],   /*!< Body mask */
        const FLOATING epsilon,        /*!< Tolerance on the Jacobian per voxel */
        const FLOATING tolerance,      /*!< Jacobian tolerance on background */
        FLOATING eta,                  /*!< Initial step length for the optimisation */
        const FLOATING alpha,          /*!< Step length increase coefficient */
        const FLOATING beta,           /*!< Step length decrease coefficient */
        const FLOATING gamma,          /*!< Armijo-Goldstein parameter */
        const FLOATING delta,          /*!< Jacobian regularisation threshold */
        const FLOATING zeta,           /*!< Jacobian regularisation weight */
        const bool strict,             /*!< Always improve maximum voxel error */
        const size_t it_max,           /*!< Maximum number of iterations */
        FLOATING field[3][nz][ny][nx]  /*!< Resulting displacement field */
        )
{
    assert(alpha > 0.0 && "alpha must be positive");
    assert(beta > 0.0 && "beta must be positive");
    assert(gamma > 0.0 && "gamma must be positive");
    assert(delta > 0.0 && "delta must be positive");
    assert(zeta > 0.0 && "zeta must be positive");
    assert(eta > 0.0 && "eta must be positive");
    assert(tolerance >= 0.0 && "Tolerance must be positive");
    assert(epsilon > 0.0 && "Epsilon must be positive");

    verbose_printf(DISPTOOLS_DEBUG,
                   "%s\n"
                   "nx:        %lu\n"
                   "ny:        %lu\n"
                   "nz:        %lu\n"
                   "dx:        %f\n"
                   "dy:        %f\n"
                   "dz:        %f\n"
                   "alpha:     %f\n"
                   "beta:      %f\n"
                   "gamma:     %f\n"
                   "delta:     %f\n"
                   "epsilon:   %f\n"
                   "zeta:      %f\n"
                   "eta:       %f\n"
                   "tolerance: %f\n"
                   "strict:    %d\n"
                   "it_max:    %lu\n",
                   __func__,
                   nx, ny, nz,
                   dx, dy, dz,
                   alpha, beta, gamma, delta,
                   epsilon, zeta, eta,
                   tolerance,
                   strict,
                   it_max);

    // Image size
    const size_t voxel_number = nx * ny * nz;
    const size_t image_size = voxel_number * sizeof (FLOATING);

    // Use two buffers that are swapped
    unsigned old_buffer = 0, new_buffer = 1;

    // Wrap arrays in data strucutres
    Image J_ = {1, nx, ny, nz, dx, dy, dz, (FLOATING*) J};
    Mask mask_ = {nx, ny, nz, (bool*) mask};

    // Create tolerance map
    Image tolerance_map = create_tolerance_map(mask_, tolerance);

    // Allocate memory for the Jacobian map of the moving field
    // Use two buffers
    Image J_field_[2] = {
        new_image(3, nx, ny, nz, dx, dy, dz),
        new_image(3, nx, ny, nz, dx, dy, dz),
    };

    // Allocate memory for the moving field
    // Use two buffers
    Image field_[2] = {
        new_image(3, nx, ny, nz, dx, dy, dz),
        new_image(3, nx, ny, nz, dx, dy, dz),
    };

    // Allocate memory for the voxel error term
    Image voxel_error = new_image(1, nx, ny, nz, dx, dy, dz);

    // Allocate memory for the gradient
    Image g = new_image(3, nx, ny, nz, dx, dy, dz);

    FLOATING last_error = DBL_MAX, error = DBL_MAX;
    FLOATING max_voxel_error = DBL_MAX, last_max_voxel_error = DBL_MAX;
    FLOATING g_norm_2 = 0.0;

    // Copy initial guess in the buffer
    memcpy(field_[old_buffer].data, field, 3 * image_size);

    // Compute the error of the initial guess
    jacobian(field_[old_buffer], J_field_[old_buffer]);
    last_error = compute_error(J_,
                               J_field_[old_buffer],
                               mask_,
                               tolerance_map,
                               voxel_error,
                               &max_voxel_error
                               );

    // Verbose feedback
    verbose_printf(true,
                   "Iteration %5ld:  "
                   "total error %6e  "
                   "max voxel error %6e  "
                   "eta %6e\n",
                   0l, last_error, max_voxel_error, eta);

    // Compute gradient
    gradient(field_[old_buffer],
             g,
             J_field_[old_buffer],
             &g_norm_2,
             voxel_error,
             delta,
             zeta
             );

    // Find an high initial eta
    do {
        eta *= alpha;

        // Update the moving displacement field
        move_field(field_[old_buffer],
                   field_[new_buffer],
                   g,
                   eta
                   );

        // Compute the Jacobian map of the moving displacement field
        jacobian(field_[new_buffer], J_field_[new_buffer]);

        // Compute the error of the moving field
        error = compute_error(J_,
                              J_field_[new_buffer],
                              mask_,
                              tolerance_map,
                              voxel_error,
                              &max_voxel_error
                              );

        // Armijo-Goldstein condition
    } while (error - last_error > -gamma * eta * g_norm_2);

    // One alpha in excess from the last iteration, one to compensate
    // for the increment before the first iteration
    eta /= alpha * alpha;

    // Recompute the initial error
    last_error = compute_error(J_,
                               J_field_[old_buffer],
                               mask_,
                               tolerance_map,
                               voxel_error,
                               &max_voxel_error
                               );

    size_t it;
    for (it = 1; it <= it_max; ++it) {

        // Compute gradient
        gradient(field_[old_buffer],
                 g,
                 J_field_[old_buffer],
                 &g_norm_2,
                 voxel_error,
                 delta,
                 zeta
                 );

        // Backtracking line search
        eta *= alpha;
        while (true) {

            // Update the moving displacement field
            move_field(field_[old_buffer],
                       field_[new_buffer],
                       g,
                       eta
                       );

            // Compute the Jacobian map of the moving displacement field
            jacobian(field_[new_buffer], J_field_[new_buffer]);

            // Compute the error of the moving field
            error = compute_error(J_,
                                  J_field_[new_buffer],
                                  mask_,
                                  tolerance_map,
                                  voxel_error,
                                  &max_voxel_error
                                  );

            // Armijo-Goldstein condition
            const bool eta_good = eta >= 1e-9;
            const bool ag_condition = error - last_error > -gamma * eta * g_norm_2;
            const bool strict_condition = strict && max_voxel_error > last_max_voxel_error;
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
                       it, error, max_voxel_error, eta);

        // Stopping conditions

        if (!isnormal(error)) {
            verbose_printf(true, "Terminating: error exploded.\n");
            break;
        }

        if (eta < 1e-9) {
            verbose_printf(true, "Error not decreasing, terminating.\n");
            break;
        }

		if (error / last_error > 0.999999) {
            verbose_printf(true, "Error not decreasing, terminating.\n");
            break;
		}

        if (!isnormal(max_voxel_error)) {
            verbose_printf(true, "Terminating: voxel error exploded.\n");
            break;
        }

        if (max_voxel_error < epsilon) {
            verbose_printf(true, "Terminating: reached desired tolerance.\n");
            break;
        }

        // Save error and swap the buffers
        last_error = error;
        last_max_voxel_error = max_voxel_error;
        XOR_SWAP(old_buffer, new_buffer);
    }

    verbose_printf(it == it_max, "Terminating: reached maximum number of iterations.\n");

    // Copy result for the caller
    memcpy(field, field_[old_buffer].data, 3 * image_size);

    // Release buffers
    delete_image(&field_[0]);
    delete_image(&field_[1]);
    delete_image(&J_field_[0]);
    delete_image(&J_field_[1]);
    delete_image(&g);
    delete_image(&voxel_error);
    delete_image(&tolerance_map);
}
