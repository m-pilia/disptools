#include "jacobian.h"
#include "jacobian_macros.h"

/*!
 * \brief Dynamic version of the function.
 */
void jacobian_dynamic(
        const size_t nx,   /*!< Image width */
        const size_t ny,   /*!< Image length */
        const size_t nz,   /*!< Image depth */
        const FLOATING dx, /*!< x spacing */
        const FLOATING dy, /*!< y spacing */
        const FLOATING dz, /*!< z spacing */
        const FLOATING *f, /*!< Vector field */
        FLOATING *J        /*!< Resulting Jacobian */
        )
{
    Image input = {3, nx, ny, nz, dx, dy, dz, (FLOATING*) f};
    Image output = {1, nx, ny, nz, dx, dy, dz, J};
    jacobian(input, output);
}

/*!
 * \brief Regularise a Jacobian map.
 *
 * Replace with `epsilon` all values of the Jacobian map that
 * fall below `epsilon`.
 */
void regularise(
        const size_t nx, /*!< Width of the image */
        const size_t ny, /*!< Length of the image */
        const size_t nz, /*!< Depth of the image */
        FLOATING *J,     /*!< Jacobian map */
        FLOATING epsilon /*!< Minimum value allowed */
        )
{
#ifdef __GNUC__
    #pragma omp parallel for
    for (size_t z = 0; z < nz; ++z) {
        for (size_t y = 0; y < ny; ++y) {
            for (size_t x = 0; x < nx; ++x) {
#else // MSVC 15 does not support OpenMP > 2.0
    int z;
    #pragma omp parallel for
    for (z = 0; z < nz; ++z) {
        for (size_t y = 0; y < ny; ++y) {
            for (size_t x = 0; x < nx; ++x) {
#endif
                if (J[z*ny*nx + y*nx + x] < epsilon) {
                    J[z*ny*nx + y*nx + x] = epsilon;
                }
            }
        }
    }
}

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
void jacobian(
        Image f,
        Image J
        )
{
    // Precompute the step for finite differences
    const FLOATING idx = 1.0 / f.dx;
    const FLOATING idy = 1.0 / f.dy;
    const FLOATING idz = 1.0 / f.dz;

    const size_t x_start = ORDER_PD / 2 ;
    const size_t y_start = ORDER_PD / 2 ;
    const size_t z_start = ORDER_PD / 2 ;
    const size_t x_stop = f.nx - ORDER_PD / 2 ;
    const size_t y_stop = f.ny - ORDER_PD / 2 ;
    const size_t z_stop = f.nz - ORDER_PD / 2 ;

    // Inner voxels
#ifdef __GNUC__
    #pragma omp parallel for collapse(3) schedule(static)
    for (size_t z = z_start; z < z_stop; ++z) {
        for (size_t y = y_start; y < y_stop; ++y) {
            for (size_t x = x_start; x < x_stop; ++x) {
#else // MSVC 15 does not support OpenMP > 2.0
    int z;
    #pragma omp parallel for
    for (z = z_start; z < z_stop; ++z) {
        for (size_t y = y_start; y < y_stop; ++y) {
            for (size_t x = x_start; x < x_stop; ++x) {
#endif
                // Approximate partial derivatives with central differences
                #if ORDER_PD == 2
                    __(J, x, y, z) = Jacobian_2(f, x, y, z, idx, idy, idz);
                #elif ORDER_PD == 4
                    __(J, x, y, z) = Jacobian_4(f, x, y, z, idx, idy, idz);
                #else
                    #error "Unsupported order for partial derivatives"
                #endif
            }
        }
    }

    // Corner cases 
    const int nx = f.nx-1;
    const int ny = f.ny-1;
    const int nz = f.nz-1;

    // 8 corners
    __(J, 0,  0,  0 ) = det3j(dfx_dx_f(0 , 0 , 0 ), dfx_dy_f(0 , 0 , 0 ), dfx_dz_f(0 , 0 , 0 ),
                              dfy_dx_f(0 , 0 , 0 ), dfy_dy_f(0 , 0 , 0 ), dfy_dz_f(0 , 0 , 0 ),
                              dfz_dx_f(0 , 0 , 0 ), dfz_dy_f(0 , 0 , 0 ), dfz_dz_f(0 , 0 , 0 ));
    __(J, 0,  0,  nz) = det3j(dfx_dx_f(0 , 0 , nz), dfx_dy_f(0 , 0 , nz), dfx_dz_b(0 , 0 , nz),
                              dfy_dx_f(0 , 0 , nz), dfy_dy_f(0 , 0 , nz), dfy_dz_b(0 , 0 , nz),
                              dfz_dx_f(0 , 0 , nz), dfz_dy_f(0 , 0 , nz), dfz_dz_b(0 , 0 , nz));
    __(J, 0,  ny, 0 ) = det3j(dfx_dx_f(0 , ny, 0 ), dfx_dy_b(0 , ny, 0 ), dfx_dz_f(0 , ny, 0 ),
                              dfy_dx_f(0 , ny, 0 ), dfy_dy_b(0 , ny, 0 ), dfy_dz_f(0 , ny, 0 ),
                              dfz_dx_f(0 , ny, 0 ), dfz_dy_b(0 , ny, 0 ), dfz_dz_f(0 , ny, 0 ));
    __(J, 0,  ny, nz) = det3j(dfx_dx_f(0 , ny, nz), dfx_dy_b(0 , ny, nz), dfx_dz_b(0 , ny, nz),
                              dfy_dx_f(0 , ny, nz), dfy_dy_b(0 , ny, nz), dfy_dz_b(0 , ny, nz),
                              dfz_dx_f(0 , ny, nz), dfz_dy_b(0 , ny, nz), dfz_dz_b(0 , ny, nz));
    __(J, nx, 0,  0 ) = det3j(dfx_dx_b(nx, 0 , 0 ), dfx_dy_f(nx, 0 , 0 ), dfx_dz_f(nx, 0 , 0 ),
                              dfy_dx_b(nx, 0 , 0 ), dfy_dy_f(nx, 0 , 0 ), dfy_dz_f(nx, 0 , 0 ),
                              dfz_dx_b(nx, 0 , 0 ), dfz_dy_f(nx, 0 , 0 ), dfz_dz_f(nx, 0 , 0 ));
    __(J, nx, 0,  nz) = det3j(dfx_dx_b(nx, 0 , nz), dfx_dy_f(nx, 0 , nz), dfx_dz_b(nx, 0 , nz),
                              dfy_dx_b(nx, 0 , nz), dfy_dy_f(nx, 0 , nz), dfy_dz_b(nx, 0 , nz),
                              dfz_dx_b(nx, 0 , nz), dfz_dy_f(nx, 0 , nz), dfz_dz_b(nx, 0 , nz));
    __(J, nx, ny, 0 ) = det3j(dfx_dx_b(nx, ny, 0 ), dfx_dy_b(nx, ny, 0 ), dfx_dz_f(nx, ny, 0 ),
                              dfy_dx_b(nx, ny, 0 ), dfy_dy_b(nx, ny, 0 ), dfy_dz_f(nx, ny, 0 ),
                              dfz_dx_b(nx, ny, 0 ), dfz_dy_b(nx, ny, 0 ), dfz_dz_f(nx, ny, 0 ));
    __(J, nx, ny, nz) = det3j(dfx_dx_b(nx, ny, nz), dfx_dy_b(nx, ny, nz), dfx_dz_b(nx, ny, nz),
                              dfy_dx_b(nx, ny, nz), dfy_dy_b(nx, ny, nz), dfy_dz_b(nx, ny, nz),
                              dfz_dx_b(nx, ny, nz), dfz_dy_b(nx, ny, nz), dfz_dz_b(nx, ny, nz));

    // 4 edges along the x axis
    for (size_t x = 1; x < f.nx-1; ++x) {
        __(J, x, 0,  0 ) = det3j(dfx_dx_c(x, 0 , 0 ), dfx_dy_f(x, 0 , 0 ), dfx_dz_f(x, 0 , 0 ),
                                 dfy_dx_c(x, 0 , 0 ), dfy_dy_f(x, 0 , 0 ), dfy_dz_f(x, 0 , 0 ),
                                 dfz_dx_c(x, 0 , 0 ), dfz_dy_f(x, 0 , 0 ), dfz_dz_f(x, 0 , 0 ));
        __(J, x, 0,  nz) = det3j(dfx_dx_c(x, 0 , nz), dfx_dy_f(x, 0 , nz), dfx_dz_b(x, 0 , nz),
                                 dfy_dx_c(x, 0 , nz), dfy_dy_f(x, 0 , nz), dfy_dz_b(x, 0 , nz),
                                 dfz_dx_c(x, 0 , nz), dfz_dy_f(x, 0 , nz), dfz_dz_b(x, 0 , nz));
        __(J, x, ny, 0 ) = det3j(dfx_dx_c(x, ny, 0 ), dfx_dy_b(x, ny, 0 ), dfx_dz_f(x, ny, 0 ),
                                 dfy_dx_c(x, ny, 0 ), dfy_dy_b(x, ny, 0 ), dfy_dz_f(x, ny, 0 ),
                                 dfz_dx_c(x, ny, 0 ), dfz_dy_b(x, ny, 0 ), dfz_dz_f(x, ny, 0 ));
        __(J, x, ny, nz) = det3j(dfx_dx_c(x, ny, nz), dfx_dy_b(x, ny, nz), dfx_dz_b(x, ny, nz),
                                 dfy_dx_c(x, ny, nz), dfy_dy_b(x, ny, nz), dfy_dz_b(x, ny, nz),
                                 dfz_dx_c(x, ny, nz), dfz_dy_b(x, ny, nz), dfz_dz_b(x, ny, nz));
    }

    for (size_t y = 1; y < f.ny-1; ++y) {

        // 4 edges along the y axis
        __(J, 0 , y, 0 ) = det3j(dfx_dx_f(0 , y, 0 ), dfx_dy_c(0 , y, 0 ), dfx_dz_f(0 , y, 0 ),
                                 dfy_dx_f(0 , y, 0 ), dfy_dy_c(0 , y, 0 ), dfy_dz_f(0 , y, 0 ),
                                 dfz_dx_f(0 , y, 0 ), dfz_dy_c(0 , y, 0 ), dfz_dz_f(0 , y, 0 ));
        __(J, 0 , y, nz) = det3j(dfx_dx_f(0 , y, nz), dfx_dy_c(0 , y, nz), dfx_dz_b(0 , y, nz),
                                 dfy_dx_f(0 , y, nz), dfy_dy_c(0 , y, nz), dfy_dz_b(0 , y, nz),
                                 dfz_dx_f(0 , y, nz), dfz_dy_c(0 , y, nz), dfz_dz_b(0 , y, nz));
        __(J, nx, y, 0 ) = det3j(dfx_dx_b(nx, y, 0 ), dfx_dy_c(nx, y, 0 ), dfx_dz_f(nx, y, 0 ),
                                 dfy_dx_b(nx, y, 0 ), dfy_dy_c(nx, y, 0 ), dfy_dz_f(nx, y, 0 ),
                                 dfz_dx_b(nx, y, 0 ), dfz_dy_c(nx, y, 0 ), dfz_dz_f(nx, y, 0 ));
        __(J, nx, y, nz) = det3j(dfx_dx_b(nx, y, nz), dfx_dy_c(nx, y, nz), dfx_dz_b(nx, y, nz),
                                 dfy_dx_b(nx, y, nz), dfy_dy_c(nx, y, nz), dfy_dz_b(nx, y, nz),
                                 dfz_dx_b(nx, y, nz), dfz_dy_c(nx, y, nz), dfz_dz_b(nx, y, nz));

        // 2 faces along the xy plane
        for (size_t x = 1; x < f.nx-1; ++x) {
            __(J, x, y, 0 ) = det3j(dfx_dx_c(x, y, 0 ), dfx_dy_c(x, y, 0 ), dfx_dz_f(x, y, 0 ),
                                    dfy_dx_c(x, y, 0 ), dfy_dy_c(x, y, 0 ), dfy_dz_f(x, y, 0 ),
                                    dfz_dx_c(x, y, 0 ), dfz_dy_c(x, y, 0 ), dfz_dz_f(x, y, 0 ));
            __(J, x, y, nz) = det3j(dfx_dx_c(x, y, nz), dfx_dy_c(x, y, nz), dfx_dz_b(x, y, nz),
                                    dfy_dx_c(x, y, nz), dfy_dy_c(x, y, nz), dfy_dz_b(x, y, nz),
                                    dfz_dx_c(x, y, nz), dfz_dy_c(x, y, nz), dfz_dz_b(x, y, nz));
        }
    }

    for (size_t z = 1; z < f.nz-1; ++z) {

        // 4 edges along the z axis
        __(J, 0 , 0 , z) = det3j(dfx_dx_f(0 , 0 , z), dfx_dy_f(0 , 0 , z), dfx_dz_c(0 , 0 , z),
                                 dfy_dx_f(0 , 0 , z), dfy_dy_f(0 , 0 , z), dfy_dz_c(0 , 0 , z),
                                 dfz_dx_f(0 , 0 , z), dfz_dy_f(0 , 0 , z), dfz_dz_c(0 , 0 , z));
        __(J, 0 , ny, z) = det3j(dfx_dx_f(0 , ny, z), dfx_dy_b(0 , ny, z), dfx_dz_c(0 , ny, z),
                                 dfy_dx_f(0 , ny, z), dfy_dy_b(0 , ny, z), dfy_dz_c(0 , ny, z),
                                 dfz_dx_f(0 , ny, z), dfz_dy_b(0 , ny, z), dfz_dz_c(0 , ny, z));
        __(J, nx, 0 , z) = det3j(dfx_dx_b(nx, 0 , z), dfx_dy_f(nx, 0 , z), dfx_dz_c(nx, 0 , z),
                                 dfy_dx_b(nx, 0 , z), dfy_dy_f(nx, 0 , z), dfy_dz_c(nx, 0 , z),
                                 dfz_dx_b(nx, 0 , z), dfz_dy_f(nx, 0 , z), dfz_dz_c(nx, 0 , z));
        __(J, nx, ny, z) = det3j(dfx_dx_b(nx, ny, z), dfx_dy_b(nx, ny, z), dfx_dz_c(nx, ny, z),
                                 dfy_dx_b(nx, ny, z), dfy_dy_b(nx, ny, z), dfy_dz_c(nx, ny, z),
                                 dfz_dx_b(nx, ny, z), dfz_dy_b(nx, ny, z), dfz_dz_c(nx, ny, z));

        // 2 faces along the xz plane
        for (size_t x = 1; x < f.nx-1; ++x) {
            __(J, x, 0 , z) = det3j(dfx_dx_c(x, 0 , z), dfx_dy_f(x, 0 , z), dfx_dz_c(x, 0 , z),
                                    dfy_dx_c(x, 0 , z), dfy_dy_f(x, 0 , z), dfy_dz_c(x, 0 , z),
                                    dfz_dx_c(x, 0 , z), dfz_dy_f(x, 0 , z), dfz_dz_c(x, 0 , z));
            __(J, x, ny, z) = det3j(dfx_dx_c(x, ny, z), dfx_dy_b(x, ny, z), dfx_dz_c(x, ny, z),
                                    dfy_dx_c(x, ny, z), dfy_dy_b(x, ny, z), dfy_dz_c(x, ny, z),
                                    dfz_dx_c(x, ny, z), dfz_dy_b(x, ny, z), dfz_dz_c(x, ny, z));
        }

        // 2 faces along the yz plane
        for (size_t y = 1; y < f.ny-1; ++y) {
            __(J, 0 , y, z) = det3j(dfx_dx_f(0 , y, z), dfx_dy_c(0 , y, z), dfx_dz_c(0 , y, z),
                                    dfy_dx_f(0 , y, z), dfy_dy_c(0 , y, z), dfy_dz_c(0 , y, z),
                                    dfz_dx_f(0 , y, z), dfz_dy_c(0 , y, z), dfz_dz_c(0 , y, z));
            __(J, nx, y, z) = det3j(dfx_dx_b(nx, y, z), dfx_dy_c(nx, y, z), dfx_dz_c(nx, y, z),
                                    dfy_dx_b(nx, y, z), dfy_dy_c(nx, y, z), dfy_dz_c(nx, y, z),
                                    dfz_dx_b(nx, y, z), dfz_dy_c(nx, y, z), dfz_dz_c(nx, y, z));
        }
    }

    // Second order central difference in the penultimate voxels
    #if ORDER_PD == 4
    for (size_t z = 1; z < f.nz-1; ++z) {
        for (size_t y = 1; y < f.ny-1; ++y) {
            __(J, 1, y, z) = Jacobian_2(f, 1, y, z, idx, idy, idz);
            __(J, f.nx-2, y, z) = Jacobian_2(f, f.nx-2, y, z, idx, idy, idz);
        }
        for (size_t x = 1; x < f.nx-1; ++x) {
            __(J, x, 1, z) = Jacobian_2(f, x, 1, z, idx, idy, idz);
            __(J, x, f.ny-2, z) = Jacobian_2(f, x, f.ny-2, z, idx, idy, idz);
        }
    }
    for (size_t y = 1; y < f.ny-1; ++y) {
        for (size_t x = 1; x < f.nx-1; ++x) {
            __(J, x, y, 1) = Jacobian_2(f, x, y, 1, idx, idy, idz);
            __(J, x, y, f.nz-2) = Jacobian_2(f, x, y, f.nz-2, idx, idy, idz);
        }
    }
    #endif
}

