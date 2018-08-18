#include "headers/jacobian.h"

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
