#include <stdlib.h>
#include <math.h>

#include "headers/shape_descriptors.h"

/*!
 * Here the minimum is not computed, the optimum is assumed to be
 * already aligned with the axis. To recover the minimum, the function
 * can be used as loss for a stochastic optimisation process managed at
 * higher level by the caller.
 *
 * Octahedroness is defined in [1] as
 *
 * \[
 *
 *      \frac{3}{8}
 *      \cdot
 *      \frac{\mu_{000}(O)^{\frac{4}{3}}}
 *           {\min_{\theta, \phi} \iiint_{R_{\theta,\phi}(O)}
 *                 \max\{|x|,|y|,|z|\} \dif x \dif y \dif z} .
 * \]
 *
 * [1] Martinez-Ortiz, Carlos A. "2D and 3D shape descriptors." (2010).
 */
FLOATING cubeness(
        const bool * restrict image,
        const size_t nx,
        const size_t ny,
        const size_t nz,
        const FLOATING xc,
        const FLOATING yc,
        const FLOATING zc,
        const FLOATING sx,
        const FLOATING sy,
        const FLOATING sz
        )
{
    (void) nz;

    double integral = 0.0;
    double volume = 0.0;
    const double dv = sx * sy * sz;

    #pragma omp parallel for reduction(+: volume, integral)
    for (size_t z = 0; z < nz; ++z) {
        for (size_t y = 0; y < ny; ++y) {
            for (size_t x = 0; x < nx; ++x) {
                if (image[z*ny*nx + y*nx + x]) {
                    volume += dv;
                    integral += max(sx * abs(x - xc),
                                    sy * abs(y - yc),
                                    sz * abs(z - zc));
                }
            }
        }
        
    }

    return (FLOATING) ((3.0 / 8.0) * pow(volume, 4.0 / 3.0) / integral);
}

/*!
 * Here the minimum is not computed, the optimum is assumed to be
 * already aligned with the axis. To recover the minimum, the function
 * can be used as loss for a stochastic optimisation process managed at
 * higher level by the caller.
 *
 * Octahedroness is defined in [1] as
 *
 * \[
 *      \frac{7}{9}
 *      \left( \frac{3}{4} \right)^{\frac{4}{3}}
 *      \cdot
 *      \frac{\mu_{000}(O)^{\frac{4}{3}}}
 *           {\min_{\theta, \phi} \iiint_{R_{\theta,\phi}(O)}
 *                 \left( |x| + |y| + |z| \right) \dif x \dif y \dif z} .
 * \]
 *
 * [1] Martinez-Ortiz, Carlos A. "2D and 3D shape descriptors." (2010).
 */
FLOATING octahedroness(
        const bool * restrict image,
        const size_t nx,
        const size_t ny,
        const size_t nz,
        const FLOATING xc,
        const FLOATING yc,
        const FLOATING zc,
        const FLOATING sx,
        const FLOATING sy,
        const FLOATING sz
        )
{
    (void) nz;

    double integral = 0.0;
    double volume = 0.0;
    const double dv = sx * sy * sz;

    #pragma omp parallel for reduction(+: volume, integral)
    for (size_t z = 0; z < nz; ++z) {
        for (size_t y = 0; y < ny; ++y) {
            for (size_t x = 0; x < nx; ++x) {
                if (image[z*ny*nx + y*nx + x]) {
                    volume += dv;
                    integral += sx * abs(x - xc) +
                                sy * abs(y - yc) +
                                sz * abs(z - zc);
                }
            }
        }
        
    }

    return (FLOATING) (pow((3.0 / 4.0) * volume, 4.0 / 3.0) / integral);
}

