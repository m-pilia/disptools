#ifndef __JACOBIAN_H_DEFINED
#define __JACOBIAN_H_DEFINED

#include "disptools.h"

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
        );

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
        );


/*!
 * \brief Compute the Jacobian determinant associated to a 
 *        tridimensional displacement field.
 */
void jacobian(
        Image f,
        Image J
        );

#endif // __JACOBIAN_H_DEFINED
