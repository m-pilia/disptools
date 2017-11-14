#ifndef _RVF_H_DEFINED
#define _RVF_H_DEFINED

#include "field.h"

/*!
 * \brief Write to file a vector field in rvf format.
 */
void write_rvf(
        const char *filename,               /*!< Filename */
        const Image image
        );

/*!
 * \brief Read from file a vector field in rvf format.
 *
 */
Image read_rvf(
        const char *filename /*!< Filename. */
        );

#endif // _RVF_H_DEFINED
