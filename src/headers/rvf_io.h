#ifndef _RVF_H_DEFINED
#define _RVF_H_DEFINED

#include "disptools.h"

/*!
 * \brief Write to file a vector field in rvf format.
 */
int write_rvf(
        const char *filename, /*!< Filename */
        const Image image     /*!< Image */
        );

/*!
 * \brief Read from file a vector field in rvf format.
 *
 */
int read_rvf(
        const char *filename, /*!< Filename. */
        Image *image          /*!< Image */ 
        );

#endif // _RVF_H_DEFINED
