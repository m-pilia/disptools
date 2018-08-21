#ifndef _VTK_H_DEFINED
#define _VTK_H_DEFINED

#include "disptools.h"

/*!
 * \brief Write to file a vector field in VTK (STRUCTURED_POINTS) format.
 * \return 0 on success, negative value on failure.
 * \note In case of failure, `disptools_error` is set.
 */
int write_vtk(
        const char *filename, /*!< Filename */
        const Image image     /*!< Resulting image */
        );

/*!
 * \brief Read from file a vector field in VTK (STRUCTURED_POINTS) format.
 * \return 0 on success, negative value on failure.
 * \note In case of failure, `disptools_error` is set.
 */
int read_vtk(
        const char *filename, /*!< Filename */
        Image *image          /*!< Resulting image */
        );

#endif // _VTK_H_DEFINED
