#ifndef _VTK_H_DEFINED
#define _VTK_H_DEFINED

#include "field.h"

/*!
 * \brief Write to file a vector field in VTK (STRUCTURED_POINTS) format.
 */
void write_vtk(
        const char *filename,                /*!< Filename */
        const Image image
        );

/*!
 * \brief Read from file a vector field in VTK (STRUCTURED_POINTS) format.
 */
Image read_vtk(
        const char *filename /*!< Filename */
        );

#endif // _VTK_H_DEFINED
