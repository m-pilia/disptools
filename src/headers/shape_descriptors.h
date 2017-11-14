#ifndef _SHAPE_DESCRIPTORS_H_DEFINED
#define _SHAPE_DESCRIPTORS_H_DEFINED

#include "field.h"

/*!
 * \brief Compute the cubeness of a binary image.
 *
 * @param image Array of boolean values with [z,y,x] indexing
 * @param nx    Image size along the x axis.
 * @param ny    Image size along the y axis.
 * @param nz    Image size along the z axis.
 * @param cx    X component of the centroid. 
 * @param cy    Y component of the centroid. 
 * @param cz    Z component of the centroid. 
 * @param sx    Image spacing along the x axis.
 * @param sy    Image spacing along the y axis.
 * @param sz    Image spacing along the z axis.
 *
 * @return A `FLOATING' value of cubeness in the interval \$[0,1]\$.
 */
FLOATING cubeness(
        const bool *image,
        const size_t nx,
        const size_t ny,
        const size_t nz,
        const FLOATING xc,
        const FLOATING yc,
        const FLOATING zc,
        const FLOATING sx,
        const FLOATING sy,
        const FLOATING sz
        );

/*!
 * \brief Compute the octahedroness of a binary image.
 *
 * @param image Array of boolean values with [z,y,x] indexing
 * @param nx    Image size along the x axis.
 * @param ny    Image size along the y axis.
 * @param nz    Image size along the z axis.
 * @param cx    X component of the centroid. 
 * @param cy    Y component of the centroid. 
 * @param cz    Z component of the centroid. 
 * @param sx    Image spacing along the x axis.
 * @param sy    Image spacing along the y axis.
 * @param sz    Image spacing along the z axis.
 *
 * @return A `FLOATING' value of octahedroness in the interval \$[0,1]\$.
 */
FLOATING octahedroness(
        const bool *image,
        const size_t nx,
        const size_t ny,
        const size_t nz,
        const FLOATING xc,
        const FLOATING yc,
        const FLOATING zc,
        const FLOATING sx,
        const FLOATING sy,
        const FLOATING sz
        );

#endif // _SHAPE_DESCRIPTORS_H_DEFINED
