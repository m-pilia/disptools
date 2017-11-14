#ifndef __JACOBIAN_H_DEFINED
#define __JACOBIAN_H_DEFINED

#include "field.h"

#define TWO_THIRDS  0.666666666666666666666666666
#define ONE_TWELFTH 0.083333333333333333333333333

// Forward differences with second order accuracy
#define dfx_dx_2f(f, x, y, z, idx) \
         ((-1.5 * _((f), (x),   (y), (z),  X) \
           +2.0 * _((f), (x)+1, (y), (z),  X) \
           -0.5 * _((f), (x)+2, (y), (z),  X) \
           ) * (idx))
#define dfx_dy_2f(f, x, y, z, idy) \
         ((-1.5 * _((f), (x), (y),   (z),  X) \
           +2.0 * _((f), (x), (y)+1, (z),  X) \
           -0.5 * _((f), (x), (y)+2, (z),  X) \
           ) * (idy))
#define dfx_dz_2f(f, x, y, z, idz) \
         ((-1.5 * _((f), (x), (y), (z),    X) \
           +2.0 * _((f), (x), (y), (z)+1,  X) \
           -0.5 * _((f), (x), (y), (z)+2,  X) \
           ) * (idz))
#define dfy_dx_2f(f, x, y, z, idx) \
         ((-1.5 * _((f), (x),   (y), (z),  Y) \
           +2.0 * _((f), (x)+1, (y), (z),  Y) \
           -0.5 * _((f), (x)+2, (y), (z),  Y) \
           ) * (idx))
#define dfy_dy_2f(f, x, y, z, idy) \
         ((-1.5 * _((f), (x), (y),   (z),  Y) \
           +2.0 * _((f), (x), (y)+1, (z),  Y) \
           -0.5 * _((f), (x), (y)+2, (z),  Y) \
           ) * (idy))
#define dfy_dz_2f(f, x, y, z, idz) \
         ((-1.5 * _((f), (x), (y), (z),    Y) \
           +2.0 * _((f), (x), (y), (z)+1,  Y) \
           -0.5 * _((f), (x), (y), (z)+2,  Y) \
           ) * (idz))
#define dfz_dx_2f(f, x, y, z, idx) \
         ((-1.5 * _((f), (x),   (y), (z),  Z) \
           +2.0 * _((f), (x)+1, (y), (z),  Z) \
           -0.5 * _((f), (x)+2, (y), (z),  Z) \
           ) * (idx))
#define dfz_dy_2f(f, x, y, z, idy) \
         ((-1.5 * _((f), (x), (y),   (z),  Z) \
           +2.0 * _((f), (x), (y)+1, (z),  Z) \
           -0.5 * _((f), (x), (y)+2, (z),  Z) \
           ) * (idy))
#define dfz_dz_2f(f, x, y, z, idz) \
         ((-1.5 * _((f), (x), (y), (z),    Z) \
           +2.0 * _((f), (x), (y), (z)+1,  Z) \
           -0.5 * _((f), (x), (y), (z)+2,  Z) \
           ) * (idz))

// Backward differences with second order accuracy
#define dfx_dx_2b(f, x, y, z, idx) \
         ((+1.5 * _((f), (x),   (y), (z),  X) \
           -2.0 * _((f), (x)-1, (y), (z),  X) \
           +0.5 * _((f), (x)-2, (y), (z),  X) \
           ) * (idx))
#define dfx_dy_2b(f, x, y, z, idy) \
         ((+1.5 * _((f), (x), (y),   (z),  X) \
           -2.0 * _((f), (x), (y)-1, (z),  X) \
           +0.5 * _((f), (x), (y)-2, (z),  X) \
           ) * (idy))
#define dfx_dz_2b(f, x, y, z, idz) \
         ((+1.5 * _((f), (x), (y), (z),    X) \
           -2.0 * _((f), (x), (y), (z)-1,  X) \
           +0.5 * _((f), (x), (y), (z)-2,  X) \
           ) * (idz))
#define dfy_dx_2b(f, x, y, z, idx) \
         ((+1.5 * _((f), (x),   (y), (z),  Y) \
           -2.0 * _((f), (x)-1, (y), (z),  Y) \
           +0.5 * _((f), (x)-2, (y), (z),  Y) \
           ) * (idx))
#define dfy_dy_2b(f, x, y, z, idy) \
         ((+1.5 * _((f), (x), (y),   (z),  Y) \
           -2.0 * _((f), (x), (y)-1, (z),  Y) \
           +0.5 * _((f), (x), (y)-2, (z),  Y) \
           ) * (idy))
#define dfy_dz_2b(f, x, y, z, idz) \
         ((+1.5 * _((f), (x), (y), (z),    Y) \
           -2.0 * _((f), (x), (y), (z)-1,  Y) \
           +0.5 * _((f), (x), (y), (z)-2,  Y) \
           ) * (idz))
#define dfz_dx_2b(f, x, y, z, idx) \
         ((+1.5 * _((f), (x),   (y), (z),  Z) \
           -2.0 * _((f), (x)-1, (y), (z),  Z) \
           +0.5 * _((f), (x)-2, (y), (z),  Z) \
           ) * (idx))
#define dfz_dy_2b(f, x, y, z, idy) \
         ((+1.5 * _((f), (x), (y),   (z),  Z) \
           -2.0 * _((f), (x), (y)-1, (z),  Z) \
           +0.5 * _((f), (x), (y)-2, (z),  Z) \
           ) * (idy))
#define dfz_dz_2b(f, x, y, z, idz) \
         ((+1.5 * _((f), (x), (y), (z),    Z) \
           -2.0 * _((f), (x), (y), (z)-1,  Z) \
           +0.5 * _((f), (x), (y), (z)-2,  Z) \
           ) * (idz))

// Central partial differences with second order accuracy 
#define dfx_dx_2(f, x, y, z, idx) \
         ((_((f), (x)+1, (y),  (z),   X) - _((f), (x)-1, (y),  (z),   X)) * (idx) * .5)
#define dfx_dy_2(f, x, y, z, idy) \
         ((_((f), (x),  (y)+1, (z),   X) - _((f), (x),  (y)-1, (z),   X)) * (idy) * .5)
#define dfx_dz_2(f, x, y, z, idz) \
         ((_((f), (x),  (y),   (z)+1, X) - _((f), (x),  (y),   (z)-1, X)) * (idz) * .5)
#define dfy_dx_2(f, x, y, z, idx) \
         ((_((f), (x)+1, (y),  (z),   Y) - _((f), (x)-1, (y),  (z),   Y)) * (idx) * .5)
#define dfy_dy_2(f, x, y, z, idy) \
         ((_((f), (x),  (y)+1, (z),   Y) - _((f), (x),  (y)-1, (z),   Y)) * (idy) * .5)
#define dfy_dz_2(f, x, y, z, idz) \
         ((_((f), (x),  (y),   (z)+1, Y) - _((f), (x),  (y),   (z)-1, Y)) * (idz) * .5)
#define dfz_dx_2(f, x, y, z, idx) \
         ((_((f), (x)+1, (y),  (z),   Z) - _((f), (x)-1, (y),  (z),   Z)) * (idx) * .5)
#define dfz_dy_2(f, x, y, z, idy) \
         ((_((f), (x),  (y)+1, (z),   Z) - _((f), (x),  (y)-1, (z),   Z)) * (idy) * .5)
#define dfz_dz_2(f, x, y, z, idz) \
         ((_((f), (x),  (y),   (z)+1, Z) - _((f), (x),  (y),   (z)-1, Z)) * (idz) * .5)

// Central partial differences with fourth order accuracy 
#define dfx_dx_4(f, x, y, z, idx) \
         (((_((f), (x)+1, (y),  (z),   X) - _((f), (x)-1, (y),  (z),   X)) * TWO_THIRDS - \
           (_((f), (x)+2, (y),  (z),   X) - _((f), (x)-2, (y),  (z),   X)) * ONE_TWELFTH) * (idx))
#define dfx_dy_4(f, x, y, z, idy) \
         (((_((f), (x),  (y)+1, (z),   X) - _((f), (x),  (y)-1, (z),   X)) * TWO_THIRDS - \
           (_((f), (x),  (y)+2, (z),   X) - _((f), (x),  (y)-2, (z),   X)) * ONE_TWELFTH) * (idy))
#define dfx_dz_4(f, x, y, z, idz) \
         (((_((f), (x),  (y),  (z)+1, X) - _((f), (x),  (y),  (z)-1, X)) * TWO_THIRDS - \
           (_((f), (x),  (y),  (z)+2, X) - _((f), (x),  (y),  (z)-2, X)) * ONE_TWELFTH) * (idz))
#define dfy_dx_4(f, x, y, z, idx) \
         (((_((f), (x)+1, (y),  (z),   Y) - _((f), (x)-1, (y),  (z),   Y)) * TWO_THIRDS - \
           (_((f), (x)+2, (y),  (z),   Y) - _((f), (x)-2, (y),  (z),   Y)) * ONE_TWELFTH) * (idx))
#define dfy_dy_4(f, x, y, z, idy) \
         (((_((f), (x),  (y)+1, (z),   Y) - _((f), (x),  (y)-1, (z),   Y)) * TWO_THIRDS - \
           (_((f), (x),  (y)+2, (z),   Y) - _((f), (x),  (y)-2, (z),   Y)) * ONE_TWELFTH) * (idy))
#define dfy_dz_4(f, x, y, z, idz) \
         (((_((f), (x),  (y),  (z)+1, Y) - _((f), (x),  (y),  (z)-1, Y)) * TWO_THIRDS - \
           (_((f), (x),  (y),  (z)+2, Y) - _((f), (x),  (y),  (z)-2, Y)) * ONE_TWELFTH) * (idz))
#define dfz_dx_4(f, x, y, z, idx) \
         (((_((f), (x)+1, (y),  (z),   Z) - _((f), (x)-1, (y),  (z),   Z)) * TWO_THIRDS - \
           (_((f), (x)+2, (y),  (z),   Z) - _((f), (x)-2, (y),  (z),   Z)) * ONE_TWELFTH) * (idx))
#define dfz_dy_4(f, x, y, z, idy) \
         (((_((f), (x),  (y)+1, (z),   Z) - _((f), (x),  (y)-1, (z),   Z)) * TWO_THIRDS - \
           (_((f), (x),  (y)+2, (z),   Z) - _((f), (x),  (y)-2, (z),   Z)) * ONE_TWELFTH) * (idy))
#define dfz_dz_4(f, x, y, z, idz) \
         (((_((f), (x),  (y),  (z)+1, Z) - _((f), (x),  (y),  (z)-1, Z)) * TWO_THIRDS - \
           (_((f), (x),  (y),  (z)+2, Z) - _((f), (x),  (y),  (z)-2, Z)) * ONE_TWELFTH) * (idz))

// Approximate the Jacobian
// Add 1.0 along the diagonal to add the identity component 
// of the transform to the displacement field.

// Second order accuracy
#define Jacobian_2(f, x, y, z, idx, idy, idz) \
   det3j(dfx_dx_2((f),(x),(y),(z),(idx)), dfx_dy_2((f),(x),(y),(z),(idy)), dfx_dz_2((f),(x),(y),(z),(idz)),  \
         dfy_dx_2((f),(x),(y),(z),(idx)), dfy_dy_2((f),(x),(y),(z),(idy)), dfy_dz_2((f),(x),(y),(z),(idz)),  \
         dfz_dx_2((f),(x),(y),(z),(idx)), dfz_dy_2((f),(x),(y),(z),(idy)), dfz_dz_2((f),(x),(y),(z),(idz)));

// Fourth order accuracy
#define Jacobian_4(f, x, y, z, idx, idy, idz) \
   det3j(dfx_dx_4((f),(x),(y),(z),(idx)), dfx_dy_4((f),(x),(y),(z),(idy)), dfx_dz_4((f),(x),(y),(z),(idz)),  \
         dfy_dx_4((f),(x),(y),(z),(idx)), dfy_dy_4((f),(x),(y),(z),(idy)), dfy_dz_4((f),(x),(y),(z),(idz)),  \
         dfz_dx_4((f),(x),(y),(z),(idx)), dfz_dy_4((f),(x),(y),(z),(idy)), dfz_dz_4((f),(x),(y),(z),(idz)));

// OBS! Unhygienic macros for short
#define dfx_dx_f(x,y,z) dfx_dx_2f((f), (x), (y), (z), (idx)) 
#define dfx_dy_f(x,y,z) dfx_dy_2f((f), (x), (y), (z), (idy)) 
#define dfx_dz_f(x,y,z) dfx_dz_2f((f), (x), (y), (z), (idz)) 
#define dfy_dx_f(x,y,z) dfy_dx_2f((f), (x), (y), (z), (idx)) 
#define dfy_dy_f(x,y,z) dfy_dy_2f((f), (x), (y), (z), (idy)) 
#define dfy_dz_f(x,y,z) dfy_dz_2f((f), (x), (y), (z), (idz)) 
#define dfz_dx_f(x,y,z) dfz_dx_2f((f), (x), (y), (z), (idx)) 
#define dfz_dy_f(x,y,z) dfz_dy_2f((f), (x), (y), (z), (idy)) 
#define dfz_dz_f(x,y,z) dfz_dz_2f((f), (x), (y), (z), (idz)) 
                                                             
#define dfx_dx_b(x,y,z) dfx_dx_2b((f), (x), (y), (z), (idx)) 
#define dfx_dy_b(x,y,z) dfx_dy_2b((f), (x), (y), (z), (idy)) 
#define dfx_dz_b(x,y,z) dfx_dz_2b((f), (x), (y), (z), (idz)) 
#define dfy_dx_b(x,y,z) dfy_dx_2b((f), (x), (y), (z), (idx)) 
#define dfy_dy_b(x,y,z) dfy_dy_2b((f), (x), (y), (z), (idy)) 
#define dfy_dz_b(x,y,z) dfy_dz_2b((f), (x), (y), (z), (idz)) 
#define dfz_dx_b(x,y,z) dfz_dx_2b((f), (x), (y), (z), (idx)) 
#define dfz_dy_b(x,y,z) dfz_dy_2b((f), (x), (y), (z), (idy)) 
#define dfz_dz_b(x,y,z) dfz_dz_2b((f), (x), (y), (z), (idz)) 
                                                             
#define dfx_dx_c(x,y,z) dfx_dx_2((f), (x), (y), (z), (idx))  
#define dfx_dy_c(x,y,z) dfx_dy_2((f), (x), (y), (z), (idy))  
#define dfx_dz_c(x,y,z) dfx_dz_2((f), (x), (y), (z), (idz))  
#define dfy_dx_c(x,y,z) dfy_dx_2((f), (x), (y), (z), (idx))  
#define dfy_dy_c(x,y,z) dfy_dy_2((f), (x), (y), (z), (idy))  
#define dfy_dz_c(x,y,z) dfy_dz_2((f), (x), (y), (z), (idz))  
#define dfz_dx_c(x,y,z) dfz_dx_2((f), (x), (y), (z), (idx))  
#define dfz_dy_c(x,y,z) dfz_dy_2((f), (x), (y), (z), (idy))  
#define dfz_dz_c(x,y,z) dfz_dz_2((f), (x), (y), (z), (idz))  

/*!
 * \brief Regularise a Jacobian map.
 *
 * Replace with `epsilon` all values of the Jacobian map that
 * fall below `epsilon`.
 */

void regularise(
        const size_t nx,        /*!< Width of the image */
        const size_t ny,        /*!< Length of the image */
        const size_t nz,        /*!< Depth of the image */
        FLOATING J[nz][ny][nx], /*!< Jacobian map */
        FLOATING epsilon        /*!< Minimum value allowed */
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
static inline void jacobian(
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
    #pragma omp parallel for collapse(3) schedule(static)
    for (size_t z = z_start; z < z_stop; ++z) {
        for (size_t y = y_start; y < y_stop; ++y) {
            for (size_t x = x_start; x < x_stop; ++x) {

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

#undef TWO_THIRDS
#undef ONE_TWELFTH

#undef dfx_dx_2f
#undef dfx_dy_2f
#undef dfx_dz_2f
#undef dfy_dx_2f
#undef dfy_dy_2f
#undef dfy_dz_2f
#undef dfz_dx_2f
#undef dfz_dy_2f
#undef dfz_dz_2f

#undef dfx_dx_2b
#undef dfx_dy_2b
#undef dfx_dz_2b
#undef dfy_dx_2b
#undef dfy_dy_2b
#undef dfy_dz_2b
#undef dfz_dx_2b
#undef dfz_dy_2b
#undef dfz_dz_2b

#undef dfx_dx_2
#undef dfx_dy_2
#undef dfx_dz_2
#undef dfy_dx_2
#undef dfy_dy_2
#undef dfy_dz_2
#undef dfz_dx_2
#undef dfz_dy_2
#undef dfz_dz_2

#undef dfx_dx_4
#undef dfx_dy_4
#undef dfx_dz_4
#undef dfy_dx_4
#undef dfy_dy_4
#undef dfy_dz_4
#undef dfz_dx_4
#undef dfz_dy_4
#undef dfz_dz_4

#undef Jacobian_2
#undef Jacobian_4

#undef dfx_dx_f
#undef dfx_dy_f
#undef dfx_dz_f
#undef dfy_dx_f
#undef dfy_dy_f
#undef dfy_dz_f
#undef dfz_dx_f
#undef dfz_dy_f
#undef dfz_dz_f

#undef dfx_dx_b
#undef dfx_dy_b
#undef dfx_dz_b
#undef dfy_dx_b
#undef dfy_dy_b
#undef dfy_dz_b
#undef dfz_dx_b
#undef dfz_dy_b
#undef dfz_dz_b

#undef dfx_dx_c
#undef dfx_dy_c
#undef dfx_dz_c
#undef dfy_dx_c
#undef dfy_dy_c
#undef dfy_dz_c
#undef dfz_dx_c
#undef dfz_dy_c
#undef dfz_dz_c

#endif // __JACOBIAN_H_DEFINED
