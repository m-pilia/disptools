#ifndef JACOBIAN_MACROS_H_INCLUDED
#define JACOBIAN_MACROS_H_INCLUDED

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

#endif /* ifndef JACOBIAN_MACROS_H_INCLUDED */
