#ifndef GENERATE_DISPLACEMENT_CUH_INCLUDED
#define GENERATE_DISPLACEMENT_CUH_INCLUDED

#ifdef __cplusplus
    extern "C" {
#endif

#include "disptools.h"


bool set_device(const int id);


void generate_displacement_greedy_cuda(
        const size_t nx,          /*!< Width of the image  */
        const size_t ny,          /*!< Length of the image */
        const size_t nz,          /*!< Depth of the image  */
        const FLOATING dx,        /*!< x spacing */
        const FLOATING dy,        /*!< y spacing */
        const FLOATING dz,        /*!< z spacing */
        const FLOATING *J,        /*!< Target Jacobian */
        const bool *mask,         /*!< Body mask */
        const FLOATING epsilon,   /*!< Tolerance on the Jacobian per voxel */
        const FLOATING tolerance, /*!< Jacobian tolerance on background */
        FLOATING eta,             /*!< Initial step length for the optimisation */
        const FLOATING eta_max,   /*!< Maximum step length allowed */
        const FLOATING alpha,     /*!< Step length increase coefficient */
        const FLOATING beta,      /*!< Step length decrease coefficient */
        const FLOATING gamma,     /*!< Armijo-Goldstein parameter */
        const FLOATING delta,     /*!< Jacobian regularisation threshold */
        const FLOATING zeta,      /*!< Jacobian regularisation weight */
        const FLOATING theta,     /*!< Termination condition based on improvement */
        const FLOATING iota,      /*!< Termination condition based on eta */
        const bool strict,        /*!< Always improve maximum voxel error */
        const size_t it_max,      /*!< Maximum number of iterations */
        FLOATING *field           /*!< Resulting displacement field */
        );

void generate_displacement_gradient_cuda(
        const size_t nx,          /*!< Width of the image  */
        const size_t ny,          /*!< Length of the image */
        const size_t nz,          /*!< Depth of the image  */
        const FLOATING dx,        /*!< x spacing */
        const FLOATING dy,        /*!< y spacing */
        const FLOATING dz,        /*!< z spacing */
        const FLOATING *J,        /*!< Target Jacobian */
        const bool *mask,         /*!< Body mask */
        const FLOATING epsilon,   /*!< Tolerance on the Jacobian per voxel */
        const FLOATING tolerance, /*!< Jacobian tolerance on background */
        FLOATING eta,             /*!< Step length for the optimisation */
        const FLOATING eta_max,   /*!< Maximum step length allowed */
        const FLOATING alpha,     /*!< Step length increase coefficient */
        const FLOATING beta,      /*!< Step length decrease coefficient */
        const FLOATING gamma,     /*!< Armijo-Goldstein parameter */
        const FLOATING delta,     /*!< Jacobian regularisation threshold */
        const FLOATING zeta,      /*!< Jacobian regularisation weight */
        const FLOATING theta,     /*!< Termination condition based on improvement */
        const FLOATING iota,      /*!< Termination condition based on eta */
        const bool strict,        /*!< Always improve maximum voxel error */
        const size_t it_max,      /*!< Maximum number of iterations */
        FLOATING *field           /*!< Resulting displacement field */
        );

#ifdef __cplusplus
}
#endif

#endif /* ifndef GENERATE_DISPLACEMENT_CUH_INCLUDED */
