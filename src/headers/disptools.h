#ifndef _FIELD_H_DEFINED
#define _FIELD_H_DEFINED

#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

// Optional features relying on GNU extensions
#if defined(__GNUC__) && defined(__linux__)
#include <execinfo.h>
#include <unistd.h>
#define BACKTRACE \
    { \
        size_t BT_BUF_SIZE = 100; \
        void *buffer[BT_BUF_SIZE]; \
        fprintf(stderr, "backtrace:\n"); \
        backtrace_symbols_fd(buffer, backtrace(buffer, BT_BUF_SIZE), STDERR_FILENO); \
    }
#else
#define BACKTRACE
#endif // defined(__GNUC__) && defined(__linux__)

/*!
 * Vector fields are represented as four-dimensional arrays, whose first index
 * is the component of the vector field (in order to allow automatic
 * vectorisation of the code), and the following indices are the spatial
 * coordinates in the volume (in order [z][y][x], to improve data locality for
 * cache performance).
 */

// Components of the vector field
#define X 0
#define Y 1
#define Z 2

// Floating point type
#ifndef FLOAT_SIZE
    #define FLOAT_SIZE 32
#endif

#if FLOAT_SIZE == 32
    #define FLOATING float
    #define FLOATING_MAX FLT_MAX
    #define FLOATING_MIN FLT_MIN
#elif FLOAT_SIZE == 64
    #define FLOATING double
    #define FLOATING_MAX DBL_MAX
    #define FLOATING_MIN DBL_MIN
#else
    #error "Invalid FLOAT_SIZE '" #FLOAT_SIZE "'. Accepted values are 32 and 64."
#endif

// Debug
#ifndef DISPTOOLS_DEBUG
    #define DISPTOOLS_DEBUG 0
#else
    #if DISPTOOLS_DEBUG
        #undef DISPTOOLS_DEBUG
        #define DISPTOOLS_DEBUG 1
    #else
        #undef DISPTOOLS_DEBUG
        #define DISPTOOLS_DEBUG 0
    #endif
#endif

// Verbose
#ifndef DISPTOOLS_VERBOSE
    #define DISPTOOLS_VERBOSE 0
#else
    #if DISPTOOLS_VERBOSE
        #undef DISPTOOLS_VERBOSE
        #define DISPTOOLS_VERBOSE 1
    #else
        #undef DISPTOOLS_VERBOSE
        #define DISPTOOLS_VERBOSE 0
    #endif
#endif

// Verbose output
#if DISPTOOLS_VERBOSE || DISPTOOLS_DEBUG
    #define verbose_printf(cond, ...) if (cond) {fprintf(stdout, __VA_ARGS__);}
#else
    #define verbose_printf(cond, ...)
#endif // DISPTOOLS_VERBOSE || DISPTOOLS_DEBUG

// Debug hash printer
#if defined(__GNUC__) && defined(__linux__) && DISPTOOLS_DEBUG
#include <openssl/md5.h>
#define hash_print(ptr, size) { \
            unsigned char md[MD5_DIGEST_LENGTH]; \
            MD5((unsigned char *) ptr, size, md); \
            for (int i = 0; i < MD5_DIGEST_LENGTH; i++) { \
                fprintf(stderr, "%02x",  md[i]); \
            } \
        }
#else
#define hash_print(...)
#endif // DISPTOOLS_DEBUG

// Parameter assertion and debug print
#define ASSERT_PARAMETERS \
    { \
        assert(alpha > 0.0 && "alpha must be positive"); \
        assert(beta > 0.0 && "beta must be positive"); \
        assert(gamma > 0.0 && "gamma must be positive"); \
        assert(delta > 0.0 && "delta must be positive"); \
        assert(zeta > 0.0 && "zeta must be positive"); \
        assert(eta > 0.0 && "eta must be positive"); \
        assert(eta_max > 0.0 && "eta_max must be positive"); \
        assert(theta >= 0.0 && "Theta must be positive"); \
        assert(iota >= 0.0 && "Iota must be positive"); \
        assert(tolerance >= 0.0 && "Tolerance must be positive"); \
        assert(epsilon > 0.0 && "Epsilon must be positive"); \
        \
        verbose_printf( \
           DISPTOOLS_DEBUG, \
           "%s\n" \
           "nx:        %zu\n" \
           "ny:        %zu\n" \
           "nz:        %zu\n" \
           "dx:        %f\n" \
           "dy:        %f\n" \
           "dz:        %f\n" \
           "alpha:     %e\n" \
           "beta:      %e\n" \
           "gamma:     %e\n" \
           "delta:     %e\n" \
           "epsilon:   %e\n" \
           "zeta:      %e\n" \
           "eta:       %e\n" \
           "eta_max:   %e\n" \
           "theta:     %e\n" \
           "iota:      %e\n" \
           "tolerance: %e\n" \
           "strict:    %d\n" \
           "it_max:    %zu\n", \
           __func__, \
           nx, ny, nz, \
           dx, dy, dz, \
           alpha, beta, gamma, delta, \
           epsilon, zeta, eta, eta_max, theta, iota, \
           tolerance, \
           strict, \
           it_max); \
    }

// XOR swap macro
#define XOR_SWAP(a, b) ((a) ^= (b), (b) ^= (a), (a) ^= (b))

// abs macro
#define abs(x) (x > 0 ? (x) : -(x))

// ternary max macro
#define max3(a, b, c) (((a) > (b) ? (a) : (b)) > (c) ? ((a) > (b) ? (a) : (b)) : (c))

// Pixel access macros
//
// This makes the code easy to read, but is looks like a lousy way to
// index an array, with a lot of repeated multiplications. However, an
// optimising compiler takes care of repeated expressions (at least on
// gcc, there is no substantial difference between this and grouping
// expressions outside the loops).
//
// Using an incrementing pointer instead of indexing can be faster for
// sequential access, but nonsequential access is required here.

#define _(img, x, y, z, d) ( \
        (img).data[(img).nz * (img).ny * (img).nx * (d) + \
                   (img).ny * (img).nx * (z) + \
                   (img).nx * (y) + \
                   (x)])

#define __(img, x, y, z) ( \
        (img).data[(img).ny * (img).nx * (z) + \
                   (img).nx * (y) + \
                   (x)])

// Determinant of a 3x3 matrix
#define det3(j11, j12, j13, \
            j21, j22, j23, \
            j31, j32, j33) \
    ( (j11) * ((j22) * (j33) - (j23) * (j32)) - \
      (j12) * ((j21) * (j33) - (j23) * (j31)) + \
      (j13) * ((j21) * (j32) - (j22) * (j31)) )

// Determinant of a 3x3 matrix composed with the identity
#define det3j(j11, j12, j13, \
              j21, j22, j23, \
              j31, j32, j33) \
        det3((j11) + 1.0, (j12),       (j13), \
             (j21),       (j22) + 1.0, (j23), \
             (j31),       (j32),       (j33) + 1.0) \


/*!
 * \brief Error state.
 */
typedef struct disptools_error_state {
    bool error;         /*!< `true` if an error was detected. */
    int line;           /*!< Line number where the error happened. */
    char file[512];     /*!< Source file where the error happened. */
    char function[512]; /*!< Function where the error happened. */
    char message[1024]; /*!< Human-readable description of the error. */
    char trace[4096];   /*!< Optional stack trace when the error was detected. */
} disptools_error_state;

/*!
 * \brief Macro to set the error state.
 *
 * \param condition_ Set the error only if the condition holds.
 * \param message_   Human-readable description of the error.
 */
#define DISPTOOLS_SET_ERROR(condition_, message_) \
    if (condition_) \
    { \
        _disptools_set_error(true, __LINE__, __FILE__, __func__, message_, ""); \
        SET_TRACE; \
    }

/*!
 * \brief Set the error state.
 *
 * \param error    Set to `true` if an error was detected.
 * \param line     Line number where the error happened.
 * \param file     Source file where the error happened.
 * \param function Function where the error happened.
 * \param message  Human-readable description of the error.
 * \param trace    Optional stack trace when the error was detected.
 */
void _disptools_set_error(
        const bool error,
        const int line,
        const char file[512],
        const char function[512],
        const char message[1024],
        const char trace[4096]);

/*!
 * \brief Get a pointer to the current error state.
 *
 * \return A pointer to the error state struct.
 */
disptools_error_state* disptools_get_error(void);

/*!
 * \brief Query the current error state.
 *
 * \return `true` if there is an error state, `false` otherwise.
 */
bool disptools_has_error(void);

/*!
 * \brief Clear the error state.
 *
 * Set to a state with no error. 
 */
void disptools_clear_error(void);

#if defined(__GNUC__) && defined(__linux__)
#include <execinfo.h>
#include <unistd.h>
#define SET_TRACE \
    { \
        size_t BT_BUF_SIZE = 100; \
        void *buffer[BT_BUF_SIZE]; \
        int traces_no = backtrace(buffer, BT_BUF_SIZE); \
        char **traces = backtrace_symbols(buffer, traces_no); \
        for (int i = 0; i < traces_no; ++i) { \
            char *dsp_trace = disptools_get_error()->trace; \
            strncat(dsp_trace, traces[i], 4096 - strlen(dsp_trace) - 1); \
            strncat(dsp_trace, "\n",  4096 - strlen(dsp_trace) - 1); \
        } \
        free(traces); \
    }
#else
#define SET_TRACE 
#endif // defined(__GNUC__) && defined(__linux__)

/*!
 * \brief Return the size (in bit) of the floating point type in use.
 */
int get_float_type_size(void);

/*!
 * \brief Data structure representing an image.
 */
typedef struct Image {
    size_t nd;   /*!< Number of dimensions of the voxel values. */
    size_t nx;   /*!< Size in voxels along the `x` axis. */
    size_t ny;   /*!< Size in voxels along the `y` axis. */
    size_t nz;   /*!< Size in voxels along the `z` axis. */
    FLOATING dx; /*!< Spacing along the `x` axis. */
    FLOATING dy; /*!< Spacing along the `y` axis. */
    FLOATING dz; /*!< Spacing along the `z` axis. */
    FLOATING * __restrict data; /*!< Image data. */
} Image;

/*!
 * \brief Data structure representing a binary mask.
 */
typedef struct Mask {
    size_t nx; /*!< Size in voxels along the `x` axis. */
    size_t ny; /*!< Size in voxels along the `y` axis. */
    size_t nz; /*!< Size in voxels along the `z` axis. */
    bool * __restrict data; /*!< Image data. */
} Mask;

/*!
 * \brief Allocate a new image.
 *
 * \note The image must be deallocated with `delete_image`.
 */
Image new_image(
        const size_t nd,   /*!< Number of dimensions of the voxel values. */
        const size_t nx,   /*!< Size in voxels along the `x` axis. */
        const size_t ny,   /*!< Size in voxels along the `y` axis. */
        const size_t nz,   /*!< Size in voxels along the `z` axis. */
        const FLOATING dx, /*!< Spacing along the `x` axis. */
        const FLOATING dy, /*!< Spacing along the `y` axis. */
        const FLOATING dz  /*!< Spacing along the `z` axis. */
        );

/*!
 * \brief Allocate a new mask.
 *
 * \note The mask must be deallocated with `delete_mask`.
 */
Mask new_mask(
        const size_t nx, /*!< Size in voxels along the `x` axis. */
        const size_t ny, /*!< Size in voxels along the `y` axis. */
        const size_t nz  /*!< Size in voxels along the `z` axis. */
        );

/*!
 * \brief Deallocate an image.
 */
void delete_image(Image *img);

/*!
 * \brief Deallocate a mask.
 */
void delete_mask(Mask *img);

/*!
 * \brief Print the image metadata.
 */
void print_image_info(const Image img);

/*!
 * \brief Obtain a binary mask from an image.
 *
 * The value of the mask will be false on voxels with the value of zero,
 * and true elsewhere.
 *
 * \note The mask must be deallocated with `delete_mask`.
 */
Mask mask_from_image(const Image img);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // _FIELD_H_DEFINED
