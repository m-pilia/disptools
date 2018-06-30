#ifndef _FIELD_H_DEFINED
#define _FIELD_H_DEFINED

#include <errno.h>
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
#elif FLOAT_SIZE == 64
    #define FLOATING double
#else
    #error "Invalid FLOAT_SIZE '" #FLOAT_SIZE "'. Accepted values are 32 and 64."
#endif

// Generic error handler that breaks the execution
#define GENERIC_ERROR_HANDLER(condition) \
    if (condition) \
    { \
        fprintf(stderr, "In %s, line %d, function %s: %s\n", __FILE__, __LINE__, __func__, strerror(errno)); \
        BACKTRACE \
        exit(EXIT_FAILURE); \
    }

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

// XOR swap macro
#define XOR_SWAP(a, b) ((a) ^= (b), (b) ^= (a), (a) ^= (b))

// abs macro
#define abs(x) (x > 0 ? (x) : -(x))

// ternary max macro
#define max(a, b, c) (((a) > (b) ? (a) : (b)) > (c) ? ((a) > (b) ? (a) : (b)) : (c))

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

int get_float_type_size(void);

typedef struct Image {
    size_t nd;
    size_t nx;
    size_t ny;
    size_t nz;
    FLOATING dx;
    FLOATING dy;
    FLOATING dz;
    FLOATING *data;
} Image;

typedef struct Mask {
    size_t nx;
    size_t ny;
    size_t nz;
    bool *data;
} Mask;

Image new_image(
        const size_t nd,
        const size_t nx,
        const size_t ny,
        const size_t nz,
        const size_t dx,
        const size_t dy,
        const size_t dz
        );

Mask new_mask(
        const size_t nx,
        const size_t ny,
        const size_t nz
        );

void delete_image(Image *img);

void delete_mask(Mask *img);

void print_image_info(const Image img);

Mask mask_from_image(const Image img);

#endif // _FIELD_H_DEFINED
