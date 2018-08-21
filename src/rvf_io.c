#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rvf_io.h"

/*!
 * \brief Write to file a vector field in rvf format.
 */
int write_rvf(
        const char *filename,               /*!< Filename */
        const Image f
        )
{
    size_t i = 0ul;

    // Allocate a buffer and copy the image data in it
    double *buffer = (double*) malloc(3 * f.nx * f.ny * f.nz * sizeof (double));
    if (!buffer) {
        DISPTOOLS_SET_ERROR(true, strerror(errno));
        return -1;
    }

    for (size_t z = 0; z < f.nz; ++z) {
        for (size_t y = 0; y < f.ny; ++y) {
            for (size_t x = 0; x < f.nx; ++x) {
                buffer[i++] = _(f, x, y, z, X);
                buffer[i++] = _(f, x, y, z, Y);
                buffer[i++] = _(f, x, y, z, Z);
            }
        }
    }

    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        DISPTOOLS_SET_ERROR(true, strerror(errno));
        free(buffer);
        return -1;
    }

    // Write the header
    fprintf(fp, "%lu %lu %lu\n%f %f %f\n", f.nx, f.ny, f.nz, f.dx, f.dy, f.dz);

    // Write the image data
    fwrite(buffer, sizeof (double), 3 * f.nx * f.ny * f.nz, fp);

    fclose(fp);
    free(buffer);

    return 0;
}

/*!
 * \brief Read from file a vector field in rvf format.
 */
int read_rvf(
        const char *filename, /*!< Filename. */
        Image *image          /*!< Image. */
        )
{
    size_t nx, ny, nz;
    float fdx, fdy, fdz;
    FLOATING dx, dy, dz;
    size_t count;
    size_t i = 0ul;

    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        DISPTOOLS_SET_ERROR(true, strerror(errno));
        return -1;
    }

    // Parse the header

    count = fscanf(fp, "%lu %lu %lu\n", &nx, &ny, &nz);
    if (3 != count) {
        DISPTOOLS_SET_ERROR(true, "read_rvf: expected size (3 integers)");
        return -1;
    }

    count = fscanf(fp, "%f %f %f\n", &fdx, &fdy, &fdz);
    if (3 != count) {
        DISPTOOLS_SET_ERROR(true, "read_rvf: expected spacing (3 floats)");
        return -1;
    }

    dx = (FLOATING) fdx;
    dy = (FLOATING) fdy;
    dz = (FLOATING) fdz;

    verbose_printf(DISPTOOLS_DEBUG, "%lu %lu %lu %f %f %f\n", nx, ny, nz, dx, dy, dz);

    // Allocate memory for a buffer and read binary data

    const size_t element_count = 3 * nx * ny * nz;

    double *buffer = (double*) malloc(element_count * sizeof (double));
    if (!buffer) {
        DISPTOOLS_SET_ERROR(true, strerror(errno));
        return -1;
    }

    fread(buffer, sizeof (double), element_count, fp);

    fclose(fp);

    // Allocate memory for the data structure
    *image = new_image(3, nx, ny, nz, dx, dy, dz);

    if (disptools_error.error) {
        free(buffer);
        return -1;
    }

    // Copy data from the buffer to the data structure
    Image ptr = *image;
    for (size_t z = 0; z < nz; ++z) {
        for (size_t y = 0; y < ny; ++y) {
            for (size_t x = 0; x < nx; ++x) {
                _(ptr, x, y, z, X) = buffer[i++];
                _(ptr, x, y, z, Y) = buffer[i++];
                _(ptr, x, y, z, Z) = buffer[i++];
            }
        }
    }

    free(buffer);

    return 0;
}
