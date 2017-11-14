#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "headers/rvf_io.h"

/*!
 * \brief Write to file a vector field in rvf format.
 */
void write_rvf(
        const char *filename,               /*!< Filename */
        const Image f
        )
{
    size_t i = 0ul;

    // Allocate a buffer and copy the image data in it
    double *buffer = (double*) malloc(3 * f.nx * f.ny * f.nz * sizeof (double));

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
        perror(strerror(errno));
        exit(errno);
    }

    // Write the header
    fprintf(fp, "%lu %lu %lu\n%f %f %f\n", f.nx, f.ny, f.nz, f.dx, f.dy, f.dz);

    // Write the image data
    fwrite(buffer, sizeof (double), 3 * f.nx * f.ny * f.nz, fp);

    fclose(fp);
    free(buffer);
}

/*!
 * \brief Read from file a vector field in rvf format.
 */
Image read_rvf(
        const char *filename /*!< Filename. */
        )
{
    size_t nx, ny, nz;
    float fdx, fdy, fdz;
    FLOATING dx, dy, dz;
    size_t count;
    size_t i = 0ul;

    FILE *fp = fopen(filename, "rb");
    GENERIC_ERROR_HANDLER(!fp);

    // Parse the header

    count = fscanf(fp, "%lu %lu %lu\n", &nx, &ny, &nz);
    GENERIC_ERROR_HANDLER(3 != count);

    count = fscanf(fp, "%f %f %f\n", &fdx, &fdy, &fdz);
    GENERIC_ERROR_HANDLER(3 != count);

    dx = (FLOATING) fdx;
    dy = (FLOATING) fdy;
    dz = (FLOATING) fdz;

    verbose_printf(DISPTOOLS_DEBUG, "%lu %lu %lu %f %f %f\n", nx, ny, nz, dx, dy, dz);

    // Allocate memory for a buffer and read binary data

    const size_t element_count = 3 * nx * ny * nz;

    double *buffer = (double*) malloc(element_count * sizeof (double));
    GENERIC_ERROR_HANDLER(!buffer);

    fread(buffer, sizeof (double), element_count, fp);

    fclose(fp);

    // Allocate memory for the data structure
    Image img = new_image(3, nx, ny, nz, dx, dy, dz);

    // Copy data from the buffer to the data structure
    for (size_t z = 0; z < nz; ++z) {
        for (size_t y = 0; y < ny; ++y) {
            for (size_t x = 0; x < nx; ++x) {
                _(img, x, y, z, X) = buffer[i++];
                _(img, x, y, z, Y) = buffer[i++];
                _(img, x, y, z, Z) = buffer[i++];
            }
        }
    }

    free(buffer);

    return img;
}
