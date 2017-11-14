#include <stdio.h>
#include <stdlib.h>

#include "../headers/vtk_io.h"
#include "../headers/rvf_io.h"
#include "../headers/jacobian.h"
#include "../headers/displacement_field_gradient.h"

int main(int argc, char *argv[])
{
    (void) argc;

    Image J = read_vtk(argv[1]);
    Image mask_image = read_vtk(argv[2]);
    Mask mask = mask_from_image(mask_image);

    print_image_info(J);

    Image field = new_image(3, J.nx, J.ny, J.nz, J.dx, J.dy, J.dz);

    generate_displacement_gradient(
            field.nx,
            field.ny,
            field.nz,
            field.dx,
            field.dy,
            field.dz,
            (void*) J.data,
            (void*) mask.data,
            .99e-3,
            .2,
            .4,
            2.0,
            0.5,
            0.3,
            0.1,
            1.0,
            false,
            10000,
            (void*) field.data
            );

    write_vtk(argv[3], J);

    delete_image(&J);
    delete_image(&mask_image);
    delete_mask(&mask);
    delete_image(&field);

    return 0;
}
