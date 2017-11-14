#include <stdio.h>
#include <stdlib.h>

#include "../headers/vtk_io.h"
#include "../headers/rvf_io.h"
#include "../headers/jacobian.h"

int main(void)
{
    /** Image field = read_rvf("../example.rvf"); */
    Image field = read_vtk("../example.vtk");

    print_image_info(field);

    Image J = new_image(1, field.nx, field.ny, field.nz, field.dx, field.dy, field.dz);

    jacobian(field, J);

    write_vtk("../jacobian_mine.vtk", J);

    delete_image(&J);
    delete_image(&field);

    return 0;
}
