#include <stdio.h>

#include "../headers/vtk_io.h"

int main(void)
{
    Image field = read_vtk("../example.vtk");

    print_image_info(field);
    printf("First component of the first voxel: %f\n", _(field, 0, 0, 0, 0));

    write_vtk("../out.vtk", field);

    delete_image(&field);

    return 0;
}
