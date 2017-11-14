#include <stdio.h>

#include "../headers/rvf_io.h"

int main(void)
{
    Image field = read_rvf("../example.rvf");

    print_image_info(field);

    write_rvf("../out.rvf", field);

    delete_image(&field);

    return 0;
}
