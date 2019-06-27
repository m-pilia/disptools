# _disptools: pure C library with a Python wrapper

# Organisation

This folder contains the sources for the C library and its Python wrapper:
+ `headers/`: contains the header files for the C library. The main header is `headers/disptools.h`, which provides the fundamental types and macros. `headers/jacobian.h` implements an inline function to compute the Jacobian of a displacement field.
+ `displacement_field_gradient.c`: implements the `gradient` algorithm.
+ `displacement_field_greedy.c`: implements the `greedy` algorithm.
+ `VolumeMatching3D.c`: file from the PREDICT tool, which implements the `matching` algorithm.
+ `rvf_io.c`: provides input and output for displacement field files in RVF format.
+ `vtk_io.c`: provides basic input and output for scalar and vector images in VTK format.
+ `shape_descriptors.c`: implements some shape descriptors.
+ `disptools.c` and `jacobian.c`: implement dynamic utility functions.

# Useful remarks

+ The C library can be compiled either in single or double precision, controlled through the macro `FLOAT_SIZE`. The NumPy type in the wrapper and in the higher level `disptools` package will be automatically selected accordingly to it.

+ The `headers/disptools.h` header provides types for image objects and macros to access the voxels: the `__(img, x, y, z)` macro for scalar images, and the `_(img, x, y, z, d)` macro for vector images. The names may look meaningless, but they are as short and unintrusive as possible, in order to keep the code compact but readable.

+ Vector images are represented in memory as arrays with indices `[d,z,y,x]` where `d` is the component of the field and `x`,`y`,`z` are the coordinates of the voxel. This differs from the ITK memory layout `[z,y,x,d]`, and the choice is motivated to allow better vectorisation of the code when SIMD instructions are enabled. 
