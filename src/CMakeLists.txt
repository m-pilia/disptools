set(DISPTOOLS_SRCS
    "displacement_field_gradient.c"
    "displacement_field_greedy.c"
    "field.c"
    "jacobian.c"
    "shape_descriptors.c"
    "VolumeMatching3D.c"
)

add_library(disptools ${DISPTOOLS_SRCS})
target_include_directories(disptools PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

