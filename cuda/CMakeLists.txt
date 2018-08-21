set(DISPTOOLS_CUDA_SRCS
    "device.cu"
    "error.cu"
    "disptools.cu"
    "jacobian.cu"
    "displacement_field_gradient.cu"
    "displacement_field_greedy.cu"
)

cuda_add_library(disptools_cuda STATIC ${DISPTOOLS_CUDA_SRCS})
target_include_directories(disptools_cuda PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/headers
    ${PROJECT_SOURCE_DIR}/src/headers
    )
