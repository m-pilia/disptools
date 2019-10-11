set(PYTHON_C_EXTENSION_SRCS
    "_disptools.c"
)

add_library(python_c_extension SHARED ${PYTHON_C_EXTENSION_SRCS})

set_target_properties(
    python_c_extension
    PROPERTIES
        PREFIX ""
        OUTPUT_NAME ${DISPTOOLS_PYTHON_C_MODULE_NAME}
        LINKER_LANGUAGE C
)

if(WIN32)
    set_target_properties(
        python_c_extension
        PROPERTIES
        SUFFIX ".pyd"
    )
endif()

target_include_directories(python_c_extension PUBLIC
    ${PYTHON_INCLUDE_DIRS}
    ${NUMPY_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/src/headers
    )

target_link_libraries(python_c_extension disptools)

if(DISPTOOLS_CUDA_SUPPORT)
    target_include_directories(python_c_extension PUBLIC
        ${PROJECT_SOURCE_DIR}/cuda/headers
        )

    target_link_libraries(python_c_extension disptools_cuda)
endif()

# On Windows, it is required to link to the Python libraries
if(WIN32)
    target_link_libraries(python_c_extension ${PYTHON_LIBRARIES})
endif()

