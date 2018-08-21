cmake_minimum_required(VERSION 3.10)

project(disptools)

option(DISPTOOLS_DEBUG "Debug build" OFF)
option(DISPTOOLS_OPT "Enable platform-specific optimisations" OFF)
option(DISPTOOLS_VERBOSE "Enable verbose output" ON)
option(DISPTOOLS_LOW_ORDER_PD "Lower order for finite differences" OFF)
option(DISPTOOLS_DOUBLE "Use double precision" OFF)
option(DISPTOOLS_FAST_MATH "Enable unsafe optimisations (non IEEE 754 compliant)" OFF)
option(DISPTOOLS_OMP_SUPPORT "Build with OpenMP support" OFF)
option(DISPTOOLS_CUDA_SUPPORT "Build with CUDA support" OFF)
option(DISPTOOLS_CUDA_ERROR_CHECK "Error check for CUDA calls" ON)
option(DISPTOOLS_CUDA_ERROR_CHECK_SYNC "Sync device when checking errors" ON)
option(DISPTOOLS_PYTHON_SUPPORT "Build Python C extension module" OFF)

set(DISPTOOLS_PYTHON_C_MODULE_NAME "_disptools" CACHE STRING "Name of the C extension module file")

# Find OpenMP if required
if(DISPTOOLS_OMP_SUPPORT)
    find_package(OpenMP REQUIRED)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
else()
    if(NOT MSVC)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-pragmas")
    endif()
endif()

# Enable debug if required
if(DISPTOOLS_DEBUG)
    add_definitions(-DDISPTOOLS_DEBUG=1)
    if(MSVC)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Zi")
    else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g")
    endif()
else()
    add_definitions(-DDISPTOOLS_DEBUG=0)
endif()

# Enable optimisations if required
if(DISPTOOLS_OPT AND NOT DISPTOOLS_DEBUG)
    if(MSVC)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Ox")
    else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native")
    endif()
endif()

# Enable unsafe optimisations if required
if(DISPTOOLS_FAST_MATH AND NOT DISPTOOLS_DEBUG)
    if(MSVC)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /fp:fast")
    else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffast-math")
    endif()
endif()

# Enable verbosity if required
if(DISPTOOLS_VERBOSE)
    add_definitions(-DDISPTOOLS_VERBOSE=1)
else()
    add_definitions(-DDISPTOOLS_VERBOSE=0)
endif()

# Use lower order PD if required
if(DISPTOOLS_LOW_ORDER_PD)
    add_definitions(-DORDER_PD=2)
else()
    add_definitions(-DORDER_PD=4)
endif()

# Use double precision if required
if(DISPTOOLS_DOUBLE)
    add_definitions(-DFLOAT_SIZE=64)
else()
    add_definitions(-DFLOAT_SIZE=32)
endif()

# Set some compiler flags
if(MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -std=c99 -fPIC")
endif()

add_subdirectory(src)

# Add CUDA component if requested
if(DISPTOOLS_CUDA_SUPPORT)
    find_package(CUDA REQUIRED)

    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_30")
    if(NOT MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    endif()

    add_definitions(-DDISPTOOLS_HAS_CUDA=1)

    if(DISPTOOLS_CUDA_ERROR_CHECK)
        add_definitions(-DDISPTOOLS_CUDA_ERROR_CHECK)
    endif()

    if(DISPTOOLS_CUDA_ERROR_CHECK_SYNC)
        add_definitions(-DDISPTOOLS_CUDA_ERROR_CHECK_SYNC)
    endif()

    add_subdirectory(cuda)
endif()

# Add Python component if requested
if(DISPTOOLS_PYTHON_SUPPORT)
    find_package(PythonLibs REQUIRED)

    if(PYTHONLIBS_VERSION_STRING VERSION_LESS "3.5")
        message(FATAL_ERROR "Python 3.5 or higher is required")
    endif()

    find_package(PythonInterp REQUIRED)
    exec_program(${PYTHON_EXECUTABLE}
                 ARGS "-c \"import numpy; print(numpy.get_include())\""
                 OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
                 RETURN_VALUE NUMPY_NOT_FOUND
                )
    if(NUMPY_NOT_FOUND)
        message(FATAL_ERROR "NumPy headers not found")
    endif()

    add_subdirectory(python_c_extension)
endif()

