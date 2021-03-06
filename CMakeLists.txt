cmake_minimum_required(VERSION 3.11)

project(disptools LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(DISPTOOLS_DEBUG "Debug build" OFF)
option(DISPTOOLS_OPT "Enable platform-specific optimisations" OFF)
option(DISPTOOLS_VERBOSE "Enable verbose output" ON)
option(DISPTOOLS_LOW_ORDER_PD "Lower order for finite differences" OFF)
option(DISPTOOLS_DOUBLE "Use double precision" OFF)
option(DISPTOOLS_FAST_MATH "Enable unsafe optimisations (non IEEE 754 compliant)" OFF)
option(DISPTOOLS_OMP_SUPPORT "Build with OpenMP support" ON)
option(DISPTOOLS_CUDA_SUPPORT "Build with CUDA support" OFF)
option(DISPTOOLS_CUDA_ERROR_CHECK "Error check for CUDA calls" ON)
option(DISPTOOLS_CUDA_ERROR_CHECK_SYNC "Sync device when checking errors" ON)
option(DISPTOOLS_PYTHON_SUPPORT "Build Python C extension module" OFF)
option(DISPTOOLS_BUILD_TESTS "Build unit tests" OFF)
option(DISPTOOLS_COVERAGE "Measure code coverage" OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

if(DISPTOOLS_COVERAGE)
    if(MSVC)
        message(FATAL_ERROR "Code coverage not implemented for MSVC")
    endif()

    set(DISPTOOLS_BUILD_TESTS ON)
    set(DISPTOOLS_DEBUG ON)

    set(DISPTOOLS_COVERAGE_FLAGS " -coverage -fprofile-arcs -ftest-coverage ")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${DISPTOOLS_COVERAGE_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -coverage -lgcov")
endif()

if(DISPTOOLS_BUILD_TESTS)
    enable_testing()
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(DISPTOOLS_PYTHON_C_MODULE_NAME "_disptools" CACHE STRING "Name of the C extension module file")

# Find OpenMP if required
if(DISPTOOLS_OMP_SUPPORT)
    find_package(OpenMP REQUIRED)
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
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
    if(MSVC)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Ox")
    else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    endif()
endif()

# Enable non-portable optimisations if required
if(DISPTOOLS_OPT AND NOT DISPTOOLS_DEBUG)
    if(NOT MSVC)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
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
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
endif()

add_subdirectory(src)

# Add CUDA component if requested
if(DISPTOOLS_CUDA_SUPPORT)
    enable_language(CUDA)

    if(DISPTOOLS_FAST_MATH AND NOT DISPTOOLS_DEBUG)
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --use_fast_math")
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

# Build test suite
if(DISPTOOLS_BUILD_TESTS)
    add_subdirectory(test)
endif()

# Add Python component if requested
if(DISPTOOLS_PYTHON_SUPPORT)
    set(Python_FIND_VIRTUALENV FIRST)
    find_package(PythonInterp 3.5 REQUIRED)

    # Find NumPy headers
    exec_program(${PYTHON_EXECUTABLE}
        ARGS "-c \"import numpy; print(numpy.get_include())\""
        OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
        RETURN_VALUE NUMPY_NOT_FOUND
        )
    if(NUMPY_NOT_FOUND)
        message(FATAL_ERROR "NumPy headers not found")
    endif()

    # Find Python headers
    exec_program(${PYTHON_EXECUTABLE}
        ARGS "-c \"import sysconfig; print(sysconfig.get_paths()['include'])\""
        OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS
        RETURN_VALUE PYTHON_INCLUDE_DIRS_NOT_FOUND
        )
    if(PYTHON_INCLUDE_DIRS_NOT_FOUND)
        message(FATAL_ERROR "Python headers not found")
    endif()

    # This goes after, since it uses PythonInterp as a hint
    if(WIN32)
        find_package(PythonLibs 3.5 REQUIRED)
    endif()

    add_subdirectory(python_c_extension)
endif()

