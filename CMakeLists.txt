cmake_minimum_required(VERSION 3.0)

project(disptools)

option(DISPTOOLS_DEBUG "Debug build" OFF)
option(DISPTOOLS_OPT "Enable platform-specific optimisations" OFF)
option(DISPTOOLS_VERBOSE "Enable verbose output" ON)
option(DISPTOOLS_LOW_ORDER_PD "Lower order for finite differences" OFF)
option(DISPTOOLS_DOUBLE "Use double precision" OFF)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(DISPTOOLS_DEBUG)
    add_definitions(-DDISPTOOLS_DEBUG=1)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g")
else()
    add_definitions(-DDISPTOOLS_DEBUG=0)
endif()

if(DISPTOOLS_OPT)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native")
endif()

if(DISPTOOLS_VERBOSE)
    add_definitions(-DDISPTOOLS_VERBOSE=1)
else()
    add_definitions(-DDISPTOOLS_VERBOSE=0)
endif()

if(DISPTOOLS_LOW_ORDER_PD)
    add_definitions(-DORDER_PD=2)
else()
    add_definitions(-DORDER_PD=4)
endif()

if(DISPTOOLS_DOUBLE)
    add_definitions(-DFLOAT_SIZE=64)
else()
    add_definitions(-DFLOAT_SIZE=32)
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -std=c99 -fPIC")

add_subdirectory(src)

