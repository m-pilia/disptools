project(disptools_test LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(GoogleTest)
include(GoogleTest)

if(DISPTOOLS_OMP_SUPPORT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(DISPTOOLS_COVERAGE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 ${DISPTOOLS_COVERAGE_FLAGS}")
endif()

# Headers

include_directories(
    ${googletest_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/src/headers
    )

# Source files

set(TEST_SOURCE_FILES
    "c/test_jacobian_macros.cpp"
    )

add_executable(disptools_test ${TEST_SOURCE_FILES})

# Link libraries

target_link_libraries(disptools_test
    disptools
    gtest
    gtest_main
    pthread
    )

# Test discovery

gtest_discover_tests(disptools_test)
