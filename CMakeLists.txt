cmake_minimum_required(VERSION 3.0)
project(Parallel_Collision_Search)

option(BUILD_TEST "Build tests" OFF)
if(APPLE)
    SET(CMAKE_C_COMPILER clang)
    SET(CMAKE_C_FLAGS "-fopenmp")
else()
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(GMP REQUIRED)

add_subdirectory(src)
