cmake_minimum_required(VERSION 3.10)
project(membrane_tensegrity)

set(CMAKE_CXX_STANDARD 14)

# https://stackoverflow.com/questions/41361631/optimize-in-cmake-by-default
# if(NOT CMAKE_BUILD_TYPE)
#     set(CMAKE_BUILD_TYPE Release)
# endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(Eigen3 CONFIG REQUIRED)

add_executable(${PROJECT_NAME}
        # data_io.cpp
        # data_io.hpp
        # geometry.h
        # rmsd.h
        # rmsd_io.h
        # rmsd_struct.h
#        test_optimization.cpp
        optimize_structure_from_input.cpp
)

target_link_libraries(${PROJECT_NAME}
        Eigen3::Eigen
)



