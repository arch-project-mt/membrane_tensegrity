cmake_minimum_required(VERSION 3.10)    # CMake version check
project(membrane_tensegrity)            # Create project named "membrane_tensegrity"

set(CMAKE_CXX_STANDARD 14)              # set the C++ standard to C++14

# https://stackoverflow.com/questions/41361631/optimize-in-cmake-by-default
 if(NOT CMAKE_BUILD_TYPE)
     set(CMAKE_BUILD_TYPE Release)
 endif()

 # set the compile options
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")   
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# find Eigen3 package and you can add more packages here,  
# you can use #include <Eigen/Dense> in your code
# don't forget to link them later in target_link_libraries
find_package(Eigen3 CONFIG REQUIRED)

# add files that need to be compiled here
add_executable(${PROJECT_NAME}  
        data_io.cpp
        data_io.hpp        
        geometry.h
        rmsd.h
        optimize_structure_from_input.cpp
)

# add libraries that need to be linked here
target_link_libraries(${PROJECT_NAME}
        Eigen3::Eigen
)
