cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

add_definitions(-std=c++14)

set(CXX_FLAGS "-Wall Wpedantic -Wextra -E")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

project(kalman)

#find_package(PCL 1.2 REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

#include_directories(${PCL_INCLUDE_DIRS})
#link_directories(${PCL_LIBRARY_DIRS})
#add_definitions(${PCL_DEFINITIONS})
#list(REMOVE_ITEM PCL_LIBRARIES "vtkproj4")

set(CMAKE_VERBOSE_MAKEFILE ON)
add_executable (main src/kalman_filter.cpp src/main.cpp src/tracking.cpp)
target_link_libraries (main Eigen3::Eigen)




