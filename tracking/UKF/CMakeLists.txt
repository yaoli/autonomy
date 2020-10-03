cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

add_definitions(-std=c++14)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

project(kalman)

find_package(PCL 1.2 REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)


add_executable (main src/main.cpp src/ukf.cpp)
target_link_libraries (main Eigen3::Eigen)




