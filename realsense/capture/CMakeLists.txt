# License: Apache 2.0. See LICENSE file in root directory.
# Copyright(c) 2019 Intel Corporation. All Rights Reserved.
#  minimum required cmake version: 3.1.0
cmake_minimum_required(VERSION 3.1.0)

project(RealsenseExamplesCapture)
find_library(REALSENSE2_LIB realsense2)
find_library(GLFW3_LIB glfw)
find_package(OpenGL REQUIRED)
add_executable(rs-capture rs-capture.cpp ../example.hpp)
set_property(TARGET rs-capture PROPERTY CXX_STANDARD 11)
target_link_libraries(rs-capture ${REALSENSE2_LIB})
target_link_libraries(rs-capture ${GLFW3_LIB})
target_link_libraries(rs-capture ${OPENGL_LIBRARIES})
include_directories(rs-capture /home/lyao/Repos/librealsense/third-party/tclap/include)
include_directories(rs-capture /home/lyao/Repos/autonomy/realsense)
