# License: Apache 2.0. See LICENSE file in root directory.
# Copyright(c) 2019 Intel Corporation. All Rights Reserved.
#  minimum required cmake version: 3.1.0
cmake_minimum_required(VERSION 3.1.0)

project(RealsenseExamplesHelloRealSense)
find_library(REALSENSE2_LIB realsense2)

add_executable(rs-hello-realsense rs-hello-realsense.cpp)
set_property(TARGET rs-hello-realsense PROPERTY CXX_STANDARD 11)
target_link_libraries(rs-hello-realsense "${REALSENSE2_LIB}")
include_directories(rs-hello-realsense, "/usr/include/librealsense2")
