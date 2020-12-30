#  minimum required cmake version: 3.1.0
cmake_minimum_required(VERSION 3.1.0)

project(RealSenseDNNExample)

# include realsense2
find_library(REALSENSE2_LIB realsense2)
if(REALSENSE2_LIB)
  message(STATUS "Found RealSense2")
else()
  message(FATAL_ERROR "RealSense2 not found")
endif()

# include opencv
find_package( OpenCV REQUIRED )
if(OPENCV_FOUND)
  message(STATUS "Found OpenCV at ${OpenCV_INCLUDE_DIRS}. Available libs: ${OpenCV_LIBS}")
else()
  message(FATAL_ERROR "OpenCV not found")
endif()

# compile target
#include_directories( ${OpenCV_INCLUDE_DIRS} )
add_executable(rs-dnn rs-dnn.cpp cv-helpers.hpp)
set_property(TARGET rs-dnn PROPERTY CXX_STANDARD 11)
target_link_libraries(rs-dnn ${REALSENSE2_LIB} ${OpenCV_LIBS})


# set_target_properties (rs-dnn PROPERTIES
# 	FOLDER "Examples/OpenCV"
# )

# download models
set(PROTOTEXT_FILE "MobileNetSSD_deploy.prototxt")
set(PROTOTEXT_URL "https://raw.githubusercontent.com/chuanqi305/MobileNet-SSD/f5d072ccc7e3dcddaa830e9805da4bf1000b2836/MobileNetSSD_deploy.prototxt")
set(PROTOTEXT_PATH "${CMAKE_CURRENT_BINARY_DIR}/${PROTOTEXT_FILE}")

set(MODEL_FILE "MobileNetSSD_deploy.caffemodel")
set(MODEL_URL "https://drive.google.com/uc?export=download&id=0B3gersZ2cHIxRm5PMWRoTkdHdHc")
set(MODEL_PATH "${CMAKE_CURRENT_BINARY_DIR}/${MODEL_FILE}")

if(NOT EXISTS "${PROTOTEXT_PATH}")
	message("Downloading ${PROTOTEXT_FILE} into ${CMAKE_CURRENT_BINARY_DIR}")
    file(DOWNLOAD "${PROTOTEXT_URL}" "${PROTOTEXT_PATH}")
endif()

if(NOT EXISTS "${MODEL_PATH}")
	message("Downloading ${MODEL_FILE} into ${CMAKE_CURRENT_BINARY_DIR}")
    file(DOWNLOAD "${MODEL_URL}" "${MODEL_PATH}")
endif()

# install(
# 	TARGETS

# 	rs-dnn

# 	RUNTIME DESTINATION
# 	${CMAKE_INSTALL_PREFIX}/bin
# )
