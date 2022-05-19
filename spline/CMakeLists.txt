cmake_minimum_required(VERSION 3.16)
project(spline)

set(CMAKE_CXX_STANDARD 14)

include_directories("include/eigen3")
include_directories("${PROJECT_SOURCE_DIR}/include/BezierSpline")
include_directories("${PROJECT_SOURCE_DIR}/include/BSpline")

aux_source_directory(${PROJECT_SOURCE_DIR}/src DIR_SRCS)

add_executable(spline main.cpp ${DIR_SRCS})