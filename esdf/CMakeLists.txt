cmake_minimum_required(VERSION 3.16)
project(esdf)

set(CMAKE_CXX_STANDARD 14)

include_directories("${PROJECT_SOURCE_DIR}/include/eigen3")
include_directories("${PROJECT_SOURCE_DIR}/include/ESDF")

aux_source_directory(${PROJECT_SOURCE_DIR}/src DIR_SRCS)

add_executable(esdf main.cpp ${DIR_SRCS})