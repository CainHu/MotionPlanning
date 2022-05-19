cmake_minimum_required(VERSION 3.16)
project(hybrid_a_star)

set(CMAKE_CXX_STANDARD 14)

include_directories("include/eigen3")
include_directories("${PROJECT_SOURCE_DIR}/include/HybridAStar")

aux_source_directory(${PROJECT_SOURCE_DIR}/src DIR_SRCS)

add_executable(hybrid_a_star main.cpp ${DIR_SRCS})