cmake_minimum_required(VERSION 3.16)
project(graph_path_finding)

set(CMAKE_CXX_STANDARD 14)

include_directories("include/eigen3")
include_directories("${PROJECT_SOURCE_DIR}/include/AStar")
include_directories("${PROJECT_SOURCE_DIR}/include/LazyThetaStar")
include_directories("${PROJECT_SOURCE_DIR}/include/DStarLite")

#set(src_AStar src/AStar.cpp)
#add_library(AStar_static STATIC ${src_AStar})
#
#set(src_LazyThetaStar src/LazyThetaStar.cpp)
#add_library(LazyThetaStar_static STATIC ${src_LazyThetaStar})
#
#set(src_DStarLite src/DStarLite.cpp)
#add_library(DStarLite_static STATIC ${src_DStarLite})

aux_source_directory(${PROJECT_SOURCE_DIR}/src DIR_SRCS)

add_executable(graph_path_finding main.cpp ${DIR_SRCS})

#target_link_libraries(graph_path_finding PUBLIC AStar_static LazyThetaStar_static DStarLite_static)
