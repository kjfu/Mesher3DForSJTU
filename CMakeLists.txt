# Set  the minimum  required version  of cmake  for a  project.
cmake_minimum_required(VERSION 3.10)

project(mesher3d)
# SET(CMAKE_C_COMPILER clang-10)
# add_compile_options(-std=c++14 -fpermissive)
set(CMAKE_CXX_STANDARD 14)

aux_source_directory(./src DIR_SRCS)


#set(KDTREE_SRC ./src/gridlab/external/kdtree/kdtree.c)

# Add an executable to the project using the specified source files.
add_executable(mesher3d  ${DIR_SRCS})

#Set properties on a target. 
set_target_properties(mesher3d PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY)
