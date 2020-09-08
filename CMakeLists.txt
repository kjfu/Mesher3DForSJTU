# Set  the minimum  required version  of cmake  for a  project.
cmake_minimum_required(VERSION 2.6)

SET(CMAKE_C_COMPILER g++)
add_compile_options(-std=c++14)

aux_source_directory(./src DIR_SRCS)

# Add an executable to the project using the specified source files.
add_executable(mesher3d ${DIR_SRCS})

#Set properties on a target. 
set_target_properties(mesher3d PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY)
