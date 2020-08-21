cmake_minimum_required(VERSION 3.1)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../cmake/")
include(ConfigureAllModules) 

project(clipper)

add_subdirectory(clipper)

__configure_all()