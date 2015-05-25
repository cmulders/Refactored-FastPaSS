cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

MESSAGE( STATUS "Building in ${BUILDENV} Environment" )

project(FastPass)

add_definitions("-std=c++11")
set(CMAKE_C_FLAGS_DEBUG "-DDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG")


add_subdirectory(vendor)
include_directories(vendor)

add_subdirectory(src)
