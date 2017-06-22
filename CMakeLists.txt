cmake_minimum_required(VERSION 2.6)

#set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w -std=c99 -O0 -g -v")
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_MODULES "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_MODULES})

find_package(MPI REQUIRED)
find_package(P4EST REQUIRED)

# compile
add_subdirectory (src)