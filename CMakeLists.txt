cmake_minimum_required(VERSION 3.6)
project(amr)

#set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror -std=c99 -O0 -v")
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_MODULES "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_MODULES})

find_package(P4EST REQUIRED)

#add_executable(amr ${SOURCE_FILES})

# compile
add_subdirectory (src)