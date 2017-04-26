project(amr)

set(AMR_SOURCE amr.c amr.h)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test_prj)
include_directories(${P4EST_INCLUDE_DIR} ${P4EST_LIBRARY})
link_directories(${P4EST_LIBRARY})

add_executable(amr ${AMR_SOURCE})
target_link_libraries (amr p4est sc)