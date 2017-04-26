set(NS_SOURCE
        ns.c ns.h
        solver.h solver.c
        data.h
        util.h util.c)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/ns)
include_directories(${P4EST_INCLUDE_DIR} ${P4EST_LIBRARY})
link_directories(${P4EST_LIBRARY})

add_executable(ns ${NS_SOURCE})
target_link_libraries (ns p4est sc)