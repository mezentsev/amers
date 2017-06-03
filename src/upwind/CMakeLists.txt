project(upwind)

set(UPWIND_SOURCE
        upwind.c upwind.h
        solver.h solver.c
        data.h
        util.h util.c
        timestep.c
        flux.c)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/upwind)
include_directories(${P4EST_INCLUDE_DIR} ${P4EST_LIBRARY})
link_directories(${P4EST_LIBRARY})

add_executable(upwind ${UPWIND_SOURCE})
target_link_libraries (upwind p4est sc)