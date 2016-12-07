set(SOURCE main.c hw32.h)

include_directories(${P4EST_INCLUDE_DIRS})

add_executable(amr ${SOURCE})
target_link_libraries (amr ${P4EST_LIBRARIES})

# find required packages
#find_package (P4EST REQUIRED)
#find_package (SC REQUIRED)

# install them
install (TARGETS amr RUNTIME DESTINATION bin)