set(SOURCE main.c)

include_directories(${P4EST_INCLUDE_DIRS})

add_executable(amr ${SOURCE})
target_link_libraries (amr ${P4EST_LIBRARIES})

# install them
install (TARGETS amr RUNTIME DESTINATION bin)