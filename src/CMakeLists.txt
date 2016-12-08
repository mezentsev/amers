set(SOURCE main.c hw32.h)

include_directories(${P4EST_INCLUDE_DIR} ${P4EST_LIBRARY})
link_directories(${P4EST_LIBRARY})

add_executable(amr main.c)
target_link_libraries (amr p4est sc)

# install them
install (TARGETS amr RUNTIME DESTINATION bin)