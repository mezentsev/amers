set(ADVACTION_SOURCE advaction.c)
set(AMR_SOURCE amr.c amr.h)

include_directories(${P4EST_INCLUDE_DIR} ${P4EST_LIBRARY})
link_directories(${P4EST_LIBRARY})

add_executable(advaction ${ADVACTION_SOURCE})
target_link_libraries (advaction p4est sc)

add_executable(amr ${AMR_SOURCE})
target_link_libraries (amr p4est sc)

# install them
#install (TARGETS advaction RUNTIME DESTINATION bin)