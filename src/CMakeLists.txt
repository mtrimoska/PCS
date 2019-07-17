set(PCS_SRC pcs_exec.c pcs.c pcs_storage.c pcs_pollard_rho.c pcs_elliptic_curve_operations.c pcs_struct_hash.c pcs_struct_hash_UNIX.c pcs_struct_PRTL.c pcs_vect_bin.c)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR})

add_executable(pcs_exec ${PCS_SRC})
target_link_libraries(pcs_exec m)
target_link_libraries(pcs_exec gmp)
