message("Processing source directory")

### Utils
file(GLOB utils_files utils/*.cc)
add_library(utils "${utils_files}")
target_include_directories(utils PUBLIC
                          "${CMAKE_CURRENT_SOURCE_DIR}"
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
)
target_link_libraries(utils Kokkos::kokkos MPI::MPI_CXX)

### Mesh
file(GLOB mesh_files mesh/*.cc)
add_library(mesh "${mesh_files}")
target_include_directories(mesh PUBLIC
                          "${CMAKE_CURRENT_SOURCE_DIR}"
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
)
target_link_libraries(mesh Kokkos::kokkos MPI::MPI_CXX)


### Disc
file(GLOB disc_files disc/*.cc)
add_library(disc "${disc_files}")
target_include_directories(disc PUBLIC
                          "${CMAKE_CURRENT_SOURCE_DIR}"
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
)
target_link_libraries(disc PUBLIC mesh Kokkos::kokkos MPI::MPI_CXX)

### Advection
file(GLOB advection_files physics/advection/*.cc)
add_library(advection "${advection_files}")
target_include_directories(advection PUBLIC
                          "${CMAKE_CURRENT_SOURCE_DIR}"
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
)
target_link_libraries(advection PUBLIC disc utils)

### Nonlinear solvers
file(GLOB nlsolver_files nonlinear_solvers/*.cc)
add_library(nlsolvers "${nlsolver_files}")
target_include_directories(nlsolvers PUBLIC
                          "${CMAKE_CURRENT_SOURCE_DIR}"
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
)
target_link_libraries(nlsolvers PUBLIC disc utils)





#set(ALL_LIBS ${ALL_LIBS} structured_partition PARENT_SCOPE)
#message("ALL_LIBS = ${ALL_LIBS}")

#add_executable(structured_part main.cc)

#target_include_directories(structured_part PUBLIC
#                          "${PROJECT_BINARY_DIR}"
#                          "${PROJECT_SOURCE_DIR}/include"
#                          )
#target_link_libraries(heat_calc PUBLIC heat_lib)


#set(install_exes structured_part)
#install(TARGETS ${install_exes} DESTINATION bin)
