message("Processing source directory")

### Utils
file(GLOB utils_files utils/*.cc)
add_library(utils "${utils_files}")
target_include_directories(utils PUBLIC
                          "${CMAKE_CURRENT_SOURCE_DIR}"
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
)
target_link_libraries(utils PUBLIC Kokkos::kokkos MPI::MPI_CXX Petsc)

### Mesh
file(GLOB mesh_files mesh/*.cc)
add_library(mesh "${mesh_files}")
target_include_directories(mesh PUBLIC
                          "${CMAKE_CURRENT_SOURCE_DIR}"
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
)
target_link_libraries(mesh Kokkos::kokkos MPI::MPI_CXX utils)


### Disc
file(GLOB disc_files disc/*.cc)
add_library(disc "${disc_files}")
target_include_directories(disc PUBLIC
                          "${CMAKE_CURRENT_SOURCE_DIR}"
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
)
target_link_libraries(disc PUBLIC mesh Kokkos::kokkos MPI::MPI_CXX)

### Physics model
file(GLOB physics_files physics/*.cc)
add_library(physics "${physics_files}")
target_include_directories(physics PUBLIC
                          "${CMAKE_CURRENT_SOURCE_DIR}"
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
)
target_link_libraries(physics PUBLIC disc utils)

### Physics common
file(GLOB physics_common_files physics/common/*.cc)
add_library(physics_common "${physics_common_files}")
target_include_directories(physics_common PUBLIC
                          "${CMAKE_CURRENT_SOURCE_DIR}"
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
)
target_link_libraries(physics_common PUBLIC utils)

### Advection
file(GLOB advection_files physics/advection/*.cc)
add_library(advection "${advection_files}")
target_include_directories(advection PUBLIC
                          "${CMAKE_CURRENT_SOURCE_DIR}"
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
)
target_link_libraries(advection PUBLIC physics physics_common disc utils linear_system)


### Euler
file(GLOB euler_files physics/euler/*.cc)
add_library(euler "${euler_files}")
target_include_directories(euler PUBLIC
                          "${CMAKE_CURRENT_SOURCE_DIR}"
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
)
target_link_libraries(euler PUBLIC physics physics_common disc utils)


### Linear System
file(GLOB linear_system_files linear_system/*.cc)
add_library(linear_system "${linear_system_files}")
target_include_directories(linear_system PUBLIC
                          "${CMAKE_CURRENT_SOURCE_DIR}"
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
)
target_link_libraries(linear_system PUBLIC disc utils BLAS::BLAS LAPACK::LAPACK Petsc)


### Nonlinear solvers
file(GLOB nlsolver_files nonlinear_solvers/*.cc)
add_library(nlsolvers "${nlsolver_files}")
target_include_directories(nlsolvers PUBLIC
                          "${CMAKE_CURRENT_SOURCE_DIR}"
                          "${PROJECT_BINARY_DIR}/include"  # needed for configured header
)
target_link_libraries(nlsolvers PUBLIC physics disc utils)


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
