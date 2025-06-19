cmake_minimum_required(VERSION 3.18.0)

# set project name
project(StructuredFiniteVolume)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

find_package(MPI REQUIRED)


# find BLAS
find_package(BLAS REQUIRED)

get_target_property(BLAS_IFACE_LIBS BLAS::BLAS INTERFACE_LINK_LIBRARIES)
list(GET BLAS_IFACE_LIBS 0 FIRST_BLAS_LIBRARY)
get_filename_component(BLAS_INCLUDE_DIR ${FIRST_BLAS_LIBRARY} DIRECTORY)
set(BLAS_INCLUDE_DIR "${BLAS_INCLUDE_DIR}/../include")
target_include_directories(BLAS::BLAS INTERFACE "${BLAS_INCLUDE_DIR}")

# find LAPACK
find_package(LAPACK REQUIRED)

# find Kokkos
find_package(Kokkos REQUIRED)

# find PETSc
find_package(PETSc)
if (PETSC_FOUND)
  add_library(Petsc SHARED IMPORTED GLOBAL)
  target_include_directories(Petsc INTERFACE ${PETSC_INCLUDES})
  set_property(TARGET Petsc PROPERTY IMPORTED_LOCATION ${PETSC_LIBRARIES})
else()
  message("Unable to locate PETSc")
endif()


set(ALL_LIBS "")

add_subdirectory(src)
add_subdirectory(test)
