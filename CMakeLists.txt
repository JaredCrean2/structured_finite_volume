cmake_minimum_required(VERSION 3.12.0)

# set project name
project(StructuredFiniteVolume)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

find_package(MPI REQUIRED)


# find BLAS
find_package(BLAS)
if (BLAS_FOUND)
  message("Found BLAS libraries: ${BLAS_LIBRARIES}")
  message("Linker flags: ${BLAS_LINKER_FLAGS}")
else()
  message(FATAL_ERROR, "Unable to locate BLAS")
endif()

if (BLA_STATIC)
  add_library(Blas STATIC IMPORTED GLOBAL)
else()
  add_library(Blas SHARED IMPORTED GLOBAL)
  target_link_options(Blas INTERFACE ${BLAS_LINKER_FLAGS})
endif()

list(GET BLAS_LIBRARIES 0 FIRST_BLAS_LIBRARY)
get_filename_component(BLAS_INCLUDE_DIR ${FIRST_BLAS_LIBRARY} DIRECTORY)
set(BLAS_INCLUDE_DIR "${BLAS_INCLUDE_DIR}/../include")
message("BLAS_INCLUDE_DIR = ${BLAS_INCLUDE_DIR}")
target_include_directories(Blas INTERFACE ${BLAS_INCLUDE_DIR})
message("blas libs = ${BLAS_LIBRARIES}")
foreach(lib ${BLAS_LIBRARIES})
  message("lib = ${lib}")
  set_property(TARGET Blas APPEND PROPERTY
               IMPORTED_LOCATION ${lib})
endforeach()


# find LAPACK
find_package(LAPACK)
if (LAPACK_FOUND)
  message("Found LAPACK libraries: ${LAPACK_LIBRARIES}")
else()
  message(FATAL_ERROR, "Unable to locate LAPACK")
endif()

if (BLA_STATIC)
  add_library(Lapack STATIC IMPORTED GLOBAL)
else()
  add_library(Lapack SHARED IMPORTED GLOBAL)
  target_link_options(Lapack INTERFACE ${LAPACK_LINKER_FLAGS})
endif()

list(GET LAPACK_LIBRARIES 0 FIRST_LAPACK_LIBRARY)
message("first lapack lib = ${FIRST_LAPACK_LIBRARY}")
get_filename_component(LAPACK_INCLUDE_DIR ${FIRST_LAPACK_LIBRARY} DIRECTORY)
set(LAPACK_INCLUDE_DIR "${LAPACK_INCLUDE_DIR}/../include")
target_include_directories(Lapack INTERFACE ${LAPACK_INCLUDE_DIR})

foreach(lib ${LAPACK_LIBRARIES})
  message("lib = ${lib}")
  set_property(TARGET Lapack APPEND PROPERTY
               IMPORTED_LOCATION ${lib})
endforeach()


# find Kokkos
find_package(Kokkos REQUIRED)

# find PETSc
#find_package(PETSc)
#if (PETSC_FOUND)
#  add_library(Petsc SHARED IMPORTED GLOBAL)
#  target_include_directories(Petsc INTERFACE ${PETSC_INCLUDES})
#  set_property(TARGET Petsc PROPERTY IMPORTED_LOCATION ${PETSC_LIBRARIES})
#else()
#  message("Unable to locate PETSc")
#endif()


#find_package(PkgConfig REQUIRED)
#message("about to search for petsc")
#pkg_check_modules(PETSc REQUIRED IMPORTED_TARGET)
##pkg_check_modules(PETSc REQUIRED)
#message("Found PETSc = ${PETSc_FOUND}")
#message("PETSc_LIBRARIES = ${PETSc_LIBRARIES}")
#message("PETSc_INCLUDE_DIRS = ${PETSc_INCLUDE_DIRS}")
#message("Found PETSc libraries: ${PETSc_LIBRARY_DIRS}")
#get_target_property(PETSC_INCLUDE_DIRS PkgConfig::PETSc INCLUDE_DIRECTORIES)
#message("include dirs = ${PETSC_INCLUDE_DIRS}")





set(ALL_LIBS "")

add_subdirectory(src)
add_subdirectory(test)
