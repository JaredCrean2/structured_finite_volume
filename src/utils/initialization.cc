#include "mpi.h"
#include "initialization.h"
#include "Kokkos_Core.hpp"
#include "petscsys.h"

namespace structured_fv {

bool& initializedMPI()
{
  static bool initailized_mpi = true;
  return initailized_mpi;
}

void initialize(int& argc, char* argv[])
{
  int flag = false;
  MPI_Initialized(&flag);
  if (!flag)
  {
    initializedMPI() = true;
    MPI_Init(&argc, &argv);
  } else
  {
    initializedMPI() = false;
  }

  if (!Kokkos::is_initialized())
    Kokkos::initialize(argc, argv);

  PetscBool is_petsc_initialized;
  PetscInitialized(&is_petsc_initialized);
  if (!is_petsc_initialized)
    PetscInitialize(&argc, &argv, nullptr, nullptr);
}

void finalize()
{
  PetscBool is_petsc_finalized;
  PetscFinalized(&is_petsc_finalized);
  if (!is_petsc_finalized)
    PetscFinalize();

  if (!Kokkos::is_finalized())
    Kokkos::finalize();

  int flag = false;
  MPI_Finalized(&flag);  
  if (initializedMPI() && !flag)
  {
    MPI_Finalize();
  }
}
}