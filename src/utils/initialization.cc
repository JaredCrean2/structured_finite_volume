#include "mpi.h"
#include "initialization.h"
#include "Kokkos_Core.hpp"

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

  Kokkos::initialize(argc, argv);
}

void finalize()
{
  int flag = false;
  MPI_Finalized(&flag);  
  if (initializedMPI() && !flag)
  {
    MPI_Finalize();
  }

  Kokkos::finalize();
}
}