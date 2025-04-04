#ifndef STRUCTURED_FINITE_VOLUME_UTILS_PROJECT_DEFS_H
#define STRUCTURED_FINITE_VOLUME_UTILS_PROJECT_DEFS_H

#include <cstddef>

#include "Kokkos_Core.hpp"


namespace structured_fv {

using UInt = unsigned int;
using Int = int;
using Real = double;
using GlobalDof = long long;

using HostExecutionSpace = Kokkos::DefaultHostExecutionSpace;
using HostMemorySpace = HostExecutionSpace::memory_space;

}

#endif