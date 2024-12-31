#ifndef STRUCTURED_FINITE_VOLUME_UTILS_INITIALIZATION_H
#define STRUCTURED_FINITE_VOLUME_UTILS_INITIALIZATION_H

#include "mpi.h"

namespace structured_fv {

bool& initializedMPI();

void initialize(int& argc, char* argv[]);

void finalize();

}

#endif