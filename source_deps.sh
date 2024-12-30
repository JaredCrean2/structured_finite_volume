#!/bin/bash

# source all the scripts to be able to find the dependencies

export OPENBLAS_NUM_THREADS=1

#source ~/build/gcc-7.3.0_install/use_gcc.sh
source ~/build/llvm/use_clang.sh
source ~/build/mpich-3.4.1_clang_install/use_mpich.sh
source ~/build/kokkos_install/use_kokkos.sh
#source ~/build/core_install/use_core.sh
#source ~/build/boost_1_76_0_install/use_boost.sh
source ~/build/googletest/install/use_gtest.sh
#source ~/build/petsc_3.17.0/debug/use_petsc.sh
#source ~/build/petsc_3.17.0/opt/use_petsc.sh
source ~/build/OpenBLAS_install/use_openblas.sh
source ~/build/paraview/ParaView-5.10.1-MPI-Linux-Python3.9-x86_64/use_paraview.sh
