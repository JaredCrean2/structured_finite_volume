#ifndef STRUCTURED_FINITE_VOLUME_MESH_BLOCK_SPEC_H
#define STRUCTURED_FINITE_VOLUME_MESH_BLOCK_SPEC_H

#include "utils/project_defs.h"

namespace structured_fv {
namespace mesh {


struct MeshBlockSpec
{
  // Note: num_cells_x and num_cells_y are the number of cells in each
  //       direction in a rotation=0 coordinate system.
  //       The rotation argument changes the block-local coordinate
  //       system but does not change the number of cells in each direction.
  //       This allows specifying a multi-block mesh using a rotation=0
  //       global coordinate system
  //       The coord_func should also be defined in the rotation=0 global
  //       coordinate system, and should map from the unit square into
  //       the desired shape
  MeshBlockSpec(UInt num_cells_x=0, UInt num_cells_y=0, UInt rotation=0,
                std::function<std::array<double, 2>(double x, double y)> coord_func = 
                  [](double x, double y) {return std::array<double, 2>{x, y};}) :
    num_cells_x(num_cells_x),
    num_cells_y(num_cells_y),
    rotation(rotation),
    coord_func(coord_func)
  {}

  UInt num_cells_x;
  UInt num_cells_y;

  UInt rotation;  // 0 = no rotation, 1 = rotate 90 degrees counter-clockwise, 2 = 180 degrees, 3 = 270 degrees
                  // this allows testing the transform between block interfaces
  std::function<std::array<double, 2>(double x, double y)> coord_func;
};


struct MeshSpec
{
  MeshSpec(UInt num_blocks_x, UInt num_blocks_y, UInt num_bc_ghost_cells=2):
    blocks("mesh_spec_blocks", num_blocks_x, num_blocks_y),
    num_bc_ghost_cells(num_bc_ghost_cells)
  {}

  Kokkos::View<MeshBlockSpec**, HostMemorySpace> blocks;
  UInt num_bc_ghost_cells;
};

}
}

#endif