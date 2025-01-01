#ifndef STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_BLOCK_H
#define STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_BLOCK_H

#include "utils/project_defs.h"
#include "utils/range.h"

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
  MeshBlockSpec(UInt num_cells_x, UInt num_cells_y, UInt rotation,
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

enum class BlockType
{
  Regular,
  GhostBC
};

// given an index (i, j) in the standard orientation (ie rotation=0), compute
// the location of the corresponding point under the rotation
std::pair<UInt, UInt> applyRotation(UInt rotation, UInt i, UInt j, UInt dimx,  UInt dimy)
{
  switch (rotation)
  {
    case 0: { return {i,   j}; }
    case 1: { return {j,  dimx -1 - i}; }
    case 2: { return {dimx -1 - i, dimy -1 - j}; }
    case 3: { return {dimy -1 - j, i}; }
    default: { throw std::runtime_error("invalid rotation"); }
  }
}

class StructuredBlock
{
  public:
    StructuredBlock(const MeshBlockSpec& spec, UInt block_id, BlockType block_type) :
      m_block_id(block_id),
      m_block_type(block_type),
      m_owned_cell_range(0, spec.rotation % 2 == 0 ? spec.num_cells_x-1 : spec.num_cells_y-1,
                         0, spec.rotation % 2 == 0 ? spec.num_cells_y-1 : spec.num_cells_x-1),
      m_owned_vert_range(0, spec.rotation % 2 == 0 ? spec.num_cells_x : spec.num_cells_y, 
                         0, spec.rotation % 2 == 0 ? spec.num_cells_y : spec.num_cells_x),
      m_owned_vert_coords("block_coords", spec.rotation % 2 == 0 ? spec.num_cells_x+1 : spec.num_cells_y+1,
                                          spec.rotation % 2 == 0 ? spec.num_cells_y+1 : spec.num_cells_x+1),
      m_offset_into_block{0, 0},
      m_all_block_size{spec.rotation % 2 == 0 ? spec.num_cells_x : spec.num_cells_y,
                       spec.rotation % 2 == 0 ? spec.num_cells_y : spec.num_cells_x}
    { 
      UInt dimx = spec.num_cells_x + 1;
      UInt dimy = spec.num_cells_y + 1;

      for (UInt i : Range(0, dimx-1))
        for (UInt j : Range(0, dimy-1))
        {
          auto [iprime, jprime] = applyRotation(spec.rotation, i, j, dimx, dimy);
          double x = double(i) / (dimx - 1);
          double y = double(j) / (dimy - 1);
          auto [xprime, yprime] = spec.coord_func(x, y);

          m_owned_vert_coords(iprime, jprime, 0) = xprime;
          m_owned_vert_coords(iprime, jprime, 1) = yprime;
        }
    }

    UInt getBlockId() const { return m_block_id;};  // block IDs are global

    BlockType getBlockType() const { return m_block_type; };

    // these are the verts of owned cells
    // slightly misleading, because we dont define ownership of verts
    Range2D getOwnedVerts() const { return m_owned_vert_range; }

    Range2D getOwnedCells() const { return m_owned_cell_range; };

    const Kokkos::View<Real**[2], HostMemorySpace>& getOwnedVertCoords() const { return m_owned_vert_coords; };

    std::array<UInt, 2> getOffsetIntoBlock() { return m_offset_into_block; };  // owned idx + offset = block idx

    // returns size of the entire block (including non-owned cells)
    std::array<UInt, 2> getAllBlockSize() const { return m_all_block_size; };

  private:
    UInt m_block_id;
    BlockType m_block_type;
    Range2D m_owned_cell_range;
    Range2D m_owned_vert_range;
    Kokkos::View<Real**[2], HostMemorySpace> m_owned_vert_coords;
    std::array<UInt, 2> m_offset_into_block;
    std::array<UInt, 2> m_all_block_size;
};


class StructuredBlockInterface
{
  public:
    UInt getLeftBlockId() const { return m_left_block_id; }

    UInt getRightBlockId() const { return m_right_block_id; }

    //TODO: indexer

  private:
    UInt m_left_block_id;
    UInt m_right_block_id;
    std::array<Int, 2> m_transform;  // for each direction in the left block, gives the
                                     // corresponding direction in the right block
                                     // ex. if the +x direction in the left block is
                                     // the negative y direction in right block, then
                                     // m_transform[0] = -1;
};


}
}

#endif

