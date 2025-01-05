#ifndef STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_BLOCK_H
#define STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_BLOCK_H

#include "utils/project_defs.h"
#include "utils/range.h"
#include <iostream>
#include "block_spec.h"

namespace structured_fv {
namespace mesh {

enum class BlockType
{
  Regular,
  GhostBC
};

// given an index (i, j) in the standard orientation (ie rotation=0), compute
// the location of the corresponding point under the rotation
inline std::pair<UInt, UInt> applyRotation(UInt rotation, UInt i, UInt j, UInt dimx,  UInt dimy)
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
      m_owned_cell_range(0, spec.rotation % 2 == 0 ? spec.num_cells_x : spec.num_cells_y,
                         0, spec.rotation % 2 == 0 ? spec.num_cells_y : spec.num_cells_x),
      m_owned_vert_range(0, spec.rotation % 2 == 0 ? spec.num_cells_x+1 : spec.num_cells_y+1, 
                         0, spec.rotation % 2 == 0 ? spec.num_cells_y+1 : spec.num_cells_x+1),
      m_owned_vert_coords("block_coords", spec.rotation % 2 == 0 ? spec.num_cells_x+1 : spec.num_cells_y+1,
                                          spec.rotation % 2 == 0 ? spec.num_cells_y+1 : spec.num_cells_x+1),
      m_offset_into_block{0, 0},
      m_all_block_size{spec.rotation % 2 == 0 ? spec.num_cells_x : spec.num_cells_y,
                       spec.rotation % 2 == 0 ? spec.num_cells_y : spec.num_cells_x}
    { 
      UInt dimx = spec.num_cells_x + 1;
      UInt dimy = spec.num_cells_y + 1;

      for (UInt i : Range(0, dimx))
        for (UInt j : Range(0, dimy))
        {
          auto [iprime, jprime] = applyRotation(spec.rotation, i, j, dimx, dimy);
          double x = double(i) / (dimx - 1);
          double y = double(j) / (dimy - 1);
          auto [xprime, yprime] = spec.coord_func(x, y);

          m_owned_vert_coords(iprime, jprime, 0) = xprime;
          m_owned_vert_coords(iprime, jprime, 1) = yprime;
        }
    }

    StructuredBlock(Kokkos::View<Real**[2], HostMemorySpace> ghost_coords, UInt block_id) :
      m_block_id(block_id),
      m_block_type(BlockType::GhostBC),
      m_owned_cell_range(0, ghost_coords.extent(0)-1, 0, ghost_coords.extent(1)-1),
      m_owned_vert_range(0, ghost_coords.extent(0),    0, ghost_coords.extent(1)),
      m_owned_vert_coords(ghost_coords),
      m_offset_into_block{0, 0},
      m_all_block_size{static_cast<UInt>(ghost_coords.extent(0)-1), static_cast<UInt>(ghost_coords.extent(1)-1)}
    {}

    StructuredBlock(const StructuredBlock& rhs) = delete;

    StructuredBlock& operator=(const StructuredBlock& rhs) = delete;

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


}
}

#endif

