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
inline std::pair<double, double> applyRotation(UInt rotation, double x, double y)
{
  switch (rotation)
  {
    case 0: { return {x,   y}; }
    case 1: { return {1 - y,  x}; }
    case 2: { return {1 - x, 1 - y}; }
    case 3: { return {y, 1-x}; }
    default: { throw std::runtime_error("invalid rotation"); }
  }
}

class StructuredBlock
{
  public:
    //TODO: the datatype should be const
    using CoordsHostView = Kokkos::View<Real**[2], HostMemorySpace>;
    StructuredBlock(const MeshBlockSpec& spec, UInt block_id, BlockType block_type) :
      m_block_id(block_id),
      m_block_type(block_type),
      m_owned_cell_range(0,  spec.num_cells_x, 0,spec.num_cells_y),
      m_owned_vert_range(0, spec.num_cells_x+1, 0, spec.num_cells_y+1),
      m_owned_vert_coords("block_coords", spec.num_cells_x+1, spec.num_cells_y+1),
      m_offset_into_block{0, 0},
      m_all_block_size{spec.num_cells_x, spec.num_cells_y}
    {
      UInt dimx = spec.num_cells_x + 1;
      UInt dimy = spec.num_cells_y + 1;

      for (UInt i : Range(0, dimx))
        for (UInt j : Range(0, dimy))
        {
          auto [x, y] = applyRotation(spec.rotation, double(i) / (dimx - 1), double(j) / (dimy - 1));
          auto [xprime, yprime] = spec.coord_func(x, y);
          std::cout << "x, y = " << x << ", " << y << ", xprime, yprime = " << xprime << ", " << yprime << std::endl;

          m_owned_vert_coords(i, j, 0) = xprime;
          m_owned_vert_coords(i, j, 1) = yprime;
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

    StructuredBlock(StructuredBlock&& rhs) = default;

    StructuredBlock& operator=(const StructuredBlock& rhs) = delete;

    StructuredBlock& operator=(StructuredBlock&& rhs) = default;

    UInt getBlockId() const { return m_block_id;};  // block IDs are global

    BlockType getBlockType() const { return m_block_type; };

    // these are the verts of owned cells
    // slightly misleading, because we dont define ownership of verts
    Range2D getOwnedVerts() const { return m_owned_vert_range; }

    Range2D getOwnedCells() const { return m_owned_cell_range; };

    const CoordsHostView& getOwnedVertCoords() const { return m_owned_vert_coords; };

    std::array<UInt, 2> getOffsetIntoBlock() const { return m_offset_into_block; };  // owned idx + offset = block idx

    // returns size of the entire block (including non-owned cells)
    std::array<UInt, 2> getAllBlockSize() const { return m_all_block_size; };

  private:
    UInt m_block_id;
    BlockType m_block_type;
    Range2D m_owned_cell_range;
    Range2D m_owned_vert_range;
    CoordsHostView m_owned_vert_coords;
    std::array<UInt, 2> m_offset_into_block;
    std::array<UInt, 2> m_all_block_size;
};


}
}

#endif

