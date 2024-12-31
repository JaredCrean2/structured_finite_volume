#ifndef STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_BLOCK_H
#define STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_BLOCK_H

#include "utils/project_defs.h"
#include "utils/range.h"

namespace structured_fv {
namespace mesh {


struct MeshBlockSpec
{
  MeshBlockSpec(UInt num_cells_x, UInt num_cells_y, UInt rotation,
                double xmin, double xmax, double ymin, double ymax) :
    num_cells_x(num_cells_x),
    num_cells_y(num_cells_y),
    rotation(rotation),
    coord_func([=](UInt i, UInt j, UInt num_verts_x, UInt num_verts_y)
    { return std::array<double, 2>{i * (xmax - xmin)/(num_verts_x - 1),
                                   j * (ymax - ymin)/(num_verts_y - 1)};
    })
  {}

  UInt num_cells_x;
  UInt num_cells_y;

  UInt rotation;  // 0 = no rotation, 1 = rotate 90 degrees counter-clockwise, 2 = 180 degrees, 3 = 270 degrees
                  // this allows testing the transform between block interfaces
  std::function<std::array<double, 2>(UInt i, UInt j, UInt numVertsX, UInt numVertsY)> coord_func;
};

enum class BlockType
{
  Regular,
  GhostBC
};

class StructuredBlock
{
  public:
    StructuredBlock(const MeshBlockSpec& spec, UInt block_id, BlockType block_type) :
      m_block_id(block_id),
      m_block_type(block_type),
      m_owned_cell_range(0, spec.num_cells_x-1, 0, spec.num_cells_y-1),
      m_owned_vert_range(0, spec.num_cells_x, 0, spec.num_cells_y),
      m_owned_vert_coords("block_coords", spec.num_cells_x+1, spec.num_cells_y+1),
      m_offset_into_block{0, 0},
      m_all_block_size{spec.num_cells_x, spec.num_cells_y}
    {
      for (UInt i : m_owned_vert_range.getXRange())
        for (UInt j : m_owned_vert_range.getYRange())
        {
          std::array<double, 2> coords = spec.coord_func(i, j, m_owned_vert_range.getXRange().size(), m_owned_vert_range.getYRange().size());
          m_owned_vert_coords(i, j, 0) = coords[0];
          m_owned_vert_coords(i, j, 1) = coords[1];
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

