#ifndef STRUCTURED_FINITE_VOLUME_MESH_ADJACENT_BLOCK_INDEXER_H
#define STRUCTURED_FINITE_VOLUME_MESH_ADJACENT_BLOCK_INDEXER_H

#include "utils/project_defs.h"
#include "neighbor_direction.h"

namespace structured_fv {
namespace mesh {


class AdjacentBlockIndexer
{
  public:
    // transform is the transformation of the coordinate system from the
    // left block to the right block, according to the CGNS definition:
    // transform[i] gives the image of the i'th axis of the left block in
    // the right block.  For example, if the +x direction of the left block
    // is the -y direction of the right block, then transform[0] = -2
    // Note: this is a little weird because CGNS uses 1 based indexing, but
    // std::array uses zero based.
    // left_block_min_cell is the minimum cell on the left block that is adjacent
    // to a cell on the right block
    // right_block_min_cell is the index of the adjacent cell in the right block
    AdjacentBlockIndexer(const std::array<Int, 2>& transform,
                         const std::array<UInt, 2>& left_block_min_cell,
                         NeighborDirection direction,
                         const std::array<UInt, 2>& right_block_min_cell) :
      m_transform_matrix("transform_matrix"),
      m_left_block_min_cell(left_block_min_cell),
      m_direction(direction),
      m_right_block_min_cell(right_block_min_cell)
    {
      assert(transform[0] != 0 && transform[1] != 0);

      m_transform_matrix(0, 0) = sgn(transform[0])*del(transform[0], 1);
      m_transform_matrix(0, 1) = sgn(transform[1])*del(transform[1], 1);
      m_transform_matrix(1, 0) = sgn(transform[0])*del(transform[0], 2);
      m_transform_matrix(1, 1) = sgn(transform[1])*del(transform[1], 2);
    }

    AdjacentBlockIndexer() :
      AdjacentBlockIndexer({1, 2}, {0, 0}, NeighborDirection::East, {0, 0})
    {}

    // i and j are used to index into the adjacent block as though it
    // were a continuation of the current block.  For example, if
    // the current block is connected to the adjacent block to the East
    // then (nx, 0) will give the indices of the cell in the adjacent
    // block that is next to the bottom right cell in the current block.
    // (where nx is the number of cells in the x direction of the current block)
    // If the blocks are joined to the West, then (-1, 0) gives the indices of
    // the cell in the adjacent block next to the lower left cell of the current
    // block
    std::array<UInt, 2> operator()(Int i, Int j) const
    {
      std::array<Int, 4> ioffset{0,  -1, 0, 1};
      std::array<Int, 4> joffset{-1,  0, 1, 0};

      i += ioffset[static_cast<int>(m_direction)];
      j += joffset[static_cast<int>(m_direction)];

      Int iprime = m_transform_matrix(0, 0) * (i - m_left_block_min_cell[0]) +
                   m_transform_matrix(0, 1) * (j - m_left_block_min_cell[1]);
      Int jprime = m_transform_matrix(1, 0) * (i - m_left_block_min_cell[0]) +
                   m_transform_matrix(1, 1) * (j - m_left_block_min_cell[1]);

      return {iprime + m_right_block_min_cell[0], jprime + m_right_block_min_cell[1]};
    }

    std::array<UInt, 2> operator()(const std::array<UInt, 2>& ij) const
    {
      return this->operator()(ij[0], ij[1]);
    }

    std::array<UInt, 2> operator()(const std::array<Int, 2>& ij) const
    {
      return this->operator()(ij[0], ij[1]);
    }

    const std::array<UInt, 2>& getRightBlockMinEntity() const { return m_right_block_min_cell; }

  private:
    // TODO: put this on stack?
    Kokkos::View<Int[2][2], HostMemorySpace> m_transform_matrix; 
    std::array<UInt, 2> m_left_block_min_cell;
    NeighborDirection m_direction;
    std::array<UInt, 2> m_right_block_min_cell;
};

}
}
#endif