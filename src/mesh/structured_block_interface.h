#ifndef STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_BLOCK_INTERFACE_H
#define STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_BLOCK_INTERFACE_H

#include "adjacent_block_indexer.h"
#include "utils/range.h"

namespace structured_fv {
namespace mesh {


class StructuredBlockInterface
{
  public:

  StructuredBlockInterface(UInt left_block_id, UInt left_block_constant_index,
                           Range left_block_variable_index,
                           NeighborDirection direction, const std::array<Int, 2>& transform,
                           UInt right_block_id, const std::array<UInt, 2>& right_block_min_cell) :
      m_left_block_id(left_block_id),
      m_left_block_direction(direction),
      m_right_block_id(right_block_id),
      m_indexer()
    {
      if (static_cast<int>(direction) % 2 == 0)
      {
        m_left_block_boundary_cells = Range2D(*(left_block_variable_index.begin()), *(left_block_variable_index.end()),
                                              left_block_constant_index, left_block_constant_index+1);
      } else
      {
        m_left_block_boundary_cells = Range2D(left_block_constant_index, left_block_constant_index+1,
                                              *(left_block_variable_index.begin()), *(left_block_variable_index.end()));
      }

      std::array<UInt, 2> left_block_min_cell = {m_left_block_boundary_cells.getXRange()(0), m_left_block_boundary_cells.getYRange()(0)};
      m_indexer = AdjacentBlockIndexer(transform, left_block_min_cell, direction, right_block_min_cell);
    }
      
    UInt getLeftBlockId() const { return m_left_block_id; }

    UInt getRightBlockId() const { return m_right_block_id; }

    NeighborDirection getNeighborDirection() const { return m_left_block_direction; }

    const Range2D& getLeftBlockBoundaryCells() const { return m_left_block_boundary_cells; }

    const AdjacentBlockIndexer& getAdjacentBlockIndexer() const { return m_indexer; }

  private:
    UInt m_left_block_id;
    NeighborDirection m_left_block_direction;  // gives the side of the block the
                                               // interface is one, from the left block perspective
    Range2D m_left_block_boundary_cells;                                                 

    UInt m_right_block_id;
    AdjacentBlockIndexer m_indexer;

};


}
}

#endif