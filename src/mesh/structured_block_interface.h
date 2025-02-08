#ifndef STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_BLOCK_INTERFACE_H
#define STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_BLOCK_INTERFACE_H

#include "adjacent_block_indexer.h"
#include "utils/range.h"

namespace structured_fv {
namespace mesh {


class StructuredBlockInterface
{
  public:

  // rangeL and ranger give corresponding the ranges of indices along the boundary
  // this allows supporting T junctions
  StructuredBlockInterface(const StructuredBlock& blockL,
                          NeighborDirection dirL,
                          const Range& rangeL,
                          const std::array<Int, 2>& transformL,
                          const StructuredBlock& blockR,
                          const Range& rangeR);

    UInt getBlockIdL() const { return m_left_block_id; }

    UInt getBlockIdR() const { return m_right_block_id; }

    NeighborDirection getNeighborDirectionL() const { return m_dirL; }

    NeighborDirection getNeighborDirectionR() const { return m_dirR; }

    const Range2D& getBoundaryCellsL() const { return m_boundary_cellsL; }

    const Range2D& getBoundaryCellsR() const { return m_boundary_cellsR; };

    const AdjacentBlockIndexer& getAdjacentBlockIndexerL() const { return m_indexerL; }

    const AdjacentBlockIndexer& getAdjacentBlockIndexerR() const { return m_indexerR; };

  private:

    std::pair<Range2D, AdjacentBlockIndexer> createAdjacentBlockIndexer( const StructuredBlock& blockL, NeighborDirection dirL, 
       const Range& rangeL, const std::array<Int, 2>& transformL,
       const StructuredBlock& blockR, NeighborDirection dirR, const Range& rangeR);

    UInt m_left_block_id;
    NeighborDirection m_dirL;  // gives the side of the block the
                               // interface is on, from the left block perspective
    std::array<Int, 2> m_transformL;
    Range2D m_boundary_cellsL;
    AdjacentBlockIndexer m_indexerL;

    UInt m_right_block_id;
    NeighborDirection m_dirR;
    std::array<Int, 2> m_transformR;
    Range2D m_boundary_cellsR;
    AdjacentBlockIndexer m_indexerR;
};


}
}

#endif