#ifndef STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_BLOCK_INTERFACE_H
#define STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_BLOCK_INTERFACE_H

#include "adjacent_block_indexer.h"
#include "utils/range.h"

namespace structured_fv {
namespace mesh {

class StructuredBlock;

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

    UInt getBlockIdL() const;

    UInt getBlockIdR() const;

    NeighborDirection getNeighborDirectionL() const;

    NeighborDirection getNeighborDirectionR() const;

    const std::array<Int, 2> getTransformL() const;

    const std::array<Int, 2> getTransformR() const;

    const Range2D& getBoundaryVertsL() const;

    const Range2D& getBoundaryVertsR() const;

    const Range2D& getBoundaryCellsL() const;

    const Range2D& getBoundaryCellsR() const;

    const AdjacentBlockIndexer& getAdjBlockVertIndexerL() const;

    const AdjacentBlockIndexer& getAdjBlockVertIndexerR() const;

    const AdjacentBlockIndexer& getAdjBlockCellIndexerL() const;

    const AdjacentBlockIndexer& getAdjBlockCellIndexerR() const;

  private:

    std::pair<Range2D, AdjacentBlockIndexer> createAdjacentBlockIndexer( const Range2D& block_rangeL, NeighborDirection dirL, 
       const Range& rangeL, const std::array<Int, 2>& transformL,
       const Range2D& block_rangeR, NeighborDirection dirR, const Range& rangeR);

    UInt m_left_block_id;
    NeighborDirection m_dirL;  // gives the side of the block the
                               // interface is on, from the left block perspective
    std::array<Int, 2> m_transformL;
    Range2D m_boundary_vertsL;
    Range2D m_boundary_cellsL;
    AdjacentBlockIndexer m_vert_indexerL;
    AdjacentBlockIndexer m_cell_indexerL;

    UInt m_right_block_id;
    NeighborDirection m_dirR;
    std::array<Int, 2> m_transformR;
    Range2D m_boundary_vertsR;
    Range2D m_boundary_cellsR;
    AdjacentBlockIndexer m_vert_indexerR;
    AdjacentBlockIndexer m_cell_indexerR;
};


}
}

#endif