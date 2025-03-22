#ifndef STRUCTURED_FINITE_VOLUME_DISC_INTERFACE_H
#define STRUCTURED_FINITE_VOLUME_DISC_INTERFACE_H

#include "utils/project_defs.h"
#include "utils/neighbor_direction.h"
#include "mesh/structured_block_interface.h"
#include "disc_block.h"

namespace structured_fv {
namespace disc {


class StructuredBlockInterface
{
  public:
    StructuredBlockInterface(const StructuredBlock& blockL, const StructuredBlock& blockR,
                             const mesh::StructuredBlockInterface& mesh_iface);

    UInt getBlockIdL() const;

    UInt getBlockIdR() const;

    const Range2D& getOwnedBoundaryVertsL() const;

    const Range2D& getOwnedBoundaryVertsR() const;

    const Range2D& getOwnedBoundaryCellsL() const;

    const Range2D& getOwnedBoundaryCellsR() const;

    NeighborDirection getNeighborDirectionL() const;

    NeighborDirection getNeighborDirectionR() const;

    const mesh::AdjacentBlockIndexer& getAdjBlockVertIndexerL() const;

    const mesh::AdjacentBlockIndexer& getAdjBlockVertIndexerR() const;

    const mesh::AdjacentBlockIndexer& getAdjBlockCellIndexerL() const;

    const mesh::AdjacentBlockIndexer& getAdjBlockCellIndexerR() const;

  private:
    const StructuredBlock& m_blockL;
    const StructuredBlock& m_blockR;
    const mesh::StructuredBlockInterface m_mesh_iface;

    Range2D m_owned_boundary_vertsL;
    Range2D m_owned_boundary_vertsR;
    Range2D m_owned_boundary_cellsL;
    Range2D m_owned_boundary_cellsR;
    mesh::AdjacentBlockIndexer m_cell_indexerL;
    mesh::AdjacentBlockIndexer m_cell_indexerR;
    mesh::AdjacentBlockIndexer m_vert_indexerL;
    mesh::AdjacentBlockIndexer m_vert_indexerR;    
};

}
}
#endif