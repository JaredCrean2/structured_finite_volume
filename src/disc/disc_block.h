#ifndef STRUCTURED_FINITE_VOLUME_DISC_BLOCK_H
#define STRUCTURED_FINITE_VOLUME_DISC_BLOCK_H

#include "utils/project_defs.h"
#include "mesh/structured_mesh.h"

namespace structured_fv {
namespace disc {


// A rectangular block of cells + n layers of ghost cells in
// the cardinal directions.  Currently ghost cells to the
// diagonal are not supported
class StructuredBlock
{
  public:
    using CoordsHostView = mesh::StructuredBlock::CoordsHostView;

    StructuredBlock(const mesh::StructuredMesh& mesh, const mesh::StructuredBlock& block, int nghost);

    UInt getBlockId() const;

    mesh::BlockType getBlockType() const;

    std::array<UInt, 2> getCellDimensions() const;

    std::array<UInt, 2> getVertDimensions() const;

    const std::array<UInt, 4> getNumGhostCellsPerDirection() const;

    Range2D getOwnedVerts() const;

    Range2D getOwnedCells() const;

    Range2D getOwnedAndGhostVerts() const;

    Range2D getOwnedAndGhostCells() const;

    std::pair<UInt, UInt> meshVertToBlockVert(UInt i, UInt j) const;

    std::pair<UInt, UInt> blockVertToMeshVert(UInt i, UInt j) const;

  private:

    std::array<UInt, 4> computeNumGhostsPerDirection(const mesh::StructuredMesh& mesh, int nghost);

    const mesh::StructuredBlock& m_mesh_block;
    std::array<UInt, 4> m_num_ghost_cells_per_direction;
};

inline std::array<Real, 2> computeCellCentroid(StructuredBlock::CoordsHostView vertCoords, int i, int j)
{
  return {(vertCoords(i, j, 0) + vertCoords(i+1, j, 0) + 
           vertCoords(i+1, j+1, 0) + vertCoords(i+1, j, 0))/4,
          (vertCoords(i, j, 1) + vertCoords(i+1, j, 1) + 
           vertCoords(i+1, j+1, 1) + vertCoords(i+1, j, 1))/4};
}


}
}

#endif
