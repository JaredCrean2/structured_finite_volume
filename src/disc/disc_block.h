#ifndef STRUCTURED_FINITE_VOLUME_DISC_BLOCK_H
#define STRUCTURED_FINITE_VOLUME_DISC_BLOCK_H

#include "utils/project_defs.h"
#include "mesh/structured_mesh.h"
#include "utils/face_iter_per_direction.h"

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

    // owned faces is a bit misleading: it means faces of owned cells.
    // if include_boundary is true, the faces on the boundary (ie. faces shared
    // between an owned cell and a non-owned cell) are included
    FaceRangePerDirection getOwnedFaces(bool include_boundary=true) const;

    // gets all x direction faces for cells in the rectangle [ownedOrGhost x owned]
    FaceRangePerDirection getOwnedAndGhostXFaces(bool include_boundary=true) const;
    
    // gets all y direction faces for cells in the rectangel [owned x ownedOrGhost]
    FaceRangePerDirection getOwnedAndGhostYFaces(bool include_boundary=true) const;

    // gets all x direction faces for cells in the rectangle [ownedOrGhost x ownedAndGhost]
    FaceRangePerDirection getOwnedAndGhostFacesWithCorners(bool include_boundary=true) const;

    std::pair<UInt, UInt> meshVertToBlockVert(UInt i, UInt j) const;

    std::pair<UInt, UInt> blockVertToMeshVert(UInt i, UInt j) const;

  private:

    std::array<UInt, 4> computeNumGhostsPerDirection(const mesh::StructuredMesh& mesh, int nghost);

    const mesh::StructuredBlock& m_mesh_block;
    std::array<UInt, 4> m_num_ghost_cells_per_direction;
};

inline std::array<Real, 2> computeCellCentroid(StructuredBlock::CoordsHostView vert_coords, int i, int j)
{
  return {(vert_coords(i, j, 0) + vert_coords(i+1, j, 0) + 
           vert_coords(i+1, j+1, 0) + vert_coords(i, j+1, 0))/4,
          (vert_coords(i, j, 1) + vert_coords(i+1, j, 1) + 
           vert_coords(i+1, j+1, 1) + vert_coords(i+1, j+1, 1))/4};
}


}
}

#endif
