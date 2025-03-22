#ifndef STRUCTURED_FINITE_VOLUME_DISC_DISCRETIZATION_H
#define STRUCTURED_FINITE_VOLUME_DISC_DISCRETIZATION_H

//#include "mesh/structured_block.h"
//#include "mesh/structured_block_interface.h"
//#include "mesh/structured_mesh.h"
#include "utils/neighbor_direction.h"
#include "utils/project_defs.h"
#include "disc_block.h"
#include "disc_interface.h"

namespace structured_fv {
namespace disc {


template <typename T>
class ElementField;

template <typename T>
using ElementFieldPtr = std::shared_ptr<ElementField<T>>;

template <typename T>
class VertexField;

template <typename T>
using VertexFieldPtr = std::shared_ptr<VertexField<T>>;

template <typename T>
class FaceField;

template <typename T>
using FaceFieldPtr = std::shared_ptr<FaceField<T>>;

class StructuredDisc
{
  public:
    StructuredDisc(std::shared_ptr<mesh::StructuredMesh> mesh, UInt num_ghost_cells);
            
    UInt getNumBlocks() const { return m_mesh->getNumBlocks(); }

    UInt getNumRegularBlocks() const { return m_mesh->getNumRegularBlocks(); }

    UInt getNumGhostBCBlocks() const { return m_mesh->getNumGhostBCBlocks(); }

    const StructuredBlock& getBlock(UInt i) const { return m_blocks[i]; }

    const StructuredBlock& getGhostBlock(UInt i) const { return m_blocks[i + getNumRegularBlocks()]; }

    UInt getNumBlockInterfaces() const { return m_mesh->getNumBlockInterfaces(); }

    UInt getNumRegularBlockInterfaces() const { return m_mesh->getNumRegularBlockInterfaces(); }

    UInt getNumGhostBCBlockInterfaces() const { return m_mesh->getNumGhostBCBlockInterfaces(); }

    const StructuredBlockInterface& getBlockInterface(UInt i) const { return m_ifaces.at(i); }

    VertexFieldPtr<Real> getCoordField() const { return m_coordField; }

    FaceFieldPtr<Real> getNormalField() const { return m_normalField; }

    ElementFieldPtr<Real> getInvCellVolumeField() const { return m_invCellVolumeField; }


  private:

    VertexFieldPtr<Real> createCoordField();

    ElementFieldPtr<Real> createInvCellVolumeField();

    FaceFieldPtr<Real> createFaceNormalField();

    std::shared_ptr<mesh::StructuredMesh> m_mesh;
    std::vector<StructuredBlock> m_blocks;
    std::vector<StructuredBlockInterface> m_ifaces;
    //UInt m_num_ghost_cells;
    VertexFieldPtr<Real> m_coordField;
    FaceFieldPtr<Real> m_normalField;
    ElementFieldPtr<Real> m_invCellVolumeField;

};


}
}

#endif