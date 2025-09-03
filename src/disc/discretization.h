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
    StructuredDisc(std::shared_ptr<mesh::StructuredMesh> mesh, UInt num_ghost_cells, int dofs_per_cell);

    StructuredDisc(const StructuredDisc& other) = delete;

    StructuredDisc& operator=(const StructuredDisc& other) = delete;
            
    UInt getNumBlocks() const { return m_mesh->getNumBlocks(); }

    UInt getNumRegularBlocks() const { return m_mesh->getNumRegularBlocks(); }

    UInt getNumGhostBCBlocks() const { return m_mesh->getNumGhostBCBlocks(); }

    const StructuredBlock& getBlock(UInt i) const { return m_blocks[i]; }

    const StructuredBlock& getGhostBlock(UInt i) const { return m_blocks[i + getNumRegularBlocks()]; }

    Range getRegularBlocksIds() const { return Range(0, getNumRegularBlocks()); }

    Range getGhostBCBlocksIds() const { return Range(getNumRegularBlocks(), getNumBlocks()); }

    Range getAllBlocksIds() const { return Range(0, getNumBlocks()); }

    UInt getNumBlockInterfaces() const { return m_mesh->getNumBlockInterfaces(); }

    UInt getNumRegularBlockInterfaces() const { return m_mesh->getNumRegularBlockInterfaces(); }

    UInt getNumGhostBCBlockInterfaces() const { return m_mesh->getNumGhostBCBlockInterfaces(); }

    Range getRegularBlockInterfacesIds() const { return Range(0, getNumRegularBlockInterfaces()); }

    Range getBCBlockInterfacesIds() const { return Range(getNumRegularBlockInterfaces(), getNumBlockInterfaces()); }

    Range getAllBlockInterfacesIds() const { return Range(0, getNumBlockInterfaces()); }

    const StructuredBlockInterface& getBlockInterface(UInt i) const { return m_ifaces.at(i); }

    FixedVec<Int, 4> getBlockInterfaces(UInt block_id) const { return m_mesh->getBlockInterfaces(block_id); }

    VertexFieldPtr<Real> getCoordField() const { return m_coordField; }

    FaceFieldPtr<Real> getNormalField() const { return m_normalField; }

    ElementFieldPtr<Real> getInvCellVolumeField() const { return m_invCellVolumeField; }

    ElementFieldPtr<GlobalDof> getDofNumbering() const { return m_dofNumbering; }

    GlobalDof getNumDofs() const { return m_num_dofs; }

    UInt getNumDofsPerNode() const { return m_dofs_per_cell; }

    UInt getNumGhostCells() const { return m_num_ghost_cells; }

  private:

    VertexFieldPtr<Real> createCoordField();

    ElementFieldPtr<Real> createInvCellVolumeField();

    FaceFieldPtr<Real> createFaceNormalField();

    ElementFieldPtr<GlobalDof> createDofNumbering(UInt dofs_per_cell, GlobalDof& num_dofs);

    UInt m_dofs_per_cell;
    UInt m_num_ghost_cells;
    std::shared_ptr<mesh::StructuredMesh> m_mesh;
    std::vector<StructuredBlock> m_blocks;
    std::vector<StructuredBlockInterface> m_ifaces;
    //UInt m_num_ghost_cells;
    VertexFieldPtr<Real> m_coordField;
    FaceFieldPtr<Real> m_normalField;
    ElementFieldPtr<Real> m_invCellVolumeField;
    ElementFieldPtr<GlobalDof> m_dofNumbering;
    GlobalDof m_num_dofs;
};

using StructuredDiscPtr = std::shared_ptr<StructuredDisc>;


}
}

#endif