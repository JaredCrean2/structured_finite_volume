#ifndef STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_BLOCK_H
#define STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_BLOCK_H

#include "utils/neighbor_direction.h"
#include "utils/project_defs.h"
#include "utils/range.h"
#include "block_spec.h"

namespace structured_fv {
namespace mesh {

enum class BlockType
{
  Regular,
  GhostBC
};


class StructuredBlock
{
  public:
    //TODO: the datatype should be const
    using CoordsHostView = Kokkos::View<Real**[2], HostMemorySpace>;
    
    StructuredBlock(const MeshBlockSpec& spec, UInt block_id, BlockType block_type);

    StructuredBlock(Kokkos::View<Real**[2], HostMemorySpace> ghost_coords, UInt block_id);

    StructuredBlock(const StructuredBlock& rhs) = delete;

    StructuredBlock(StructuredBlock&& rhs) = default;

    StructuredBlock& operator=(const StructuredBlock& rhs) = delete;

    StructuredBlock& operator=(StructuredBlock&& rhs) = default;

    UInt getBlockId() const;// block IDs are global

    BlockType getBlockType() const;

    // these are the verts of owned cells
    // slightly misleading, because we dont define ownership of verts
    Range2D getOwnedVerts() const;

    Range2D getOwnedCells() const;

    const CoordsHostView& getOwnedVertCoords() const;

    std::array<UInt, 2> getOffsetIntoBlock() const; // owned idx + offset = block idx

    // returns size of the entire block (including non-owned cells)
    std::array<UInt, 2> getAllBlockSize() const;

  private:
    UInt m_block_id;
    BlockType m_block_type;
    Range2D m_owned_cell_range;
    Range2D m_owned_vert_range;
    CoordsHostView m_owned_vert_coords;
    std::array<UInt, 2> m_offset_into_block;
    std::array<UInt, 2> m_all_block_size;
};

UInt getNumOwnedCells(const StructuredBlock& block, NeighborDirection dir);


}
}

#endif

