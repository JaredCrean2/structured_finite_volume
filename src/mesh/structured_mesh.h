#ifndef STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_MESH_H
#define STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_MESH_H

#include "structured_block.h"
#include "structured_block_interface.h"
#include "block_spec.h"

namespace structured_fv {
namespace mesh {

class StructuredMesh
{
  public:
    // creates a mesh with 4 BCs, in the North, East, South, and West direction
    StructuredMesh(const MeshSpec& spec);

    UInt getNumBlocks() const;

    UInt getNumRegularBlocks() const;

    UInt getNumGhostBCBlocks() const;

    const StructuredBlock& getBlock(UInt block) const;

    const StructuredBlock& getRegularBlock(UInt block) const;

    const StructuredBlock& getGhostBCBlock(UInt bc_block) const;  

    UInt getNumBlockInterfaces() const;

    UInt getNumRegularBlockInterfaces() const;

    UInt getNumGhostBCBlockInterfaces() const;

    const StructuredBlockInterface& getBlockInterface(UInt iface) const;

    const StructuredBlockInterface& getRegularBlockInterface(UInt iface) const;

    const StructuredBlockInterface& getGhostBCBlockInterface(UInt iface) const;

    UInt getNumBCs() const;

    // gives the range of ghost blocks that are part of the given BC
    Range getBCBlockRange(UInt bc) const;

    // gives the range of ghost block interfaces that are part of the BC
    // (ie. one of the blocks in a ghost BC block)
    Range getBCInterfaceRange(UInt bc) const;

    std::array<Int, 4> getBlockInterfaces(UInt block) const;

  private:
    void createBCGhosts(const MeshSpec& spec);

    void createBCGhost(const MeshBlockSpec& spec, UInt regular_block_id, NeighborDirection domain_boundary);


    std::array<UInt, 2> m_block_counts{0, 0};
    std::array<UInt, 2> m_block_iface_counts{0, 0};
    std::vector<Range> m_bc_block_ranges;
    std::vector<Range> m_bc_iface_ranges;
    std::vector<StructuredBlock> m_blocks;
    std::vector<StructuredBlockInterface> m_block_interfaces;
    std::vector<std::array<Int, 4>> m_block_interface_connectivity;
};

}
}

#endif

