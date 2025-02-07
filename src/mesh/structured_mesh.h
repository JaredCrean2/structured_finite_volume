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

    UInt getNumBlocks() const { return m_block_counts[0] + m_block_counts[1]; }

    UInt getNumRegularBlocks() const { return m_block_counts[0]; }

    UInt getNumGhostBCBlocks() const { return m_block_counts[1]; }

    const StructuredBlock& getBlock(UInt block) const
    {
      return m_blocks.at(block);
    }

    const StructuredBlock& getRegularBlock(UInt block) const
    {
      return m_blocks.at(block);
    }

    const StructuredBlock& getGhostBCBlock(UInt bc_block) const
    {
      return m_blocks.at(getNumRegularBlocks() + bc_block);
    }    

    UInt getNumBlockInterfaces() const { return m_block_iface_counts[0] + m_block_iface_counts[1]; }

    UInt getNumRegularBlockInterfaces() const { return m_block_iface_counts[0]; }

    UInt getNumGhostBCBlockInterfaces() const { return m_block_iface_counts[1]; }

    const StructuredBlockInterface& getBlockInterface(UInt iface) const
    {
      return m_block_interfaces.at(iface);
    }

    const StructuredBlockInterface& getRegularBlockInterface(UInt iface) const
    {
      return m_block_interfaces.at(iface);
    }

    const StructuredBlockInterface& getGhostBCBlockInterface(UInt iface) const
    {
      return m_block_interfaces.at(getNumRegularBlockInterfaces() + iface);
    }

    UInt getNumBCs() const { return m_bc_block_ranges.size(); }

    // gives the range of ghost blocks that are part of the given BC
    Range getBCBlockRange(UInt bc) const { return m_bc_block_ranges.at(bc); }

    // gives the range of ghost block interfaces that are part of the BC
    // (ie. one of the blocks in a ghost BC block)
    Range getBCInterfaceRange(UInt bc) const { return m_bc_iface_ranges.at(bc); }

    std::array<Int, 4> getBlockInterfaces(UInt block) const { return m_block_interface_connectivity.at(block); }

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

