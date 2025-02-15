#ifndef STRUCTURED_FINITE_VOLUME_DISC_DISCRETIZATION_H
#define STRUCTURED_FINITE_VOLUME_DISC_DISCRETIZATION_H

//#include "mesh/structured_block.h"
//#include "mesh/structured_block_interface.h"
//#include "mesh/structured_mesh.h"
#include "utils/project_defs.h"
#include "disc_block.h"
#include "disc_interface.h"

namespace structured_fv {
namespace disc {


template <typename T>
class ElementField;


class StructuredDisc
{
  public:
    StructuredDisc(std::shared_ptr<mesh::StructuredMesh> mesh, UInt num_ghost_cells) :
      m_mesh(mesh),
      m_num_ghost_cells(num_ghost_cells)
    {
      for (UInt i=0; i < mesh->getNumBlocks(); ++i)
      {
        m_blocks.emplace_back(*mesh, mesh->getBlock(i), num_ghost_cells);
      }

      for (UInt i=0; i < mesh->getNumBlockInterfaces(); ++i)
      {
        const mesh::StructuredBlockInterface& mesh_iface = mesh->getBlockInterface(i);
        m_ifaces.emplace_back(m_blocks[mesh_iface.getBlockIdL()], m_blocks[mesh_iface.getBlockIdR()], mesh_iface);
      }
    }
            
    UInt getNumBlocks() const { return m_mesh->getNumBlocks(); }

    UInt getNumRegularBlocks() const { return m_mesh->getNumRegularBlocks(); }

    UInt getNumGhostBCBlocks() const { return m_mesh->getNumGhostBCBlocks(); }

    const StructuredBlock& getBlock(UInt i) const { return m_blocks[i]; }

    const StructuredBlock& getGhostBlock(UInt i) const { return m_blocks[i + getNumRegularBlocks()]; }

    UInt getNumBlockInterfaces() const { return m_mesh->getNumBlockInterfaces(); }

    UInt getNumRegularBlockInterfaces() const { return m_mesh->getNumRegularBlockInterfaces(); }

    UInt getNumGhostBCBlockInterfaces() const { return m_mesh->getNumGhostBCBlockInterfaces(); }

    const StructuredBlockInterface& getBlockInterface(UInt i) const { return m_ifaces.at(i); }

    // overwrite ghost values with the owner values
    template <typename T>
    void updateGhostValues(ElementField<T>& field);

  private:
    std::shared_ptr<mesh::StructuredMesh> m_mesh;
    std::vector<StructuredBlock> m_blocks;
    std::vector<StructuredBlockInterface> m_ifaces;
    UInt m_num_ghost_cells;
};

// defines a field that stores n values in each element of each block (including ghosts)
template <typename T>
class ElementField
{
  public:
    using FieldData = Kokkos::View<T***, HostMemorySpace>;

    ElementField(std::shared_ptr<StructuredDisc> disc, Int nvals_per_element)
    {
      for (UInt i=0; i < disc->getNumBlocks(); ++i)
      {
        const StructuredBlock& block = disc->getBlock(i);
        auto dimensions = block.getCellDimensions();
        m_data.emplace_back("field_data", dimensions[0], dimensions[1], nvals_per_element);
      }
    }

    FieldData& getData(UInt block) { return m_data[block]; }

    FieldData& getData(UInt block) const { return m_data[block]; }

    T& operator()(UInt block, UInt i, UInt j, UInt v) { return m_data[block](i, j, v); }      

    T& operator()(UInt block, UInt i, UInt j, UInt v) const { return m_data[block](i, j, v); }      

  private:
    std::vector<FieldData> m_data;
};

template <typename T>
void StructuredDisc::updateGhostValues(ElementField<T>& field)
{
  for (UInt i=0; i < getNumBlockInterfaces(); ++i)
  {
    const StructuredBlockInterface& iface = getBlockInterface(i);
    const auto& indexerL = iface.getAdjacentBlockIndexerL();
    const auto& indexerR = iface.getAdjacentBlockIndexerR();
    auto fieldL = field.getData(iface.getBlockIdL());
    auto fieldR = field.getData(iface.getBlockIdR());
    
    for (UInt i : iface.getOwnedBoundaryCellsL().getXRange())
      for (UInt j : iface.getOwnedBoundaryCellsL().getXRange())
        for (int v=1; v <= m_num_ghost_cells; ++v)
        {
          auto [iprime, jprime] = mesh::computeIndices(iface.getNeighborDirectionL(), v, i, j);
          auto [ineighbor, jneighbor] = indexerL(iprime, jprime);

          for (UInt d=0; d < fieldL.extent(2); ++d)
          {
            fieldL(iprime, jprime, d) = fieldR(ineighbor, jneighbor, d);
          }      
        }

    //TODO: it would have better temporal locality to merge this
    //      into the above loops
    for (UInt i : iface.getOwnedBoundaryCellsR().getXRange())
      for (UInt j : iface.getOwnedBoundaryCellsR().getXRange())
        for (int v=1; v <= m_num_ghost_cells; ++v)
        {
          auto [iprime, jprime] = mesh::computeIndices(iface.getNeighborDirectionR(), v, i, j);
          auto [ineighbor, jneighbor] = indexerR(iprime, jprime);

          for (UInt d=0; d < fieldL.extent(2); ++d)
          {
            fieldR(iprime, jprime, d) = fieldL(ineighbor, jneighbor, d);
          }      
        }        
  }
}

}
}

#endif