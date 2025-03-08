#include "discretization.h"

namespace structured_fv {
namespace disc {

StructuredDisc::StructuredDisc(std::shared_ptr<mesh::StructuredMesh> mesh, UInt num_ghost_cells) :
  m_mesh(mesh)
  //m_num_ghost_cells(num_ghost_cells)
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

  //m_coordField = createCoordField();
  //m_normalField = createNormalField();
}

/*
std::shared_ptr<ElementField<Real>> StructuredDisc::createCoordField()
{
  auto coordField = ElementField<Real>(*this, 2);
  for (UInt b=0; b < getNumBlocks(); ++b)
  {
    
  }

}
*/



}
}