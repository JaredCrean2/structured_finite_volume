#include "discretization.h"
#include "disc/disc_block.h"

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

  m_coordField = createCoordField();
  //m_normalField = createNormalField();
}


VertexFieldPtr<Real> StructuredDisc::createCoordField()
{
  auto coordField = std::make_shared<VertexField<Real>>(*this, 2);
  for (UInt i=0; i < getNumBlocks(); ++i)
  {
    const StructuredBlock& block = getBlock(i);
    const mesh::StructuredBlock& mesh_block = m_mesh->getBlock(i);
    const auto& meshVertCoords = mesh_block.getOwnedVertCoords();
    auto& blockVertCoords = coordField->getData(i);
    for (UInt i : block.getOwnedVerts().getXRange())
      for (UInt j : block.getOwnedVerts().getYRange())
      {
        auto [imesh, jmesh] = block.blockVertToMeshVert(i, j);
        for (UInt d=0; d < 2; ++d)
        {
          blockVertCoords(i, j, d) = meshVertCoords(imesh, jmesh, d);
        }
      }
  }

  coordField->updateGhostValues();

  return coordField;
}


}
}