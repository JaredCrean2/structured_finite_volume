#include "disc_block.h"

namespace structured_fv {
namespace disc {


StructuredBlock::StructuredBlock(const mesh::StructuredMesh& mesh, const mesh::StructuredBlock& block, int nghost) :
  m_mesh_block(block),
  m_num_ghost_cells_per_direction(computeNumGhostsPerDirection(mesh, nghost))
{}

UInt StructuredBlock::getBlockId() const { return m_mesh_block.getBlockId(); }

mesh::BlockType StructuredBlock::getBlockType() const { return m_mesh_block.getBlockType(); }

std::array<UInt, 2> StructuredBlock::getCellDimensions() const { return getOwnedAndGhostCells().getDimensions(); }

std::array<UInt, 2> StructuredBlock::getVertDimensions() const { return getOwnedAndGhostVerts().getDimensions(); }

const std::array<UInt, 4> StructuredBlock::getNumGhostCellsPerDirection() const { return m_num_ghost_cells_per_direction; }

Range2D StructuredBlock::getOwnedVerts() const
{
  UInt xoffset = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::West)];
  UInt yoffset = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::South)];
  Range2D mesh_range = m_mesh_block.getOwnedVerts();
  return Range2D(*mesh_range.getXRange().begin() + xoffset,
                  *mesh_range.getXRange().end()   + xoffset,
                  *mesh_range.getYRange().begin() + yoffset,
                  *mesh_range.getYRange().end()   + yoffset);
}

Range2D StructuredBlock::getOwnedCells() const
{
  UInt xoffset = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::West)];
  UInt yoffset = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::South)];
  Range2D mesh_range = m_mesh_block.getOwnedCells();
  return Range2D(*mesh_range.getXRange().begin() + xoffset,
                  *mesh_range.getXRange().end()   + xoffset,
                  *mesh_range.getYRange().begin() + yoffset,
                  *mesh_range.getYRange().end()   + yoffset);      
}

Range2D StructuredBlock::getOwnedAndGhostVerts() const
{
  UInt extra_x_verts = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::West)] +
                        m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::East)];
  UInt extra_y_verts = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::South)] +
                        m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::North)];
  Range2D mesh_range = m_mesh_block.getOwnedVerts();
  return Range2D(0, mesh_range.getXRange().size() + extra_x_verts,
                  0, mesh_range.getYRange().size() + extra_y_verts);
}

Range2D StructuredBlock::getOwnedAndGhostCells() const
{
  UInt extra_x_cells = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::West)] +
                        m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::East)];
  UInt extra_y_cells = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::South)] +
                        m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::North)];
  Range2D mesh_range = m_mesh_block.getOwnedCells();
  return Range2D(0, mesh_range.getXRange().size() + extra_x_cells,
                  0, mesh_range.getYRange().size() + extra_y_cells);
}

std::pair<UInt, UInt> StructuredBlock::meshVertToBlockVert(UInt i, UInt j) const
{
  UInt xoffset = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::West)];
  UInt yoffset = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::South)];
  return {i + xoffset, j + yoffset};
}

std::pair<UInt, UInt> StructuredBlock::blockVertToMeshVert(UInt i, UInt j) const
{
  assert(in(getOwnedVerts(), i, j));
  UInt xoffset = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::West)];
  UInt yoffset = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::South)];      
  return {i - xoffset, j - yoffset};
}

std::array<UInt, 4> StructuredBlock::computeNumGhostsPerDirection(const mesh::StructuredMesh& mesh, int nghost)
{
  std::array<UInt, 4> num_ghost_cells_per_direction;
  const auto& block_iface_idxs = mesh.getBlockInterfaces(m_mesh_block.getBlockId());
  for (int i=0; i < 4; ++i)
  {
    if (block_iface_idxs[i] >= 0)
    {
      num_ghost_cells_per_direction[i] = nghost;
      const mesh::StructuredBlockInterface& iface = mesh.getBlockInterface(block_iface_idxs[i]);
      const mesh::StructuredBlock& blockL = mesh.getBlock(iface.getBlockIdL());
      const mesh::StructuredBlock& blockR = mesh.getBlock(iface.getBlockIdR());

      if (blockL.getOwnedCells().getDimensions()[mesh::to_int(iface.getNeighborDirectionL()) % 2] < nghost)
        throw std::runtime_error("block is too small for given number of ghosts");

      if (blockR.getOwnedCells().getDimensions()[mesh::to_int(iface.getNeighborDirectionR()) % 2] < nghost)
        throw std::runtime_error("block is too small for given number of ghosts");
    } else
    {
      num_ghost_cells_per_direction[i] = 0;
    }
  }

  return num_ghost_cells_per_direction;
}

}
}