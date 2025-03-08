#ifndef STRUCTURED_FINITE_VOLUME_DISC_BLOCK_H
#define STRUCTURED_FINITE_VOLUME_DISC_BLOCK_H

#include "utils/project_defs.h"
#include "mesh/structured_mesh.h"

namespace structured_fv {
namespace disc {


// A rectangular block of cells + n layers of ghost cells in
// the cardinal directions.  Currently ghost cells to the
// diagonal are not supported
class StructuredBlock
{
  public:
    using CoordsHostView = mesh::StructuredBlock::CoordsHostView;

    StructuredBlock(const mesh::StructuredMesh& mesh, const mesh::StructuredBlock& block, int nghost) :
      m_mesh_block(block),
      m_num_ghost_cells_per_direction{0, 0, 0, 0}
    {
      setNumGhostsPerDirection(mesh, nghost);
      //setVertCoords(mesh, nghost);
    }

    UInt getBlockId() const { return m_mesh_block.getBlockId(); }

    mesh::BlockType getBlockType() const { return m_mesh_block.getBlockType(); }

    std::array<UInt, 2> getCellDimensions() const { return getOwnedAndGhostCells().getDimensions(); }

    std::array<UInt, 2> getVertDimensions() const { return getOwnedAndGhostVerts().getDimensions(); }

    const std::array<UInt, 4> getNumGhostCellsPerDirection() const { return m_num_ghost_cells_per_direction; }

    Range2D getOwnedVerts() const
    {
      UInt xoffset = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::West)];
      UInt yoffset = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::South)];
      Range2D mesh_range = m_mesh_block.getOwnedVerts();
      return Range2D(*mesh_range.getXRange().begin() + xoffset,
                     *mesh_range.getXRange().end()   + xoffset,
                     *mesh_range.getYRange().begin() + yoffset,
                     *mesh_range.getYRange().end()   + yoffset);
    }

    Range2D getOwnedCells() const
    {
      UInt xoffset = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::West)];
      UInt yoffset = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::South)];
      Range2D mesh_range = m_mesh_block.getOwnedCells();
      return Range2D(*mesh_range.getXRange().begin() + xoffset,
                     *mesh_range.getXRange().end()   + xoffset,
                     *mesh_range.getYRange().begin() + yoffset,
                     *mesh_range.getYRange().end()   + yoffset);      
    }

    Range2D getOwnedAndGhostVerts() const
    {
      UInt extra_x_verts = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::West)] +
                           m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::East)];
      UInt extra_y_verts = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::South)] +
                           m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::North)];
      Range2D mesh_range = m_mesh_block.getOwnedVerts();
      return Range2D(0, mesh_range.getXRange().size() + extra_x_verts,
                     0, mesh_range.getYRange().size() + extra_y_verts);
    }

    Range2D getOwnedAndGhostCells() const
    {
      UInt extra_x_cells = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::West)] +
                           m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::East)];
      UInt extra_y_cells = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::South)] +
                           m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::North)];
      Range2D mesh_range = m_mesh_block.getOwnedCells();
      return Range2D(0, mesh_range.getXRange().size() + extra_x_cells,
                     0, mesh_range.getYRange().size() + extra_y_cells);
    }

    std::pair<UInt, UInt> meshVertToBlockVert(UInt i, UInt j) const
    {
      UInt xoffset = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::West)];
      UInt yoffset = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::South)];
      return {i + xoffset, j + yoffset};
    }

    std::pair<UInt, UInt> blockVertToMeshVert(UInt i, UInt j) const
    {
      assert(in(getOwnedVerts(), i, j));
      UInt xoffset = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::West)];
      UInt yoffset = m_num_ghost_cells_per_direction[mesh::to_int(mesh::NeighborDirection::South)];      
      return {i - xoffset, j - yoffset};
    }

  private:

    void setNumGhostsPerDirection(const mesh::StructuredMesh& mesh, int nghost)
    {
      const auto& block_iface_idxs = mesh.getBlockInterfaces(m_mesh_block.getBlockId());
      for (int i=0; i < 4; ++i)
      {
        if (block_iface_idxs[i] >= 0)
        {
          m_num_ghost_cells_per_direction[i] = nghost;
          const mesh::StructuredBlockInterface& iface = mesh.getBlockInterface(block_iface_idxs[i]);
          const mesh::StructuredBlock& blockL = mesh.getBlock(iface.getBlockIdL());
          const mesh::StructuredBlock& blockR = mesh.getBlock(iface.getBlockIdR());

          if (blockL.getOwnedCells().getDimensions()[mesh::to_int(iface.getNeighborDirectionL()) % 2] < nghost)
            throw std::runtime_error("block is too small for given number of ghosts");

          if (blockR.getOwnedCells().getDimensions()[mesh::to_int(iface.getNeighborDirectionR()) % 2] < nghost)
            throw std::runtime_error("block is too small for given number of ghosts");
        } else
        {
          m_num_ghost_cells_per_direction[i] = 0;
        }
      }
    }
/*
    void setVertCoords(const mesh::StructuredMesh& mesh, int nghost)
    {
      Kokkos::resize(m_vert_coords, getOwnedAndGhostVerts().getXRange().size(),
                                    getOwnedAndGhostVerts().getYRange().size(), 2);
      const auto& meshVertCoords = m_mesh_block.getOwnedVertCoords();                                    
      for (UInt i : getOwnedVerts().getXRange())
        for (UInt j : getOwnedVerts().getYRange())
        {
          auto [imesh, jmesh] = blockVertToMeshVert(i, j);
          for (UInt d=0; d < 2; ++d)
          {
            m_vert_coords(i, j, d) = meshVertCoords(imesh, jmesh, d);
          }
        }

      for (Int iface_idx : mesh.getBlockInterfaces(m_mesh_block.getBlockId()))
      {
        if (iface_idx != -1)
        {
          const mesh::StructuredBlockInterface& iface = mesh.getBlockInterface(iface_idx);
          bool amLeftBlock = iface.getBlockIdL() == m_mesh_block.getBlockId();
          const mesh::StructuredBlock& blockR       = amLeftBlock ? mesh.getBlock(iface.getBlockIdR()) : 
                                                                    mesh.getBlock(iface.getBlockIdL());
          const mesh::AdjacentBlockIndexer& indexer = amLeftBlock ? iface.getAdjBlockCellIndexerL() : 
                                                                    iface.getAdjBlockCellIndexerR();
          const Range2D& boundary_cells             = amLeftBlock ? iface.getBoundaryCellsL() :
                                                                    iface.getBoundaryCellsR();
          mesh::NeighborDirection dir               = amLeftBlock ? iface.getNeighborDirectionL() : 
                                                                    iface.getNeighborDirectionR();
          const auto neighborCoords = blockR.getOwnedVertCoords();

          for (UInt imesh : boundary_cells.getXRange())
            for(UInt jmesh : boundary_cells.getYRange())
            {
              auto [iblock, jblock] = meshVertToBlockVert(imesh, jmesh);
              for (int n=1; n <= nghost; ++n)
              {
                auto [ighost, jghost] = mesh::computeIndices(dir, n, iblock, jblock);
                auto [ineighbor, jneighbor] = indexer(mesh::computeIndices(dir, n, imesh, jmesh));
                for (int d=0; d < 2; ++d)
                  m_vert_coords(ighost, jghost, d) = neighborCoords(ineighbor, jneighbor, d);
              }
            }
        }        
      }
    }
*/

    const mesh::StructuredBlock& m_mesh_block;
    std::array<UInt, 4> m_num_ghost_cells_per_direction;
};

inline std::array<Real, 2> computeCellCentroid(StructuredBlock::CoordsHostView vertCoords, int i, int j)
{
  return {(vertCoords(i, j, 0) + vertCoords(i+1, j, 0) + 
           vertCoords(i+1, j+1, 0) + vertCoords(i+1, j, 0))/4,
          (vertCoords(i, j, 1) + vertCoords(i+1, j, 1) + 
           vertCoords(i+1, j+1, 1) + vertCoords(i+1, j, 1))/4};
}


}
}

#endif
