#ifndef STRUCTURED_FINITE_VOLUME_DISC_INTERFACE_H
#define STRUCTURED_FINITE_VOLUME_DISC_INTERFACE_H

#include "utils/project_defs.h"
#include "mesh/structured_block_interface.h"
#include "disc_block.h"

namespace structured_fv {
namespace disc {


class StructuredBlockInterface
{
  public:
    StructuredBlockInterface(const StructuredBlock& blockL, const StructuredBlock& blockR,
                             const mesh::StructuredBlockInterface& mesh_iface) :
      m_blockL(blockL),
      m_blockR(blockR),
      m_mesh_iface(mesh_iface)
    {
      std::cout << "\nConstructing disc interface" << std::endl;

      auto num_ghost_cells_per_directionL = blockL.getNumGhostCellsPerDirection();
      auto num_ghost_cells_per_directionR = blockR.getNumGhostCellsPerDirection();
      std::array<UInt, 2> left_block_min_cell = {*mesh_iface.getBoundaryCellsL().getXRange().begin(),
                                                 *mesh_iface.getBoundaryCellsL().getYRange().begin()};
      std::array<UInt, 2> right_block_min_cell = {*mesh_iface.getBoundaryCellsR().getXRange().begin(),
                                                  *mesh_iface.getBoundaryCellsR().getYRange().begin()};                                                 
      std::array<UInt, 2> other_block_min_cellL = mesh_iface.getAdjacentBlockIndexerL().getRightBlockMinCell();
      std::array<UInt, 2> other_block_min_cellR = mesh_iface.getAdjacentBlockIndexerL().getRightBlockMinCell();

      std::array<Int, 2> offset_directions{to_int(mesh::NeighborDirection::West),
                                           to_int(mesh::NeighborDirection::South)};
      for (int i=0; i < 2; ++i)
      {
        Int dir = offset_directions[i];
        left_block_min_cell[i]   += num_ghost_cells_per_directionL[dir];
        other_block_min_cellL[i] += num_ghost_cells_per_directionR[dir];

        right_block_min_cell[i]  += num_ghost_cells_per_directionR[dir];
        other_block_min_cellR[i] += num_ghost_cells_per_directionL[dir];
      }
    
      m_indexerL = mesh::AdjacentBlockIndexer(mesh_iface.getTransformL(), left_block_min_cell,
                                              mesh_iface.getNeighborDirectionL(), other_block_min_cellL);
      m_indexerR = mesh::AdjacentBlockIndexer(mesh_iface.getTransformR(), right_block_min_cell,
                                              mesh_iface.getNeighborDirectionR(), other_block_min_cellR);

      std::array<UInt, 2> ownedCellsXL{*mesh_iface.getBoundaryCellsL().getXRange().begin() + num_ghost_cells_per_directionL[offset_directions[0]],
                                       *mesh_iface.getBoundaryCellsL().getXRange().end()   + num_ghost_cells_per_directionL[offset_directions[0]]};

      std::array<UInt, 2> ownedCellsYL{*mesh_iface.getBoundaryCellsL().getYRange().begin() + num_ghost_cells_per_directionL[offset_directions[1]] ,
                                       *mesh_iface.getBoundaryCellsL().getYRange().end()   + num_ghost_cells_per_directionL[offset_directions[1]]};

      std::array<UInt, 2> ownedCellsXR{*mesh_iface.getBoundaryCellsR().getXRange().begin() + num_ghost_cells_per_directionR[offset_directions[0]] ,
                                       *mesh_iface.getBoundaryCellsR().getXRange().end()   + num_ghost_cells_per_directionR[offset_directions[0]]};

      std::array<UInt, 2> ownedCellsYR{*mesh_iface.getBoundaryCellsR().getYRange().begin() + num_ghost_cells_per_directionR[offset_directions[1]],
                                       *mesh_iface.getBoundaryCellsR().getYRange().end()   + num_ghost_cells_per_directionR[offset_directions[1]]};

      std::cout << "mesh::getBoundaryCellsR = " << mesh_iface.getBoundaryCellsR() << std::endl;
      std::cout << "ownedCellsXR = " << ownedCellsXR[0] << ", " << ownedCellsXR[1] << std::endl;
      m_owned_boundary_cellsL = Range2D(ownedCellsXL[0], ownedCellsXL[1], ownedCellsYL[0], ownedCellsYL[1]);
      m_owned_boundary_cellsR = Range2D(ownedCellsXR[0], ownedCellsXR[1], ownedCellsYR[0], ownedCellsYR[1]);                                                    
    }

    UInt getBlockIdL() const { return m_blockL.getBlockId(); }

    UInt getBlockIdR() const { return m_blockR.getBlockId(); }

    const Range2D& getOwnedBoundaryCellsL() const { return m_owned_boundary_cellsL; }

    const Range2D& getOwnedBoundaryCellsR() const { return m_owned_boundary_cellsR; }

    mesh::NeighborDirection getNeighborDirectionL() const { return m_mesh_iface.getNeighborDirectionL(); }

    mesh::NeighborDirection getNeighborDirectionR() const { return m_mesh_iface.getNeighborDirectionR(); }

    const mesh::AdjacentBlockIndexer& getAdjacentBlockIndexerL() const
    {
      return m_indexerL;
    }

    const mesh::AdjacentBlockIndexer& getAdjacentBlockIndexerR() const
    {
      return m_indexerR;
    }

  private:
    const StructuredBlock& m_blockL;
    const StructuredBlock& m_blockR;
    const mesh::StructuredBlockInterface m_mesh_iface;

    Range2D m_owned_boundary_cellsL;
    Range2D m_owned_boundary_cellsR;
    mesh::AdjacentBlockIndexer m_indexerL;
    mesh::AdjacentBlockIndexer m_indexerR;
};

}
}
#endif