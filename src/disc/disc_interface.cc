#include "disc_interface.h"

namespace structured_fv {
namespace disc {

StructuredBlockInterface::StructuredBlockInterface(const StructuredBlock& blockL, 
                                                   const StructuredBlock& blockR,
                                                   const mesh::StructuredBlockInterface& mesh_iface) :
  m_blockL(blockL),
  m_blockR(blockR),
  m_mesh_iface(mesh_iface)
{

  auto num_ghost_cells_per_directionL = blockL.getNumGhostCellsPerDirection();
  auto num_ghost_cells_per_directionR = blockR.getNumGhostCellsPerDirection();
  std::array<UInt, 2> left_block_min_vert = {*mesh_iface.getBoundaryVertsL().getXRange().begin(),
                                              *mesh_iface.getBoundaryVertsL().getYRange().begin()};
  std::array<UInt, 2> right_block_min_vert = {*mesh_iface.getBoundaryVertsR().getXRange().begin(),
                                              *mesh_iface.getBoundaryVertsR().getYRange().begin()};

  std::array<UInt, 2> left_block_min_cell = {*mesh_iface.getBoundaryCellsL().getXRange().begin(),
                                              *mesh_iface.getBoundaryCellsL().getYRange().begin()};
  std::array<UInt, 2> right_block_min_cell = {*mesh_iface.getBoundaryCellsR().getXRange().begin(),
                                              *mesh_iface.getBoundaryCellsR().getYRange().begin()}; 

  std::array<UInt, 2> other_block_min_vertL = mesh_iface.getAdjBlockVertIndexerL().getRightBlockMinEntity();
  std::array<UInt, 2> other_block_min_vertR = mesh_iface.getAdjBlockVertIndexerR().getRightBlockMinEntity();

  std::array<UInt, 2> other_block_min_cellL = mesh_iface.getAdjBlockCellIndexerL().getRightBlockMinEntity();
  std::array<UInt, 2> other_block_min_cellR = mesh_iface.getAdjBlockCellIndexerR().getRightBlockMinEntity();

  std::array<Int, 2> offset_directions{to_int(NeighborDirection::West),
                                        to_int(NeighborDirection::South)};
  for (int i=0; i < 2; ++i)
  {
    Int dir = offset_directions[i];
    left_block_min_vert[i]   += num_ghost_cells_per_directionL[dir];
    other_block_min_vertL[i] += num_ghost_cells_per_directionR[dir];        
    left_block_min_cell[i]   += num_ghost_cells_per_directionL[dir];
    other_block_min_cellL[i] += num_ghost_cells_per_directionR[dir];

    right_block_min_vert[i]  += num_ghost_cells_per_directionR[dir];
    other_block_min_vertR[i] += num_ghost_cells_per_directionL[dir];
    right_block_min_cell[i]  += num_ghost_cells_per_directionR[dir];
    other_block_min_cellR[i] += num_ghost_cells_per_directionL[dir];
  }

  m_vert_indexerL = mesh::AdjacentBlockIndexer(mesh_iface.getTransformL(), left_block_min_vert,
                                                mesh_iface.getNeighborDirectionL(), other_block_min_vertL);
  m_vert_indexerR = mesh::AdjacentBlockIndexer(mesh_iface.getTransformR(), right_block_min_vert,
                                                mesh_iface.getNeighborDirectionR(), other_block_min_vertR);

  m_cell_indexerL = mesh::AdjacentBlockIndexer(mesh_iface.getTransformL(), left_block_min_cell,
                                                mesh_iface.getNeighborDirectionL(), other_block_min_cellL);
  m_cell_indexerR = mesh::AdjacentBlockIndexer(mesh_iface.getTransformR(), right_block_min_cell,
                                                mesh_iface.getNeighborDirectionR(), other_block_min_cellR);


  std::array<UInt, 2> ownedVertsXL{*mesh_iface.getBoundaryVertsL().getXRange().begin() + num_ghost_cells_per_directionL[offset_directions[0]],
                                    *mesh_iface.getBoundaryVertsL().getXRange().end()   + num_ghost_cells_per_directionL[offset_directions[0]]};

  std::array<UInt, 2> ownedVertsYL{*mesh_iface.getBoundaryVertsL().getYRange().begin() + num_ghost_cells_per_directionL[offset_directions[1]] ,
                                    *mesh_iface.getBoundaryVertsL().getYRange().end()   + num_ghost_cells_per_directionL[offset_directions[1]]};

  std::array<UInt, 2> ownedVertsXR{*mesh_iface.getBoundaryVertsR().getXRange().begin() + num_ghost_cells_per_directionR[offset_directions[0]] ,
                                    *mesh_iface.getBoundaryVertsR().getXRange().end()   + num_ghost_cells_per_directionR[offset_directions[0]]};

  std::array<UInt, 2> ownedVertsYR{*mesh_iface.getBoundaryVertsR().getYRange().begin() + num_ghost_cells_per_directionR[offset_directions[1]],
                                    *mesh_iface.getBoundaryVertsR().getYRange().end()   + num_ghost_cells_per_directionR[offset_directions[1]]};


  std::array<UInt, 2> ownedCellsXL{*mesh_iface.getBoundaryCellsL().getXRange().begin() + num_ghost_cells_per_directionL[offset_directions[0]],
                                    *mesh_iface.getBoundaryCellsL().getXRange().end()   + num_ghost_cells_per_directionL[offset_directions[0]]};

  std::array<UInt, 2> ownedCellsYL{*mesh_iface.getBoundaryCellsL().getYRange().begin() + num_ghost_cells_per_directionL[offset_directions[1]] ,
                                    *mesh_iface.getBoundaryCellsL().getYRange().end()   + num_ghost_cells_per_directionL[offset_directions[1]]};

  std::array<UInt, 2> ownedCellsXR{*mesh_iface.getBoundaryCellsR().getXRange().begin() + num_ghost_cells_per_directionR[offset_directions[0]] ,
                                    *mesh_iface.getBoundaryCellsR().getXRange().end()   + num_ghost_cells_per_directionR[offset_directions[0]]};

  std::array<UInt, 2> ownedCellsYR{*mesh_iface.getBoundaryCellsR().getYRange().begin() + num_ghost_cells_per_directionR[offset_directions[1]],
                                    *mesh_iface.getBoundaryCellsR().getYRange().end()   + num_ghost_cells_per_directionR[offset_directions[1]]};

  m_owned_boundary_vertsL = Range2D(ownedVertsXL[0], ownedVertsXL[1], ownedVertsYL[0], ownedVertsYL[1]);
  m_owned_boundary_vertsR = Range2D(ownedVertsXR[0], ownedVertsXR[1], ownedVertsYR[0], ownedVertsYR[1]);

  m_owned_boundary_cellsL = Range2D(ownedCellsXL[0], ownedCellsXL[1], ownedCellsYL[0], ownedCellsYL[1]);
  m_owned_boundary_cellsR = Range2D(ownedCellsXR[0], ownedCellsXR[1], ownedCellsYR[0], ownedCellsYR[1]);                                                    
}

UInt StructuredBlockInterface::getBlockIdL() const { return m_blockL.getBlockId(); }

UInt StructuredBlockInterface::getBlockIdR() const { return m_blockR.getBlockId(); }

const Range2D& StructuredBlockInterface::getOwnedBoundaryVertsL() const { return m_owned_boundary_vertsL; }

const Range2D& StructuredBlockInterface::getOwnedBoundaryVertsR() const { return m_owned_boundary_vertsR; }

const Range2D& StructuredBlockInterface::getOwnedBoundaryCellsL() const { return m_owned_boundary_cellsL; }

const Range2D& StructuredBlockInterface::getOwnedBoundaryCellsR() const { return m_owned_boundary_cellsR; }

NeighborDirection StructuredBlockInterface::getNeighborDirectionL() const { return m_mesh_iface.getNeighborDirectionL(); }

NeighborDirection StructuredBlockInterface::getNeighborDirectionR() const { return m_mesh_iface.getNeighborDirectionR(); }

const mesh::AdjacentBlockIndexer& StructuredBlockInterface::getAdjBlockVertIndexerL() const { return m_vert_indexerL; }

const mesh::AdjacentBlockIndexer& StructuredBlockInterface::getAdjBlockVertIndexerR() const { return m_vert_indexerR; }

const mesh::AdjacentBlockIndexer& StructuredBlockInterface::getAdjBlockCellIndexerL() const { return m_cell_indexerL; }

const mesh::AdjacentBlockIndexer& StructuredBlockInterface::getAdjBlockCellIndexerR() const { return m_cell_indexerR; }

}
}