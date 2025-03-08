#include "structured_block_interface.h"
#include "mesh/adjacent_block_indexer.h"
#include "structured_block.h"
#include <sstream>

namespace structured_fv {
namespace mesh {

// rangeL and ranger give corresponding the ranges of indices along the boundary
// this allows supporting T junctions
StructuredBlockInterface::StructuredBlockInterface(const StructuredBlock& blockL,
                                                   NeighborDirection dirL,
                                                   const Range& cellRangeL,
                                                   const std::array<Int, 2>& transformL,
                                                   const StructuredBlock& blockR,
                                                   const Range& cellRangeR) :
  m_left_block_id(blockL.getBlockId()),
  m_dirL(dirL),
  m_transformL(transformL),
  m_right_block_id(blockR.getBlockId()),
  m_dirR(getNeighborImage(dirL, transformL)),
  m_transformR(getInverseTransform(transformL))
{
  m_transformR = getInverseTransform(transformL);

  if (mesh::getNumOwnedCells(blockL, dirL) != mesh::getNumOwnedCells(blockR, m_dirR))
  {
    std::stringstream ss;
    ss << "cannot create interface between block " << blockL.getBlockId() << ", dir " << dirL
       << " and block " << blockR.getBlockId() << ", dir " << m_dirL << std::endl
       << "Left block has " << mesh::getNumOwnedCells(blockL, dirL) << " cells along the interface"
       << " while the right block has " << mesh::getNumOwnedCells(blockR, m_dirR) << std::endl;
    throw std::runtime_error(ss.str());
  }


  if (cellRangeL.size() != cellRangeR.size())
    throw std::runtime_error("boundary ranges on block interface do not agree");  

  Range vertRangeL(*cellRangeL.begin(), (*cellRangeR.end())+1);
  Range vertRangeR(*cellRangeR.begin(), (*cellRangeR.end())+1);

  auto [boundary_vertsL, vert_indexerL] = createAdjacentBlockIndexer(blockL.getOwnedVerts(), m_dirL, vertRangeL, m_transformL,
                                                                     blockR.getOwnedVerts(), m_dirR, vertRangeR);
  auto [boundary_vertsR, vert_indexerR] = createAdjacentBlockIndexer(blockR.getOwnedVerts(), m_dirR, vertRangeR, m_transformR,
                                                                     blockL.getOwnedVerts(), m_dirL, vertRangeL);
  m_boundary_vertsL = boundary_vertsL;
  m_vert_indexerL   = vert_indexerL;
  m_boundary_vertsR = boundary_vertsR;
  m_vert_indexerR   = vert_indexerR;

  auto [boundary_cellsL, cell_indexerL] = createAdjacentBlockIndexer(blockL.getOwnedCells(), m_dirL, cellRangeL, m_transformL,
                                                                     blockR.getOwnedCells(), m_dirR, cellRangeR);
  auto [boundary_cellsR, cell_indexerR] = createAdjacentBlockIndexer(blockR.getOwnedCells(), m_dirR, cellRangeR, m_transformR,
                                                                     blockL.getOwnedCells(), m_dirL, cellRangeL);
  m_boundary_cellsL = boundary_cellsL;
  m_cell_indexerL   = cell_indexerL;
  m_boundary_cellsR = boundary_cellsR;
  m_cell_indexerR   = cell_indexerR;
}

std::pair<Range2D, AdjacentBlockIndexer> StructuredBlockInterface::createAdjacentBlockIndexer(
   const Range2D& block_rangeL, NeighborDirection dirL, 
   const Range& rangeL, const std::array<Int, 2>& transformL,
   const Range2D& block_rangeR, NeighborDirection dirR, const Range& rangeR)
{
  Range2D boundary_entities = getBoundaryRange(block_rangeL, dirL, rangeL);
  std::array<UInt, 2> min_cellL{*boundary_entities.getXRange().begin(),
                                *boundary_entities.getYRange().begin()};

  UInt variable_indexR = getMinEntityOnBoundary(dirL, transformL, rangeR);
  UInt constant_indexR = getConstantIndexAlongBoundary(block_rangeR, dirR);
  std::array<UInt, 2> min_entityR;
  min_entityR[0] = to_int(dirR) % 2 == 0 ? variable_indexR : constant_indexR;
  min_entityR[1] = to_int(dirR) % 2 == 0 ? constant_indexR : variable_indexR;

  AdjacentBlockIndexer indexer(transformL, min_cellL, dirL, min_entityR);

  return {boundary_entities, indexer};
}
 
}
}