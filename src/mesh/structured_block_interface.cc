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
                                                   const Range& rangeL,
                                                   const std::array<Int, 2>& transformL,
                                                   const StructuredBlock& blockR,
                                                   const Range& rangeR) :
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


  if (rangeL.size() != rangeR.size())
    throw std::runtime_error("boundary ranges on block interface do not agree");  

  auto [boundary_cellsL, indexerL] = createAdjacentBlockIndexer(blockL, m_dirL, rangeL, m_transformL,
                                                                blockR, m_dirR, rangeR);
  auto [boundary_cellsR, indexerR] = createAdjacentBlockIndexer(blockR, m_dirR, rangeR, m_transformR, blockL, m_dirL, rangeL);

  m_boundary_cellsL = boundary_cellsL;
  m_indexerL = indexerL;
  m_boundary_cellsR = boundary_cellsR;
  m_indexerR = indexerR;
}

std::pair<Range2D, AdjacentBlockIndexer> StructuredBlockInterface::createAdjacentBlockIndexer(
   const StructuredBlock& blockL, NeighborDirection dirL, 
   const Range& rangeL, const std::array<Int, 2>& transformL,
   const StructuredBlock& blockR, NeighborDirection dirR, const Range& rangeR)
{
  Range2D boundary_cells = getBoundaryRange(blockL, dirL, rangeL);
  std::array<UInt, 2> min_cellL{*boundary_cells.getXRange().begin(),
                                *boundary_cells.getYRange().begin()};

  UInt variable_indexR = getMinCellOnBoundary(dirL, transformL, rangeR);
  UInt constant_indexR = getConstantIndexAlongBoundary(blockR, dirR);
  std::array<UInt, 2> min_cellR;
  min_cellR[0] = to_int(dirR) % 2 == 0 ? variable_indexR : constant_indexR;
  min_cellR[1] = to_int(dirR) % 2 == 0 ? constant_indexR : variable_indexR;

  AdjacentBlockIndexer indexer(transformL, min_cellL, dirL, min_cellR);

  return {boundary_cells, indexer};
}
 
}
}