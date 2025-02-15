#include "neighbor_direction.h"
#include "structured_block.h"

namespace structured_fv {
namespace mesh {

UInt getConstantIndexAlongBoundary(const StructuredBlock& block, NeighborDirection dir)
{
  switch (dir)
  {
    case NeighborDirection::North: { return *(block.getOwnedCells().getYRange().end())-1; }
    case NeighborDirection::East:  { return *(block.getOwnedCells().getXRange().end())-1; }
    case NeighborDirection::South: { return *(block.getOwnedCells().getXRange().begin()); }
    case NeighborDirection::West:  { return *(block.getOwnedCells().getYRange().begin()); }
  }
}

Range getVariableIndexAlongBoundary(const StructuredBlock& block, NeighborDirection dir)
{
  if (static_cast<int>(dir) % 2 == 0)
  {
    return block.getOwnedCells().getXRange();
  } else
  {
    return block.getOwnedCells().getYRange();
  }
}

Range2D getBoundaryRange(const StructuredBlock& block, NeighborDirection dir, const Range& boundary_subset)
{
  UInt constant_index = getConstantIndexAlongBoundary(block, dir);
  if (to_int(dir) % 2 == 0)
  {
    return Range2D(*boundary_subset.begin(), *boundary_subset.end(), constant_index, constant_index+1);
  } else
  {
    return Range2D(constant_index, constant_index+1, *boundary_subset.begin(), *boundary_subset.end());
  }
}
}
}