#include "neighbor_direction.h"

namespace structured_fv {
namespace mesh {

UInt getConstantIndexAlongBoundary(const Range2D& block_range, NeighborDirection dir)
{
  switch (dir)
  {
    case NeighborDirection::North: { return *(block_range.getYRange().end())-1; }
    case NeighborDirection::East:  { return *(block_range.getXRange().end())-1; }
    case NeighborDirection::South: { return *(block_range.getXRange().begin()); }
    case NeighborDirection::West:  { return *(block_range.getYRange().begin()); }
  }
}

Range getVariableIndexAlongBoundary(const Range2D& block_range, NeighborDirection dir)
{
  if (static_cast<int>(dir) % 2 == 0)
  {
    return block_range.getXRange();
  } else
  {
    return block_range.getYRange();
  }
}

Range2D getBoundaryRange(const Range2D& block_range, NeighborDirection dir, const Range& boundary_subset)
{
  UInt constant_index = getConstantIndexAlongBoundary(block_range, dir);
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