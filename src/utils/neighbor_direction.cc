#include "neighbor_direction.h"
#include <iostream>

namespace structured_fv {


std::ostream& operator<<(std::ostream& os, NeighborDirection dir)
{
  std::array<std::string, 4> names{"North", "East", "South", "West"};
  os << names[to_int(dir)];
  return os;
}

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

std::array<Int, 2> getTransform(Int rotation_left, Int rotation_right)
{
  Int delta_r = (rotation_right - rotation_left + 4) % 4;
  std::array<std::array<Int, 2>, 4> transforms;
  transforms[0] = {1, 2};
  transforms[1] = {-2, 1};
  transforms[2] = {-1, -2};
  transforms[3] = {2, -1};

  return transforms[delta_r];
}

std::array<Int, 2> getInverseTransform(const std::array<Int, 2>& transform)
{
  std::array<Int, 2> inverse_transform;
  for (int i=0; i < 2; ++i)
    inverse_transform[std::abs(transform[i])-1] = sgn(transform[i])*(i+1);

  return inverse_transform;
}

// given an direction and a transfrom from the left block perspective, return
// the entry of rangeR that corresponds to the first cell on the boundary in blockL
UInt getMinEntityOnBoundary(NeighborDirection dirL, const std::array<Int, 2>& transformL, const Range& rangeR)
{
  bool is_reversed = (to_int(dirL) % 2 == 0 && sgn(transformL[0]) < 0) ||
                     (to_int(dirL) % 2 == 1 && sgn(transformL[1]) < 0);
  return is_reversed ? rangeR(rangeR.size()-1) : rangeR(0);                  
}

}