#ifndef STRUCTURED_FINITE_VOLUME_NEIGHBOR_DIRECTION_H
#define STRUCTURED_FINITE_VOLUME_NEIGHBOR_DIRECTION_H

#include "utils/project_defs.h"
#include "utils/range.h"
#include <iosfwd>

namespace structured_fv {


template <typename T>
constexpr Int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

template <typename T>
constexpr Int del(T val1, T val2)
{
  return std::abs(val1) == std::abs(val2) ? 1 : 0;
}

enum class NeighborDirection
{
  North,
  East,
  South,
  West
};

constexpr int to_int(NeighborDirection dir)
{
  return static_cast<int>(dir);
}

std::ostream& operator<<(std::ostream& os, NeighborDirection dir);

// given a direction for a block with rotation=0, computes the direction
// of the same side (in the global coordinate system) of the block
// after a given rotation.
// Ex. Given the East side of a block and a rotation 1,
//     this function returns South because the block was rotated
//     90 degrees counter clockwise, so the rightmost side of the block
//     is now the South side
constexpr NeighborDirection rotate(NeighborDirection dir, UInt rotation)
{
  return static_cast<NeighborDirection>((static_cast<int>(dir) + rotation) % 4);
}

// gives a 1-based axis and a sign that describes which face of
// the rectangle a given NeighborDirection is.
// Ex. East is the positive end of the x axis -> +1
constexpr Int toSignedAxis(NeighborDirection dir)
{
  constexpr std::array<Int, 4> dirs{2, 1, -2, -1};
  return dirs[to_int(dir)];
}

// given the result of toSignedAxis, compute the NeighborDirection
constexpr NeighborDirection toNeighborDirection(Int signed_axis)
{
  constexpr std::array<NeighborDirection, 5> dirs{NeighborDirection::South,
                                                  NeighborDirection::West,
                                                  NeighborDirection::South, // unused
                                                  NeighborDirection::East,
                                                  NeighborDirection::North};

  return dirs[signed_axis + 2];
}

constexpr NeighborDirection getNeighborImage(NeighborDirection dir, const std::array<Int, 2>& transform)
{
  Int signed_axis = toSignedAxis(dir);
  Int other_block_axis = -sgn(signed_axis) *transform[std::abs(signed_axis)-1];
  return toNeighborDirection(other_block_axis);
}

constexpr std::array<Int, 2> computeIndices(NeighborDirection dir, Int offset, Int i, Int j)
{
  std::array<Int, 2> idxs = {i, j};
  //Int sign = to_int(dir) < 2 ? 1 : -1;
  Int sign = 1 - 2*(to_int(dir)/2);
  idxs[0] += sign * (to_int(dir) % 2) * offset;
  idxs[1] += sign * ((to_int(dir) + 1) % 2) * offset;

  return idxs;
}

UInt getConstantIndexAlongBoundary(const Range2D& block_range, NeighborDirection dir);

Range getVariableIndexAlongBoundary(const Range2D& block_range, NeighborDirection dir);

Range2D getBoundaryRange(const Range2D& block_range, NeighborDirection dir, const Range& boundary_subset);

std::array<Int, 2> getTransform(Int rotation_left, Int rotation_right);

std::array<Int, 2> getInverseTransform(const std::array<Int, 2>& transform);

// given an direction and a transfrom from the left block perspective, return
// the entry of rangeR that corresponds to the first cell on the boundary in blockL
UInt getMinEntityOnBoundary(NeighborDirection dirL, const std::array<Int, 2>& transformL, const Range& rangeR);

}

#endif