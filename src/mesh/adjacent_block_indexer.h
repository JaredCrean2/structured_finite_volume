#ifndef STRUCTURED_FINITE_VOLUME_MESH_ADJACENT_BLOCK_INDEXER_H
#define STRUCTURED_FINITE_VOLUME_MESH_ADJACENT_BLOCK_INDEXER_H

#include <cmath>
#include "utils/project_defs.h"
#include "structured_block.h"

namespace structured_fv {
namespace mesh {


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

inline std::ostream& operator<<(std::ostream& os, NeighborDirection dir)
{
  std::array<std::string, 4> names{"North", "East", "South", "West"};
  os << names[to_int(dir)];
  return os;
}

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

constexpr std::array<Int, 2> computeIndices(NeighborDirection dir, UInt offset, Int i, Int j)
{
  std::array<Int, 2> idxs = {i, j};
  //Int sign = to_int(dir) < 2 ? 1 : -1;
  Int sign = 1 - 2*(to_int(dir)/2);
  idxs[0] += sign * (to_int(dir) % 2) * offset;
  idxs[1] += sign * ((to_int(dir) + 1) % 2) * offset;

  return idxs;
}

inline UInt getConstantIndexAlongBoundary(const StructuredBlock& block, NeighborDirection dir)
{
  switch (dir)
  {
    case NeighborDirection::North: { return *(block.getOwnedCells().getYRange().end())-1; }
    case NeighborDirection::East:  { return *(block.getOwnedCells().getXRange().end())-1; }
    case NeighborDirection::South: { return *(block.getOwnedCells().getXRange().begin()); }
    case NeighborDirection::West:  { return *(block.getOwnedCells().getYRange().begin()); }
  }
}

inline Range getVariableIndexAlongBoundary(const StructuredBlock& block, NeighborDirection dir)
{
  if (static_cast<int>(dir) % 2 == 0)
  {
    return block.getOwnedCells().getXRange();
  } else
  {
    return block.getOwnedCells().getYRange();
  }
}

inline Range2D getBoundaryRange(const StructuredBlock& block, NeighborDirection dir, const Range& boundary_subset)
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

inline std::array<Int, 2> getTransform(Int rotation_left, Int rotation_right)
{
  Int delta_r = (rotation_right - rotation_left + 4) % 4;
  std::array<std::array<Int, 2>, 4> transforms;
  transforms[0] = {1, 2};
  transforms[1] = {-2, 1};
  transforms[2] = {-1, -2};
  transforms[3] = {2, -1};

  return transforms[delta_r];
}

inline std::array<Int, 2> getInverseTransform(const std::array<Int, 2>& transform)
{
  std::array<Int, 2> inverse_transform;
  for (int i=0; i < 2; ++i)
    inverse_transform[std::abs(transform[i])-1] = sgn(transform[i])*(i+1);

  return inverse_transform;
}

// given an direction and a transfrom from the left block perspective, return
// the entry of rangeR that corresponds to the first cell on the boundary in blockL
inline UInt getMinCellOnBoundary(NeighborDirection dirL, const std::array<Int, 2>& transformL, const Range& rangeR)
{
  bool is_reversed = (to_int(dirL) % 2 == 0 && sgn(transformL[1]) < 0) ||
                     (to_int(dirL) % 2 == 1 && sgn(transformL[0]) < 0);
  return is_reversed ? rangeR(rangeR.size()-1) : rangeR(0);                  
}


class AdjacentBlockIndexer
{
  public:
    // transform is the transformation of the coordinate system from the
    // left block to the right block, according to the CGNS definition:
    // transform[i] gives the image of the i'th axis of the left block in
    // the right block.  For example, if the +x direction of the left block
    // is the -y direction of the right block, then transform[0] = -2
    // Note: this is a little weird because CGNS uses 1 based indexing, but
    // std::array uses zero based.
    // left_block_min_cell is the minimum cell on the left block that is adjacent
    // to a cell on the right block
    // right_block_min_cell is the index of the adjacent cell in the right block
    AdjacentBlockIndexer(const std::array<Int, 2>& transform,
                         const std::array<UInt, 2>& left_block_min_cell,
                         NeighborDirection direction,
                         const std::array<UInt, 2>& right_block_min_cell) :
      m_transform_matrix("transform_matrix"),
      m_left_block_min_cell(left_block_min_cell),
      m_direction(direction),
      m_right_block_min_cell(right_block_min_cell)
    {
      assert(transform[0] != 0 && transform[1] != 0);

      m_transform_matrix(0, 0) = sgn(transform[0])*del(transform[0], 1);
      m_transform_matrix(0, 1) = sgn(transform[1])*del(transform[1], 1);
      m_transform_matrix(1, 0) = sgn(transform[0])*del(transform[0], 2);
      m_transform_matrix(1, 1) = sgn(transform[1])*del(transform[1], 2);
    }

    AdjacentBlockIndexer() :
      AdjacentBlockIndexer({1, 2}, {0, 0}, NeighborDirection::East, {0, 0})
    {}

    // i and j are used to index into the adjacent block as though it
    // were a continuation of the current block.  For example, if
    // the current block is connected to the adjacent block to the East
    // then (nx, 0) will give the indices of the cell in the adjacent
    // block that is next to the bottom right cell in the current block.
    // (where nx is the number of cells in the x direction of the current block)
    // If the blocks are joined to the West, then (-1, 0) gives the indices of
    // the cell in the adjacent block next to the lower left cell of the current
    // block
    std::array<UInt, 2> operator()(Int i, Int j) const
    {
      std::array<Int, 4> ioffset{0,  -1, 0, 1};
      std::array<Int, 4> joffset{-1,  0, 1, 0};

      i += ioffset[static_cast<int>(m_direction)];
      j += joffset[static_cast<int>(m_direction)];

      Int iprime = m_transform_matrix(0, 0) * (i - m_left_block_min_cell[0]) +
                   m_transform_matrix(0, 1) * (j - m_left_block_min_cell[1]);
      Int jprime = m_transform_matrix(1, 0) * (i - m_left_block_min_cell[0]) +
                   m_transform_matrix(1, 1) * (j - m_left_block_min_cell[1]);

      return {iprime + m_right_block_min_cell[0], jprime + m_right_block_min_cell[1]};
    }

    std::array<UInt, 2> operator()(const std::array<UInt, 2>& ij) const
    {
      return this->operator()(ij[0], ij[1]);
    }

    std::array<UInt, 2> operator()(const std::array<Int, 2>& ij) const
    {
      return this->operator()(ij[0], ij[1]);
    }

    const std::array<UInt, 2>& getRightBlockMinCell() const { return m_right_block_min_cell; }

  private:
    // TODO: put this on stack?
    Kokkos::View<Int[2][2], HostMemorySpace> m_transform_matrix; 
    std::array<UInt, 2> m_left_block_min_cell;
    NeighborDirection m_direction;
    std::array<UInt, 2> m_right_block_min_cell;
};

}
}
#endif