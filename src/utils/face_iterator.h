#ifndef STRUCTURED_FINITE_VOLUME_UTILS_FACE_ITER_H
#define STRUCTURED_FINITE_VOLUME_UTILS_FACE_ITER_H

#include "mesh/neighbor_direction.h"
#include "range.h"
#include "bitwise.h"
#include <iosfwd>

#include <iostream>
#include <bitset>

namespace structured_fv {

using BitsetU = std::bitset<8*sizeof(UInt)>;  //TODO: DEBUGGING


struct FaceId
{
  constexpr FaceId() = default;

  constexpr FaceId(UInt cell_i_left_, UInt cell_j_left_, UInt cell_i_right_, UInt cell_j_right_,
                   mesh::NeighborDirection dirL_, mesh::NeighborDirection dirR_):
    cell_i_left(cell_i_left_),
    cell_j_left(cell_j_left_),
    cell_i_right(cell_i_right_),
    cell_j_right(cell_j_right_),
    dirL(dirL_),
    dirR(dirR_)
  {}

  UInt cell_i_left = 0;
  UInt cell_j_left = 0;
  UInt cell_i_right = 0;
  UInt cell_j_right = 0;
  mesh::NeighborDirection dirL = mesh::NeighborDirection::North;
  mesh::NeighborDirection dirR = mesh::NeighborDirection::North;
};

constexpr bool operator==(const FaceId& lhs, const FaceId& rhs)
{
  return (lhs.cell_i_left == rhs.cell_i_left) &&
         (lhs.cell_j_left == rhs.cell_j_left) &&
         (lhs.cell_i_right == rhs.cell_i_right) &&
         (lhs.cell_j_right == rhs.cell_j_right) &&
         (lhs.dirL == rhs.dirL) &&
         (lhs.dirR == rhs.dirR);
}

constexpr bool operator!=(const FaceId& lhs, const FaceId& rhs)
{
  return !(lhs == rhs);
}

std::ostream& operator<<(std::ostream& os, const FaceId& faceid);


// this works, but it is more than 10x slower than doing 2 
// separate loops over the x faces and then the y faces
class FaceIter
{
  public:

    struct PastTheEndTag {};

    FaceIter(const Range2D& cell_range) :
      m_cell_i(*cell_range.getXRange().begin()),
      m_cell_j(*cell_range.getYRange().begin()),
      m_is_x_direction(UIntMax),
      m_cell_range(cell_range),
      m_xmin(*cell_range.getXRange().begin()),
      m_ymin(*m_cell_range.getYRange().begin()),
      m_xmax_for_x_direction(*m_cell_range.getXRange().end()-2),
      m_xmax_for_y_direction(*m_cell_range.getXRange().end()-1),
      m_ymax_for_x_direction(*m_cell_range.getYRange().end()-1)
    {}


    FaceIter(const Range2D& cell_range, PastTheEndTag) :
      m_cell_i(0),
      m_cell_j(*cell_range.getYRange().end()-1),
      m_is_x_direction(0),
      m_cell_range(cell_range),
      m_xmin(*cell_range.getXRange().begin()),
      m_ymin(*m_cell_range.getYRange().begin()),
      m_xmax_for_x_direction(*m_cell_range.getXRange().end()-2),
      m_xmax_for_y_direction(*m_cell_range.getXRange().end()-1),
      m_ymax_for_x_direction(*m_cell_range.getYRange().end()-1)      
    {}

    constexpr FaceIter& operator++()
    {
      // version 3
      auto [cell_i_if_x_face, cell_j_if_x_face] = getUpdatedIndices(m_cell_i, m_cell_j, m_xmin, m_xmax_for_x_direction);
      auto [cell_i_if_y_face, cell_j_if_y_face] = getUpdatedIndices(m_cell_i, m_cell_j, m_xmin, m_xmax_for_y_direction);
      UInt at_end_of_x_faces = m_is_x_direction & mask_if_equal(m_cell_i, m_xmax_for_x_direction) & 
                                                  mask_if_equal(m_cell_j, m_ymax_for_x_direction);
      UInt x_update_mask = m_is_x_direction & ~at_end_of_x_faces;  // if m_is_x_direction && !at_end_of_x_faces
      UInt reset_to_beginning_mask = at_end_of_x_faces;
      UInt y_update_mask = ~m_is_x_direction;

      m_cell_i = (x_update_mask & cell_i_if_x_face) |
                 (reset_to_beginning_mask & m_xmin) |
                 (y_update_mask & cell_i_if_y_face);

      m_cell_j = (x_update_mask & cell_j_if_x_face) |
                 (reset_to_beginning_mask & m_ymin) |
                 (y_update_mask & cell_j_if_y_face);
      
      m_is_x_direction = m_is_x_direction & ~at_end_of_x_faces;

      return *this;
    }

    constexpr FaceId operator*() const
    {     
      if (m_is_x_direction)
      {
        return {m_cell_i, m_cell_j, m_cell_i+1, m_cell_j,
                mesh::NeighborDirection::East, mesh::NeighborDirection::West};
      } else
      {
        return {m_cell_i, m_cell_j, m_cell_i, m_cell_j+1,
                mesh::NeighborDirection::North, mesh::NeighborDirection::South};
      }    
    }  

  private:
    static constexpr std::pair<UInt, UInt> getUpdatedIndices(UInt cell_i, UInt cell_j, UInt xmin, UInt xmax)
    {
      UInt delta_y = cell_i == xmax;

      UInt new_i_if_at_xmax = xmin;
      UInt new_i_if_not_at_xmax = cell_i+1;
      UInt all_ones_if_at_xmax = mask_if_equal(cell_i, xmax);
      UInt new_i = (new_i_if_at_xmax & all_ones_if_at_xmax) + (new_i_if_not_at_xmax & ~all_ones_if_at_xmax);

      return {new_i, cell_j + delta_y};
    }

    // index of left/bottom cell of the interface
    UInt m_cell_i;
    UInt m_cell_j;
    UInt m_is_x_direction;
    Range2D m_cell_range;

    UInt m_xmin;
    UInt m_ymin;
    UInt m_xmax_for_x_direction;
    UInt m_xmax_for_y_direction;
    UInt m_ymax_for_x_direction;

    friend constexpr bool operator==(const FaceIter& lhs, const FaceIter& rhs);
};

constexpr bool operator==(const FaceIter& lhs, const FaceIter& rhs)
{
  return //lhs.m_cell_range     == rhs.m_cell_range &&
         lhs.m_cell_i         == rhs.m_cell_i &&
         lhs.m_cell_j         == rhs.m_cell_j &&
         lhs.m_is_x_direction == rhs.m_is_x_direction;
}

constexpr bool operator!=(const FaceIter& lhs, const FaceIter& rhs)
{
  return !(lhs == rhs);
}

// iterates over all the internal faces of the given cell range.
// internal faces means all faces that have a cell within the range
// on both sides
class FaceRange
{
  public:
    FaceRange(const Range2D& cell_range) :
      m_cell_range(cell_range)
    {}

    FaceIter begin() const { return FaceIter(m_cell_range); }

    FaceIter end() const { return FaceIter(m_cell_range, FaceIter::PastTheEndTag{}); }

  private:
    Range2D m_cell_range;
};

}

#endif