#ifndef STRUCTURED_FINITE_VOLUME_UTILS_FACE_ITER_H
#define STRUCTURED_FINITE_VOLUME_UTILS_FACE_ITER_H

#include "mesh/neighbor_direction.h"
#include "range.h"
#include <iosfwd>

namespace structured_fv {

struct FaceId
{
  UInt cell_i_left;
  UInt cell_j_left;
  UInt cell_i_right;
  UInt cell_j_right;
  mesh::NeighborDirection dirL;
  mesh::NeighborDirection dirR;
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


class FaceIter
{
  public:

    struct PastTheEndTag {};

    FaceIter(const Range2D& cell_range) :
      m_cell_range(cell_range),
      m_cell_i(*cell_range.getXRange().begin()),
      m_cell_j(*cell_range.getYRange().begin()),
      m_is_x_direction(true)
    {}


    FaceIter(const Range2D& cell_range, PastTheEndTag) :
      m_cell_range(cell_range),
      m_cell_i(0),
      m_cell_j(*cell_range.getYRange().end()-1),
      m_is_x_direction(false)
    {}

    constexpr FaceIter& operator++()
    {      
      ++m_cell_i;
      if (m_is_x_direction)
      {
        if (m_cell_i == (*m_cell_range.getXRange().end()-1) && 
            m_cell_j == (*m_cell_range.getYRange().end()-1))
        {
          m_is_x_direction = false;
          m_cell_i = *m_cell_range.getXRange().begin();
          m_cell_j = *m_cell_range.getYRange().begin();
        } else if (m_cell_i == (*m_cell_range.getXRange().end()-1))
        {
          m_cell_i = *m_cell_range.getXRange().begin();
          ++m_cell_j;
        }
      } else
      {
        if (m_cell_i == (*m_cell_range.getXRange().end()))
        {
          m_cell_i = *m_cell_range.getXRange().begin();
          ++m_cell_j;
        }

#ifndef NDEBUG
        if (m_cell_j >= (*m_cell_range.getYRange().end()-1) && 
            m_cell_i > *m_cell_range.getXRange().begin())
        {
          throw std::runtime_error("attempted to increment past-the-end FaceIter");
        }
#endif
      }
      return *this;
    }

    FaceId operator*() const
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
    Range2D m_cell_range;
    // index of left/bottom cell of the interface
    UInt m_cell_i;
    UInt m_cell_j;
    bool m_is_x_direction;

    friend constexpr bool operator==(const FaceIter& lhs, const FaceIter& rhs);
};

constexpr bool operator==(const FaceIter& lhs, const FaceIter& rhs)
{
  return lhs.m_cell_range     == rhs.m_cell_range &&
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