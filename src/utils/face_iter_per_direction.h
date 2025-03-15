#ifndef STRUCTURED_FINITE_VOLUME_UTILS_FACE_ITER_PER_DIRECTION_H
#define STRUCTURED_FINITE_VOLUME_UTILS_FACE_ITER_PER_DIRECTION_H

#include "range.h"
#include "face_iterator.h"

namespace structured_fv {

struct XDirTag {};
struct YDirTag {};

class FaceRangePerDirection
{
  public:
    // iterates over the faces of the given cell range.
    // if include_boundary is false, only includes interior faces
    // if include_boundary is true, include all faces
    // This requires there to be cells to the left and below the
    // bottom left cell of the cell_range
    FaceRangePerDirection(const Range2D& cell_range, bool include_boundary=false)
    {
      if (include_boundary)
      {
        if (*cell_range.getXRange().begin() == 0 || *cell_range.getYRange().begin() == 0)
          throw std::runtime_error("cannot iterate over cell faces with boundaries if the first cell has index 0");
        
        assert(*cell_range.getYRange().begin() > 0);
        m_face_range_x = Range2D(*cell_range.getXRange().begin()-1, *cell_range.getXRange().end(),
                                 *cell_range.getYRange().begin(),   *cell_range.getYRange().end()),
        m_face_range_y = Range2D(*cell_range.getXRange().begin(),   *cell_range.getXRange().end(),
                                 *cell_range.getYRange().begin()-1, *cell_range.getYRange().end());
      } else
      {
        m_face_range_x = Range2D(*cell_range.getXRange().begin(), *cell_range.getXRange().end()-1,
                                 *cell_range.getYRange().begin(), *cell_range.getYRange().end()),
        m_face_range_y = Range2D(*cell_range.getXRange().begin(), *cell_range.getXRange().end(),
                                 *cell_range.getYRange().begin(), *cell_range.getYRange().end()-1);
      }
    }

    Range getXRange(XDirTag) const { return m_face_range_x.getXRange(); }

    Range getYRange(XDirTag) const { return m_face_range_x.getYRange(); }

    FaceId getFaceId(XDirTag, UInt cell_i, UInt cell_j) const
    {
      return {cell_i, cell_j, cell_i+1, cell_j,
              mesh::NeighborDirection::East, mesh::NeighborDirection::West};
    }

    Range getXRange(YDirTag) const { return m_face_range_y.getXRange(); }

    Range getYRange(YDirTag) const { return m_face_range_y.getYRange(); }

    FaceId getFaceId(YDirTag, UInt cell_i, UInt cell_j) const
    {
      return {cell_i, cell_j, cell_i, cell_j+1,
              mesh::NeighborDirection::North, mesh::NeighborDirection::South};
    }    

  private:
    Range2D m_face_range_x;
    Range2D m_face_range_y;
};


}

#endif