#ifndef STRUCTURED_FINITE_VOLUME_DISC_FACE_FIELD_H
#define STRUCTURED_FINITE_VOLUME_DISC_FACE_FIELD_H

#include "utils/neighbor_direction.h"
#include "utils/face_iter_per_direction.h"
#include "utils/project_defs.h"
#include "utils/traits.h"
#include "disc/discretization.h"
#include "disc/vert_field.h"

namespace structured_fv {
namespace disc {

template <typename T>
class FaceField
{
  public:
    using FieldData = Kokkos::View<T***, HostMemorySpace>;
    using ConstFieldData = Kokkos::View<const T***, HostMemorySpace>;

    FaceField(const StructuredDisc& disc, Int nvals_per_element) :
      m_disc(disc),
      m_nvals_per_element(nvals_per_element)
    {
      for (UInt i=0; i < disc.getNumBlocks(); ++i)
      {
        const StructuredBlock& block = disc.getBlock(i);
        FaceRangePerDirection faces = block.getOwnedAndGhostFacesWithCorners();

        m_data[getDirection(NeighborDirection::East)].emplace_back("field_data",
          faces.getXRange(XDirTag()).size(), faces.getYRange(XDirTag()).size(), nvals_per_element);

        m_data[getDirection(NeighborDirection::North)].emplace_back("field_data",
          faces.getXRange(YDirTag()).size(), faces.getYRange(YDirTag()).size(), nvals_per_element);
      }
    }

    FaceField(const FaceField& other) = delete;

    FaceField& operator=(const FaceField& other) = delete;

    UInt getNumBlocks() const { return m_data.size(); }

    // if dir == North or South, returns the Y direction faces (parallel to x axis)
    // if dir == East or West, returns the X direction faces (parallel to y axis)
    FieldData& getData(UInt block, NeighborDirection dir)
    {
      return m_data[getDirection(dir)][block];
    }

    const ConstFieldData& getData(UInt block, NeighborDirection dir) const
    {
      return m_data[getDirection(dir)][block];
    }

    T& operator()(UInt block, NeighborDirection dir, UInt i, UInt j, UInt v)
    {

      return m_data[getDirection(dir)][block](i, j, v);
    }      

    const T& operator()(UInt block, NeighborDirection dir, UInt i, UInt j, UInt v) const
    { 
      return m_data[getDirection(dir)][block](i, j, v);
    }

    void set(const T& val);

    // Func is a callable object (Real x, Real y) -> FixedVec<T, num_vals_per_element>
    // return type can be anything of the correct length that supports operator[]
    // x and y are the centroid of each face
    template <typename Func, IsFuncXY_t<Func> = true>
    void set(Func func, bool include_corners=false);

  private:
    static constexpr UInt getDirection(NeighborDirection dir)
    {
      return to_int(dir) % 2;
    }

    const StructuredDisc& m_disc;
    const UInt m_nvals_per_element;
    FixedVec<std::vector<FieldData>, 2> m_data;
};

template <typename T>
void FaceField<T>::set(const T& val)
{
  for (UInt dir=0; dir < 2; ++dir)
    for (FieldData& field_data : m_data[dir])
      Kokkos::deep_copy(field_data, val);
}


template <typename T>
template <typename Func, IsFuncXY_t<Func>>
void FaceField<T>::set(Func func, bool include_corners)
{
  for (UInt block_id=0; block_id < getNumBlocks(); ++block_id)
  {
    const StructuredBlock& block = m_disc.getBlock(block_id);
    auto coords = m_disc.getCoordField()->getData(block_id);
    auto east_data = getData(block_id, NeighborDirection::East);
    auto north_data = getData(block_id, NeighborDirection::North);

    auto set_values = [&](const FaceRangePerDirection& x_faces, const FaceRangePerDirection& y_faces)
    {
      for (UInt i : x_faces.getXRange(XDirTag()))
        for (UInt j : x_faces.getYRange(XDirTag()))
        {
          double x = (coords(i, j, 0) + coords(i, j+1, 0))/2;
          double y = (coords(i, j, 1) + coords(1, j+1, 1))/2;
          auto vals = func(x, y);
          for (UInt v=0; v < m_nvals_per_element; ++v)
            east_data(i, j, v) = vals[v];
        }

      for (UInt i : y_faces.getXRange(YDirTag()))
        for (UInt j : y_faces.getYRange(YDirTag()))
        {
          double x = (coords(i, j, 0) + coords(i+1, j, 0))/2;
          double y = (coords(i, j, 1) + coords(i+1, j, 1))/2;
          auto vals = func(x, y);
          for (UInt v=0; v < m_nvals_per_element; ++v)
            north_data(i, j, v) = vals[v];
        }
    };

    if (include_corners)
    {
      FaceRangePerDirection all_faces = block.getOwnedAndGhostFacesWithCorners();
      set_values(all_faces, all_faces);
    } else
    {
      FaceRangePerDirection x_faces = block.getOwnedAndGhostXFaces();
      FaceRangePerDirection y_faces = block.getOwnedAndGhostYFaces(); 
      set_values(x_faces, y_faces);
    }
  }
}

}
}

#endif