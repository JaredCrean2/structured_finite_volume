#include "gtest/gtest.h"
#include <bitset>

#include "mesh/neighbor_direction.h"
#include "utils/face_iter_per_direction.h"
#include "utils/project_defs.h"
#include "utils/face_iterator.h"
#include "utils/cpu_timer.h"

using namespace structured_fv;

namespace {
using HostView = Kokkos::View<double***, HostMemorySpace>;

HostView create_field1(const Range2D& cell_range)
{
  HostView field1("field", cell_range.getXRange().size(), cell_range.getYRange().size(), 4);
  for (UInt i : cell_range.getXRange())
    for (UInt j : cell_range.getYRange())
      for (UInt v=0; v < 4; ++v)
        field1(i, j, v) = (v+1)*(i+j);

  return field1;
}

HostView create_field2(const Range2D& cell_range)
{
  HostView field2("field", cell_range.getXRange().size(), cell_range.getYRange().size(), 4);
  for (UInt i : cell_range.getXRange())
    for (UInt j : cell_range.getYRange())
      for (UInt v=0; v < 4; ++v)
        field2(i, j, v) = 0;

  return field2 ; 
}

void work(const FaceId& face, Kokkos::View<double***, HostMemorySpace>& field1,
                         Kokkos::View<double***, HostMemorySpace>& field2)
{
  for (UInt i=0; i < 4; ++i)
  {
    field2(face.cell_i_left, face.cell_i_right, i) = field1(face.cell_i_left, face.cell_j_left, i) + 
                                                     field1(face.cell_i_left, face.cell_j_right, i);
  }
}

__attribute((noinline)) void kernel(const Range2D& cell_range)
{
  std::cout << "Testing with cell range " << cell_range << std::endl;
  HostView field1 = create_field1(cell_range);
  HostView field2 = create_field2(cell_range);

  {
    ScopedCPUTimer timer("explicit loop time");

    for (UInt i : Range(*cell_range.getXRange().begin(), *cell_range.getXRange().end()-1))
      for (UInt j : cell_range.getYRange())
      {
        FaceId face{i, j, i+1, j, mesh::NeighborDirection::East, mesh::NeighborDirection::West};
        work(face, field1, field2);
      }

    for (UInt i : cell_range.getXRange())
      for (UInt j : Range(*cell_range.getYRange().begin(), *cell_range.getYRange().end()-1))
      {
        FaceId face{i, j, i, j+1, mesh::NeighborDirection::North, mesh::NeighborDirection::South};
        work(face, field1, field2);
      }
  }

}

/*
__attribute((noinline)) void kernel2(const Range2D& cell_range)
{
  std::cout << "Testing with cell range " << cell_range << std::endl;
  HostView field1 = create_field1(cell_range);
  HostView field2 = create_field2(cell_range);

  {
    ScopedCPUTimer timer("iterator time");
  
    for (const FaceId& face : FaceRange(cell_range))
    {
      work(face, field1, field2);
    }
  }  
}
*/


__attribute((noinline)) void kernel4(const Range2D& cell_range)
{
  std::cout << "Testing with cell range " << cell_range << std::endl;
  HostView field1 = create_field1(cell_range);
  HostView field2 = create_field2(cell_range);

  {
    ScopedCPUTimer timer("2 iterator single class loop time");
    FaceRangePerDirection faces(cell_range);

    for (UInt i : faces.getXRange(XDirTag()))
      for (UInt j : faces.getYRange(XDirTag()))
      {
        FaceId face = faces.getFaceId(XDirTag(), i, j);
        work(face, field1, field2);
      }

    for (UInt i : faces.getXRange(YDirTag()))
      for (UInt j : faces.getYRange(YDirTag()))
      {
        FaceId face = faces.getFaceId(YDirTag(), i, j);
        work(face, field1, field2);
      }
  }
}

TEST(FaceIter, Performance)
{
  Range2D range(0, 10000, 0, 10000);
  kernel(range);
  //kernel2(range);
  kernel4(range);
}

}