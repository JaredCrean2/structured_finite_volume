#include "gtest/gtest.h"
#include <bitset>

#include "mesh/neighbor_direction.h"
#include "utils/project_defs.h"
#include "utils/face_iterator.h"
#include "utils/cpu_timer.h"

using namespace structured_fv;

namespace {
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
  Kokkos::View<double***, HostMemorySpace> field1("field", cell_range.getXRange().size(), cell_range.getYRange().size(), 4);
  for (UInt i : cell_range.getXRange())
    for (UInt j : cell_range.getYRange())
      for (UInt v=0; v < 4; ++v)
        field1(i, j, v) = (v+1)*(i+j);

  Kokkos::View<double***, HostMemorySpace> field2("field", cell_range.getXRange().size(), cell_range.getYRange().size(), 4);
  for (UInt i : cell_range.getXRange())
    for (UInt j : cell_range.getYRange())
      for (UInt v=0; v < 4; ++v)
        field2(i, j, v) = 0;

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

__attribute((noinline)) void kernel2(const Range2D& cell_range)
{
  std::cout << "Testing with cell range " << cell_range << std::endl;
  Kokkos::View<double***, HostMemorySpace> field1("field", cell_range.getXRange().size(), cell_range.getYRange().size(), 4);
  for (UInt i : cell_range.getXRange())
    for (UInt j : cell_range.getYRange())
      for (UInt v=0; v < 4; ++v)
        field1(i, j, v) = (v+1)*(i+j);

  Kokkos::View<double***, HostMemorySpace> field2("field", cell_range.getXRange().size(), cell_range.getYRange().size(), 4);
  for (UInt i : cell_range.getXRange())
    for (UInt j : cell_range.getYRange())
      for (UInt v=0; v < 4; ++v)
        field2(i, j, v) = 0;

  {
    ScopedCPUTimer timer("iterator time");
  
    for (const FaceId& face : FaceRange(cell_range))
    {
      work(face, field1, field2);
    }
  }  

}

TEST(FaceIter, Performance)
{
  Range2D range(0, 10000, 0, 10000);
  kernel(range);
  kernel2(range);
}

TEST(FaceIter, Bits)
{
  using BitsetU = std::bitset<8*sizeof(UInt)>;
  for (UInt i=0; i < 8; ++i)
    std::cout << i << " = " << BitsetU(UInt(i)) << std::endl;

  std::cout << "equality = " << BitsetU(1==1) << std::endl;
  std::cout << "inequality = " << BitsetU(1==2) << std::endl;


}

}