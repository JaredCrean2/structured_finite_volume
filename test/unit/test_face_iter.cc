#include "gtest/gtest.h"
#include "utils/face_iterator.h"

using namespace structured_fv;

namespace {

FaceId makeFaceId(UInt cell_i_left, UInt cell_j_left, UInt cell_i_right, UInt cell_j_right,
                  NeighborDirection dirL, NeighborDirection dirR)
{
  return FaceId{cell_i_left, cell_j_left, cell_i_right, cell_j_right, dirL, dirR};
}

}

TEST(FaceIter, FaceIdEquality)
{
  FaceId face1{0, 1, 2, 3, NeighborDirection::North, NeighborDirection::South};
  FaceId face2{0, 1, 2, 3, NeighborDirection::North, NeighborDirection::South};
  FaceId face3{0, 2, 2, 3, NeighborDirection::North, NeighborDirection::South};

  EXPECT_TRUE(face1 == face2);
  EXPECT_FALSE(face1 != face2);
  EXPECT_FALSE(face1 == face3);
  EXPECT_TRUE(face1 != face3);
}


TEST(FaceIter, Iteration)
{
  FaceIter iter(Range2D(0, 3, 0, 4));

  // East-West
  EXPECT_EQ(*iter, makeFaceId(0, 0, 1, 0, NeighborDirection::East, NeighborDirection::West));
  ++iter;
  EXPECT_EQ(*iter, makeFaceId(1, 0, 2, 0, NeighborDirection::East, NeighborDirection::West));

  ++iter;
  EXPECT_EQ(*iter, makeFaceId(0, 1, 1, 1, NeighborDirection::East, NeighborDirection::West));
  ++iter;
  EXPECT_EQ(*iter, makeFaceId(1, 1, 2, 1, NeighborDirection::East, NeighborDirection::West));

  ++iter;
  EXPECT_EQ(*iter, makeFaceId(0, 2, 1, 2, NeighborDirection::East, NeighborDirection::West));
  ++iter;
  EXPECT_EQ(*iter, makeFaceId(1, 2, 2, 2, NeighborDirection::East, NeighborDirection::West));

  ++iter;
  EXPECT_EQ(*iter, makeFaceId(0, 3, 1, 3, NeighborDirection::East, NeighborDirection::West));
  ++iter;
  EXPECT_EQ(*iter, makeFaceId(1, 3, 2, 3, NeighborDirection::East, NeighborDirection::West));

  // North-South
  ++iter;
  EXPECT_EQ(*iter, makeFaceId(0, 0, 0, 1, NeighborDirection::North, NeighborDirection::South));
  ++iter;
  EXPECT_EQ(*iter, makeFaceId(1, 0, 1, 1, NeighborDirection::North, NeighborDirection::South));
  ++iter;
  EXPECT_EQ(*iter, makeFaceId(2, 0, 2, 1, NeighborDirection::North, NeighborDirection::South));

  ++iter;
  EXPECT_EQ(*iter, makeFaceId(0, 1, 0, 2, NeighborDirection::North, NeighborDirection::South));
  ++iter;
  EXPECT_EQ(*iter, makeFaceId(1, 1, 1, 2, NeighborDirection::North, NeighborDirection::South));
  ++iter;
  EXPECT_EQ(*iter, makeFaceId(2, 1, 2, 2, NeighborDirection::North, NeighborDirection::South));  

  ++iter;
  EXPECT_EQ(*iter, makeFaceId(0, 2, 0, 3, NeighborDirection::North, NeighborDirection::South));
  ++iter;
  EXPECT_EQ(*iter, makeFaceId(1, 2, 1, 3, NeighborDirection::North, NeighborDirection::South));
  ++iter;
  EXPECT_EQ(*iter, makeFaceId(2, 2, 2, 3, NeighborDirection::North, NeighborDirection::South));  
}

TEST(FaceIter, Range)
{
  FaceRange range(Range2D(0, 3, 0, 4));
  EXPECT_EQ(*(range.end()), makeFaceId(0, 3, 0, 4, NeighborDirection::North, NeighborDirection::South));

  int nfaces = 0;
  for (const FaceId& face : range)
  {
    if (nfaces < 8)
    {
      EXPECT_EQ(face.dirL, NeighborDirection::East);
    } else
    {
      EXPECT_EQ(face.dirL, NeighborDirection::North);
    }
    nfaces++;
  }

  EXPECT_EQ(nfaces, 2*4 + 3*3);
}
