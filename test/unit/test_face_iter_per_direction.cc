#include "gtest/gtest.h"
#include "utils/neighbor_direction.h"
#include "utils/face_iter_per_direction.h"
#include "utils/face_iterator.h"

using namespace structured_fv;

namespace {

FaceId makeFaceId(UInt cell_i_left, UInt cell_j_left, UInt cell_i_right, UInt cell_j_right,
                  NeighborDirection dirL, NeighborDirection dirR)
{
  return FaceId{cell_i_left, cell_j_left, cell_i_right, cell_j_right, dirL, dirR};
}

}

TEST(FaceIterPerDirection, XDirection)
{
  Range2D cell_range(0, 3, 0, 4);
  FaceRangePerDirection face_range(cell_range);
  EXPECT_EQ(face_range.getXRange(XDirTag()), Range(1, 3));
  EXPECT_EQ(face_range.getYRange(XDirTag()), Range(0, 4));

  EXPECT_EQ(face_range.getFaceId(XDirTag(), 1, 0), makeFaceId(0, 0, 1, 0, NeighborDirection::East, NeighborDirection::West));
  EXPECT_EQ(face_range.getFaceId(XDirTag(), 2, 1), makeFaceId(1, 1, 2, 1, NeighborDirection::East, NeighborDirection::West));
}

TEST(FaceIterPerDirection, XDirectionFaceCount)
{
  Range2D cell_range(0, 3, 0, 4);
  FaceRangePerDirection face_range(cell_range);

  int nfaces = 0;
  for ([[maybe_unused]] UInt i : face_range.getXRange(XDirTag()))
    for ([[maybe_unused]] UInt j : face_range.getYRange(XDirTag()))
      nfaces++;

  EXPECT_EQ(nfaces, 8);
}

TEST(FaceIterPerDirection, XDirectionWithGhosts)
{
  Range2D cell_range(0, 3, 0, 4);
  FaceRangePerDirection face_range(cell_range, true);
  EXPECT_EQ(face_range.getXRange(XDirTag()), Range(0, 4));
  EXPECT_EQ(face_range.getYRange(XDirTag()), Range(0, 4));

  EXPECT_EQ(face_range.getFaceId(XDirTag(), 1, 0), makeFaceId(0, 0, 1, 0, NeighborDirection::East, NeighborDirection::West));
  EXPECT_EQ(face_range.getFaceId(XDirTag(), 1, 1), makeFaceId(0, 1, 1, 1, NeighborDirection::East, NeighborDirection::West));
}

TEST(FaceIterPerDirection, XDirectionWithGhostsFaceCount)
{
  Range2D cell_range(0, 3, 0, 4);
  FaceRangePerDirection face_range(cell_range, true);

  int nfaces = 0;
  for ([[maybe_unused]] UInt i : face_range.getXRange(XDirTag()))
    for ([[maybe_unused]] UInt j : face_range.getYRange(XDirTag()))
      nfaces++;

  EXPECT_EQ(nfaces, 16);
}

TEST(FaceIterPerDirection, YDirection)
{
  Range2D cell_range(0, 3, 0, 4);
  FaceRangePerDirection face_range(cell_range);
  EXPECT_EQ(face_range.getXRange(YDirTag()), Range(0, 3));
  EXPECT_EQ(face_range.getYRange(YDirTag()), Range(1, 4));

  EXPECT_EQ(face_range.getFaceId(YDirTag(), 0, 1), makeFaceId(0, 0, 0, 1, NeighborDirection::North, NeighborDirection::South));
  EXPECT_EQ(face_range.getFaceId(YDirTag(), 1, 2), makeFaceId(1, 1, 1, 2, NeighborDirection::North, NeighborDirection::South));
}

TEST(FaceIterPerDirection, YDirectionFaceCount)
{
  Range2D cell_range(0, 3, 0, 4);
  FaceRangePerDirection face_range(cell_range);

  int nfaces = 0;
  for ([[maybe_unused]] UInt i : face_range.getXRange(YDirTag()))
    for ([[maybe_unused]] UInt j : face_range.getYRange(YDirTag()))
      nfaces++;

  EXPECT_EQ(nfaces, 9);
}

TEST(FaceIterPerDirection, YDirectionWithGhosts)
{
  Range2D cell_range(0, 3, 0, 4);
  FaceRangePerDirection face_range(cell_range, true);
  EXPECT_EQ(face_range.getXRange(YDirTag()), Range(0, 3));
  EXPECT_EQ(face_range.getYRange(YDirTag()), Range(0, 5));

  EXPECT_EQ(face_range.getFaceId(YDirTag(), 0, 1), makeFaceId(0, 0, 0, 1, NeighborDirection::North, NeighborDirection::South));
  EXPECT_EQ(face_range.getFaceId(YDirTag(), 1, 2), makeFaceId(1, 1, 1, 2, NeighborDirection::North, NeighborDirection::South));
}

TEST(FaceIterPerDirection, YDirectionWithGhostsFaceCount)
{
  Range2D cell_range(0, 3, 0, 4);
  FaceRangePerDirection face_range(cell_range, true);

  int nfaces = 0;
  for ([[maybe_unused]] UInt i : face_range.getXRange(YDirTag()))
    for ([[maybe_unused]] UInt j : face_range.getYRange(YDirTag()))
      nfaces++;

  EXPECT_EQ(nfaces, 15);
}

TEST(FaceIterPerDirection, Increment)
{
  UInt i = 2;
  UInt j = 3;
  {
    auto [i2, j2] = increment(XDirTag(), i, j, 2);
    EXPECT_EQ(i2, 4);
    EXPECT_EQ(j2, j);
  }

  {
    auto [i2, j2] = increment(XDirTag(), i, j, -2);
    EXPECT_EQ(i2, 0);
    EXPECT_EQ(j2, j);
  }

  {
    auto [i2, j2] = increment(YDirTag(), i, j, 2);
    EXPECT_EQ(i2, i);
    EXPECT_EQ(j2, 5);
  }

  {
    auto [i2, j2] = increment(YDirTag(), i, j, -2);
    EXPECT_EQ(i2, i);
    EXPECT_EQ(j2, 1);
  }  
}

TEST(FaceIterPerDirection, EmptyRangeWithBoundary)
{
  Range2D range(0, 0, 0, 0);
  FaceRangePerDirection face_range(range, true);
  int nfaces = 0;
  for ([[maybe_unused]] UInt i : face_range.getXRange(XDirTag{}))
    for ([[maybe_unused]] UInt j : face_range.getYRange(XDirTag{}))
      nfaces++;

  for ([[maybe_unused]] UInt i : face_range.getXRange(YDirTag{}))
    for ([[maybe_unused]] UInt j : face_range.getYRange(YDirTag{}))
      nfaces++;
  
  EXPECT_EQ(nfaces, 0);
}

TEST(FaceIterPerDirection, EmptyRangeWithOutBoundary)
{
  Range2D range(0, 0, 0, 0);
  FaceRangePerDirection face_range(range, false);
  int nfaces = 0;
  for ([[maybe_unused]] UInt i : face_range.getXRange(XDirTag{}))
    for ([[maybe_unused]] UInt j : face_range.getYRange(XDirTag{}))
      nfaces++;

  for ([[maybe_unused]] UInt i : face_range.getXRange(YDirTag{}))
    for ([[maybe_unused]] UInt j : face_range.getYRange(YDirTag{}))
      nfaces++;
  
  EXPECT_EQ(nfaces, 0);
}


TEST(FaceIterPerDirection, BoundaryFlagsAllIncluded)
{
  Range2D range(0, 3, 0, 4);
  FaceRangeBoundaryFlags boundary_flags;
  FaceRangePerDirection face_range(range, boundary_flags);

  EXPECT_EQ(face_range.getRange(XDirTag{}), Range2D(0, 4, 0, 4));
  EXPECT_EQ(face_range.getRange(YDirTag{}), Range2D(0, 3, 0, 5));
}

TEST(FaceIterPerDirection, BoundaryFlagsAllNoBottom)
{
  Range2D range(0, 3, 0, 4);
  FaceRangeBoundaryFlags boundary_flags;
  boundary_flags.include_bottom = false;
  FaceRangePerDirection face_range(range, boundary_flags);

  EXPECT_EQ(face_range.getRange(XDirTag{}), Range2D(0, 4, 0, 4));
  EXPECT_EQ(face_range.getRange(YDirTag{}), Range2D(0, 3, 1, 5));
}

TEST(FaceIterPerDirection, BoundaryFlagsAllNoLeft)
{
  Range2D range(0, 3, 0, 4);
  FaceRangeBoundaryFlags boundary_flags;
  boundary_flags.include_left = false;
  FaceRangePerDirection face_range(range, boundary_flags);

  EXPECT_EQ(face_range.getRange(XDirTag{}), Range2D(1, 4, 0, 4));
  EXPECT_EQ(face_range.getRange(YDirTag{}), Range2D(0, 3, 0, 5));
}

TEST(FaceIterPerDirection, BoundaryFlagsAllExcluded)
{
  Range2D range(0, 3, 0, 4);
  FaceRangeBoundaryFlags boundary_flags;
  boundary_flags.include_bottom = false;
  boundary_flags.include_top    = false;
  boundary_flags.include_left   = false;
  boundary_flags.include_right  = false;
  FaceRangePerDirection face_range(range, boundary_flags);

  EXPECT_EQ(face_range.getRange(XDirTag{}), Range2D(1, 3, 0, 4));
  EXPECT_EQ(face_range.getRange(YDirTag{}), Range2D(0, 3, 1, 4));
}

TEST(FaceIterPerDirection, BoundaryFlagsEmptyRangeWithOutBoundary)
{
  Range2D range(0, 0, 0, 0);
  FaceRangeBoundaryFlags boundary_flags;
  boundary_flags.include_bottom = false;
  boundary_flags.include_top    = false;
  boundary_flags.include_left   = false;
  boundary_flags.include_right  = false;  
  FaceRangePerDirection face_range(range, boundary_flags);
  int nfaces = 0;
  for ([[maybe_unused]] UInt i : face_range.getXRange(XDirTag{}))
    for ([[maybe_unused]] UInt j : face_range.getYRange(XDirTag{}))
      nfaces++;

  for ([[maybe_unused]] UInt i : face_range.getXRange(YDirTag{}))
    for ([[maybe_unused]] UInt j : face_range.getYRange(YDirTag{}))
      nfaces++;
  
  EXPECT_EQ(nfaces, 0);
}

//TODO: is first constructor still needed