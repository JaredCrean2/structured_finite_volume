#include "gtest/gtest.h"
#include "mesh/neighbor_direction.h"
#include "utils/face_iter_per_direction.h"
#include "utils/face_iterator.h"

using namespace structured_fv;

namespace {

FaceId makeFaceId(UInt cell_i_left, UInt cell_j_left, UInt cell_i_right, UInt cell_j_right,
                  mesh::NeighborDirection dirL, mesh::NeighborDirection dirR)
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

  EXPECT_EQ(face_range.getFaceId(XDirTag(), 1, 0), makeFaceId(0, 0, 1, 0, mesh::NeighborDirection::East, mesh::NeighborDirection::West));
  EXPECT_EQ(face_range.getFaceId(XDirTag(), 2, 1), makeFaceId(1, 1, 2, 1, mesh::NeighborDirection::East, mesh::NeighborDirection::West));
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

  EXPECT_EQ(face_range.getFaceId(XDirTag(), 1, 0), makeFaceId(0, 0, 1, 0, mesh::NeighborDirection::East, mesh::NeighborDirection::West));
  EXPECT_EQ(face_range.getFaceId(XDirTag(), 1, 1), makeFaceId(0, 1, 1, 1, mesh::NeighborDirection::East, mesh::NeighborDirection::West));
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

  EXPECT_EQ(face_range.getFaceId(YDirTag(), 0, 1), makeFaceId(0, 0, 0, 1, mesh::NeighborDirection::North, mesh::NeighborDirection::South));
  EXPECT_EQ(face_range.getFaceId(YDirTag(), 1, 2), makeFaceId(1, 1, 1, 2, mesh::NeighborDirection::North, mesh::NeighborDirection::South));
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

  EXPECT_EQ(face_range.getFaceId(YDirTag(), 0, 1), makeFaceId(0, 0, 0, 1, mesh::NeighborDirection::North, mesh::NeighborDirection::South));
  EXPECT_EQ(face_range.getFaceId(YDirTag(), 1, 2), makeFaceId(1, 1, 1, 2, mesh::NeighborDirection::North, mesh::NeighborDirection::South));
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