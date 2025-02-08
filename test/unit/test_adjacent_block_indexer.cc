#include "gtest/gtest.h"
#include "mesh/adjacent_block_indexer.h"

using namespace structured_fv;
using namespace structured_fv::mesh;

namespace {

std::array<UInt, 2> make_array(const std::array<UInt, 2>& arr)
{
  return arr;
}

std::array<Int, 2> make_array_signed(const std::array<Int, 2>& arr)
{
  return arr;
}
}

TEST(AdjacentBlockIndexer, SameCoordSystemEast)
{
  // left block is 4 x 3, // right block is 5 x 3
  std::array<Int, 2> transform = {1, 2};
  std::array<UInt, 2> left_block_min_cell = {3, 0};
  NeighborDirection direction = NeighborDirection::East;
  std::array<UInt, 2> right_block_min_cell = {0, 0};
  AdjacentBlockIndexer indexer(transform, left_block_min_cell, direction, right_block_min_cell);

  EXPECT_EQ(indexer(4, 0), make_array({0, 0}));
  EXPECT_EQ(indexer(5, 0), make_array({1, 0}));
  EXPECT_EQ(indexer(4, 1), make_array({0, 1}));

  std::array<Int, 2> idx_signed{4, 0};
  std::array<UInt, 2> idx_unsigned{4, 0};
  EXPECT_EQ(indexer(idx_signed), make_array({0, 0}));
  EXPECT_EQ(indexer(idx_unsigned), make_array({0, 0}));

}

TEST(AdjacentBlockIndexer, SameCoordSystemNorth)
{
  // bottom block is 4 x 3, // top block is 4 x 5
  std::array<Int, 2> transform = {1, 2};
  std::array<UInt, 2> left_block_min_cell = {0, 2};
  NeighborDirection direction = NeighborDirection::North;
  std::array<UInt, 2> right_block_min_cell = {0, 0};
  AdjacentBlockIndexer indexer(transform, left_block_min_cell, direction, right_block_min_cell);

  EXPECT_EQ(indexer(0, 3), make_array({0, 0}));
  EXPECT_EQ(indexer(1, 3), make_array({1, 0}));
  EXPECT_EQ(indexer(0, 4), make_array({0, 1}));
}

TEST(AdjacentBlockIndexer, SameCoordSystemSouth)
{
  // top block is 4 x 3, // bottom block is 4 x 5
  std::array<Int, 2> transform = {1, 2};
  std::array<UInt, 2> left_block_min_cell = {0, 0};
  NeighborDirection direction = NeighborDirection::South;
  std::array<UInt, 2> right_block_min_cell = {0, 4};
  AdjacentBlockIndexer indexer(transform, left_block_min_cell, direction, right_block_min_cell);

  EXPECT_EQ(indexer(0, -1), make_array({0, 4}));
  EXPECT_EQ(indexer(1, -1), make_array({1, 4}));
  EXPECT_EQ(indexer(0, -2), make_array({0, 3}));
}

TEST(AdjacentBlockIndexer, SameCoordSystemWest)
{
  // left block is 4 x 3, // right block is 5 x 3
  std::array<Int, 2> transform = {1, 2};
  std::array<UInt, 2> left_block_min_cell = {0, 0};
  NeighborDirection direction = NeighborDirection::West;
  std::array<UInt, 2> right_block_min_cell = {3, 0};
  AdjacentBlockIndexer indexer(transform, left_block_min_cell, direction, right_block_min_cell);

  EXPECT_EQ(indexer(-1, 0), make_array({3, 0}));
  EXPECT_EQ(indexer(-2, 0), make_array({2, 0}));
  EXPECT_EQ(indexer(-1, 1), make_array({3, 1}));
}

TEST(AdjacentBlockIndexer, ReversedCoordSystemEast)
{
  // right block has origin in top right corner, 
  // left block is 4 x 3, // right block is 3 x 5
  std::array<Int, 2> transform = {-2, -1};
  std::array<UInt, 2> left_block_min_cell = {3, 0};
  NeighborDirection direction = NeighborDirection::East;
  std::array<UInt, 2> right_block_min_cell = {2, 4};
  AdjacentBlockIndexer indexer(transform, left_block_min_cell, direction, right_block_min_cell);

  EXPECT_EQ(indexer(4, 0), make_array({2, 4}));
  EXPECT_EQ(indexer(5, 0), make_array({2, 3}));
  EXPECT_EQ(indexer(4, 1), make_array({1, 4}));
}


TEST(AdjacentBlockIndexer, ReversedCoordSystemNorth)
{
  // top block has origin in top right corner, 
  // bottom block is 4 x 3, // right block is 5 x 4
  std::array<Int, 2> transform = {-2, -1};
  std::array<UInt, 2> left_block_min_cell = {0, 2};
  NeighborDirection direction = NeighborDirection::North;
  std::array<UInt, 2> right_block_min_cell = {4, 3};
  AdjacentBlockIndexer indexer(transform, left_block_min_cell, direction, right_block_min_cell);

  EXPECT_EQ(indexer(0, 3), make_array({4, 3}));
  EXPECT_EQ(indexer(0, 4), make_array({3, 3}));
  EXPECT_EQ(indexer(1, 3), make_array({4, 2}));
}

TEST(AdjacentBlockIndexer, ReversedCoordSystemWest)
{
  // right block has origin in top right corner, 
  // left block is 4 x 3, // right block is 3 x 5
  std::array<Int, 2> transform = {-2, -1};
  std::array<UInt, 2> left_block_min_cell = {0, 0};
  NeighborDirection direction = NeighborDirection::West;
  std::array<UInt, 2> right_block_min_cell = {2, 0};
  AdjacentBlockIndexer indexer(transform, left_block_min_cell, direction, right_block_min_cell);

  EXPECT_EQ(indexer(-1, 0), make_array({2, 0}));
  EXPECT_EQ(indexer(-2, 0), make_array({2, 1}));
  EXPECT_EQ(indexer(-1, 1), make_array({1, 0}));
}

TEST(AdjacentBlockIndexer, ReversedCoordSystemSouth)
{
  // top block has origin in top right corner, 
  // bottom block is 4 x 3, // right block is 5 x 4
  std::array<Int, 2> transform = {-2, -1};
  std::array<UInt, 2> left_block_min_cell = {0, 0};
  NeighborDirection direction = NeighborDirection::South;
  std::array<UInt, 2> right_block_min_cell = {0, 3};
  AdjacentBlockIndexer indexer(transform, left_block_min_cell, direction, right_block_min_cell);

  EXPECT_EQ(indexer(0, -1), make_array({0, 3}));
  EXPECT_EQ(indexer(0, -2), make_array({1, 3}));
  EXPECT_EQ(indexer(1, -1), make_array({0, 2}));
}

TEST(AdjacentBlockIndexer, ToSignedAxis)
{
  EXPECT_EQ(toSignedAxis(NeighborDirection::North),  2);
  EXPECT_EQ(toSignedAxis(NeighborDirection::East),   1);
  EXPECT_EQ(toSignedAxis(NeighborDirection::South), -2);
  EXPECT_EQ(toSignedAxis(NeighborDirection::West),  -1);
}

TEST(AdjacentBlockIndexer, ToNeighborDirection)
{
  EXPECT_EQ(toNeighborDirection(2),  NeighborDirection::North);
  EXPECT_EQ(toNeighborDirection(1),  NeighborDirection::East);
  EXPECT_EQ(toNeighborDirection(-2), NeighborDirection::South);
  EXPECT_EQ(toNeighborDirection(-1), NeighborDirection::West);
}

TEST(AdjacentBlockIndexer, getNeighborImage)
{
  EXPECT_EQ(getNeighborImage(NeighborDirection::East, { 1,  2}), NeighborDirection::West);
  EXPECT_EQ(getNeighborImage(NeighborDirection::East, {-2,  1}), NeighborDirection::North);
  EXPECT_EQ(getNeighborImage(NeighborDirection::East, {-1, -2}), NeighborDirection::East);
  EXPECT_EQ(getNeighborImage(NeighborDirection::East, { 2, -1}), NeighborDirection::South);  

  EXPECT_EQ(getNeighborImage(NeighborDirection::West, { 1,  2}), NeighborDirection::East);
  EXPECT_EQ(getNeighborImage(NeighborDirection::West, {-2,  1}), NeighborDirection::South);
  EXPECT_EQ(getNeighborImage(NeighborDirection::West, {-1, -2}), NeighborDirection::West);
  EXPECT_EQ(getNeighborImage(NeighborDirection::West, { 2, -1}), NeighborDirection::North);  
}

TEST(AdjacentBlockIndexer, ComputeIndices)
{
  EXPECT_EQ(computeIndices(NeighborDirection::North, 2, 2, 3),  make_array_signed({2, 5}));
  EXPECT_EQ(computeIndices(NeighborDirection::East,  2, 2, 3),  make_array_signed({4, 3}));
  EXPECT_EQ(computeIndices(NeighborDirection::South, 2, 2, 3),  make_array_signed({2, 1}));
  EXPECT_EQ(computeIndices(NeighborDirection::West,  2, 2, 3),  make_array_signed({0, 3}));
}

TEST(AdjacentBlockIndexer, GetInverseTransform)
{
  EXPECT_EQ(getInverseTransform({1, 2}), make_array_signed({1, 2}));
  EXPECT_EQ(getInverseTransform({-1, 2}), make_array_signed({-1, 2}));
  EXPECT_EQ(getInverseTransform({1, -2}), make_array_signed({1, -2}));
  EXPECT_EQ(getInverseTransform({2, 1}), make_array_signed({2, 1}));
  EXPECT_EQ(getInverseTransform({2, -1}), make_array_signed({-2, 1}));
}