#include "gtest/gtest.h"
#include "mesh/adjacent_block_indexer.h"
#include "utils/project_defs.h"
#include "mesh/structured_block_interface.h"

using namespace structured_fv;
using namespace structured_fv::mesh;

namespace {

std::array<UInt, 2> make_array(const std::array<UInt, 2>& arr)
{
  return arr;
}
}

TEST(StructuredBlockInterface, East)
{
  UInt left_block_id = 0;
  UInt left_block_constant_index = 3;
  Range left_block_variable_index(0, 6);
  NeighborDirection dir = NeighborDirection::East;
  std::array<Int, 2> transform{1, 2};

  UInt right_block_id = 1;
  std::array<UInt, 2> right_block_min_cell{0, 0};

  StructuredBlockInterface iface(left_block_id, left_block_constant_index, left_block_variable_index,
                                 dir, transform, right_block_id, right_block_min_cell);

  EXPECT_EQ(iface.getLeftBlockId(), left_block_id);
  EXPECT_EQ(iface.getRightBlockId(), right_block_id);
  EXPECT_EQ(iface.getNeighborDirection(), dir);
  EXPECT_EQ(iface.getLeftBlockBoundaryCells(), Range2D(3, 4, 0, 6));

  const AdjacentBlockIndexer& indexer = iface.getAdjacentBlockIndexer();
  EXPECT_EQ(indexer(left_block_constant_index+1, left_block_variable_index(0)), make_array({0, 0}));
  EXPECT_EQ(indexer(left_block_constant_index+1, left_block_variable_index(1)), make_array({0, 1}));
  EXPECT_EQ(indexer(left_block_constant_index+2, left_block_variable_index(0)), make_array({1, 0}));  
}

TEST(StructuredBlockInterface, North)
{
  UInt left_block_id = 0;
  UInt left_block_constant_index = 5;
  Range left_block_variable_index(0, 4);
  NeighborDirection dir = NeighborDirection::North;
  std::array<Int, 2> transform{1, 2};

  UInt right_block_id = 1;
  std::array<UInt, 2> right_block_min_cell{0, 0};

  StructuredBlockInterface iface(left_block_id, left_block_constant_index, left_block_variable_index,
                                 dir, transform, right_block_id, right_block_min_cell);

  EXPECT_EQ(iface.getLeftBlockId(), left_block_id);
  EXPECT_EQ(iface.getRightBlockId(), right_block_id);
  EXPECT_EQ(iface.getNeighborDirection(), dir);
  EXPECT_EQ(iface.getLeftBlockBoundaryCells(), Range2D(0, 4, 5, 6));

  const AdjacentBlockIndexer& indexer = iface.getAdjacentBlockIndexer();
  EXPECT_EQ(indexer(left_block_variable_index(0), left_block_constant_index+1), make_array({0, 0}));
  EXPECT_EQ(indexer(left_block_variable_index(1), left_block_constant_index+1), make_array({1, 0}));
  EXPECT_EQ(indexer(left_block_variable_index(0), left_block_constant_index+2), make_array({0, 1}));  
}