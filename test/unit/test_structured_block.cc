#include "gtest/gtest.h"
#include "mesh/structured_block.h"

using namespace structured_fv;
using namespace structured_fv::mesh;

namespace {

std::array<UInt, 2> make_array(const std::array<UInt, 2>& arr)
{
  return arr;
}
}

TEST(StructuredBlock, SimpleBlock)
{
  MeshBlockSpec spec(2, 4, 0, 0.0, 1.0, 0.0, 2.0);
  StructuredBlock block(spec, 2, BlockType::Regular);

  EXPECT_EQ(block.getBlockId(), 2);
  EXPECT_EQ(block.getBlockType(), BlockType::Regular);
  EXPECT_EQ(block.getOwnedVerts(), Range2D(0, 2, 0, 4));
  EXPECT_EQ(block.getOwnedCells(), Range2D(0, 1, 0, 3));
  EXPECT_EQ(block.getOffsetIntoBlock(), make_array({0, 0}));
  EXPECT_EQ(block.getAllBlockSize(), make_array({2, 4}));

  auto coords = block.getOwnedVertCoords();
  EXPECT_EQ(coords.extent(0), 3);
  EXPECT_EQ(coords.extent(1), 5);
  EXPECT_EQ(coords.extent(2), 2);
  
  EXPECT_DOUBLE_EQ(coords(0, 0, 0), 0);
  EXPECT_DOUBLE_EQ(coords(0, 0, 1), 0);

  EXPECT_DOUBLE_EQ(coords(1, 0, 0), 0.5);
  EXPECT_DOUBLE_EQ(coords(1, 0, 1), 0); 

  EXPECT_DOUBLE_EQ(coords(2, 0, 0), 1.0);
  EXPECT_DOUBLE_EQ(coords(2, 0, 1), 0);  


  EXPECT_DOUBLE_EQ(coords(0, 1, 0), 0);
  EXPECT_DOUBLE_EQ(coords(0, 1, 1), 0.5);

  EXPECT_DOUBLE_EQ(coords(1, 1, 0), 0.5);
  EXPECT_DOUBLE_EQ(coords(1, 1, 1), 0.5); 

  EXPECT_DOUBLE_EQ(coords(2, 1, 0), 1.0);
  EXPECT_DOUBLE_EQ(coords(2, 1, 1), 0.5);    
}