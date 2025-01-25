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
  auto coord_func = [](double x, double y) { return std::array<double, 2>{x, 4*y}; };
  MeshBlockSpec spec(2, 4, 0, coord_func);
  StructuredBlock block(spec, 2, BlockType::Regular);

  EXPECT_EQ(block.getBlockId(), 2);
  EXPECT_EQ(block.getBlockType(), BlockType::Regular);
  EXPECT_EQ(block.getOwnedVerts(), Range2D(0, 3, 0, 5));
  EXPECT_EQ(block.getOwnedCells(), Range2D(0, 2, 0, 4));
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
  EXPECT_DOUBLE_EQ(coords(0, 1, 1), 1.0);

  EXPECT_DOUBLE_EQ(coords(1, 1, 0), 0.5);
  EXPECT_DOUBLE_EQ(coords(1, 1, 1), 1.0); 

  EXPECT_DOUBLE_EQ(coords(2, 1, 0), 1.0);
  EXPECT_DOUBLE_EQ(coords(2, 1, 1), 1.0);    
}

/*
TEST(StructuredBlock, applyRotation0)
{
  UInt rotation = 0;
  UInt dimx = 5;
  UInt dimy = 4;
  EXPECT_EQ(applyRotation(rotation, 0, 0, dimx, dimy), std::make_pair(UInt(0), UInt(0)));
  EXPECT_EQ(applyRotation(rotation, 4, 0, dimx, dimy), std::make_pair(UInt(4), UInt(0)));
  EXPECT_EQ(applyRotation(rotation, 0, 3, dimx, dimy), std::make_pair(UInt(0), UInt(3)));
  EXPECT_EQ(applyRotation(rotation, 2, 1, dimx, dimy), std::make_pair(UInt(2), UInt(1)));
}

TEST(StructuredBlock, applyRotation1)
{
  UInt rotation = 1;
  UInt dimx = 5;
  UInt dimy = 4;
  EXPECT_EQ(applyRotation(rotation, 0, 0, dimx, dimy), std::make_pair(UInt(0), UInt(4)));
  EXPECT_EQ(applyRotation(rotation, 4, 0, dimx, dimy), std::make_pair(UInt(0), UInt(0)));
  EXPECT_EQ(applyRotation(rotation, 0, 3, dimx, dimy), std::make_pair(UInt(3), UInt(4)));
  EXPECT_EQ(applyRotation(rotation, 2, 1, dimx, dimy), std::make_pair(UInt(1), UInt(2)));
}

TEST(StructuredBlock, applyRotation2)
{
  UInt rotation = 2;
  UInt dimx = 5;
  UInt dimy = 4;
  EXPECT_EQ(applyRotation(rotation, 0, 0, dimx, dimy), std::make_pair(UInt(4), UInt(3)));
  EXPECT_EQ(applyRotation(rotation, 4, 0, dimx, dimy), std::make_pair(UInt(0), UInt(3)));
  EXPECT_EQ(applyRotation(rotation, 0, 3, dimx, dimy), std::make_pair(UInt(4), UInt(0)));
  EXPECT_EQ(applyRotation(rotation, 2, 1, dimx, dimy), std::make_pair(UInt(2), UInt(2)));
}

TEST(StructuredBlock, applyRotation3)
{
  UInt rotation = 3;
  UInt dimx = 5;
  UInt dimy = 4;
  EXPECT_EQ(applyRotation(rotation, 0, 0, dimx, dimy), std::make_pair(UInt(3), UInt(0)));
  EXPECT_EQ(applyRotation(rotation, 4, 0, dimx, dimy), std::make_pair(UInt(3), UInt(4)));
  EXPECT_EQ(applyRotation(rotation, 0, 3, dimx, dimy), std::make_pair(UInt(0), UInt(0)));
  EXPECT_EQ(applyRotation(rotation, 2, 1, dimx, dimy), std::make_pair(UInt(2), UInt(2)));
}
*/

TEST(StructuredBlock, SimpleBlockTransform1)
{
  auto coord_func = [](double x, double y) { return std::array<double, 2>{x, 4*y}; };
  MeshBlockSpec spec(2, 4, 1, coord_func);
  StructuredBlock block(spec, 2, BlockType::Regular);

  EXPECT_EQ(block.getBlockId(), 2);
  EXPECT_EQ(block.getBlockType(), BlockType::Regular);
  EXPECT_EQ(block.getOwnedVerts(), Range2D(0, 3, 0, 5));
  EXPECT_EQ(block.getOwnedCells(), Range2D(0, 2, 0, 4));
  EXPECT_EQ(block.getOffsetIntoBlock(), make_array({0, 0}));
  EXPECT_EQ(block.getAllBlockSize(), make_array({2, 4}));

  auto coords = block.getOwnedVertCoords();
  EXPECT_EQ(coords.extent(0), 3);
  EXPECT_EQ(coords.extent(1), 5);
  EXPECT_EQ(coords.extent(2), 2);

  EXPECT_DOUBLE_EQ(coords(0, 0, 0), 1);
  EXPECT_DOUBLE_EQ(coords(0, 0, 1), 0);

  EXPECT_DOUBLE_EQ(coords(2, 4, 0), 0);
  EXPECT_DOUBLE_EQ(coords(2, 4, 1), 4);  
}

TEST(StructuredBlock, SimpleBlockTransform2)
{
  auto coord_func = [](double x, double y) { return std::array<double, 2>{x, 4*y}; };
  MeshBlockSpec spec(2, 4, 2, coord_func);
  StructuredBlock block(spec, 2, BlockType::Regular);

  EXPECT_EQ(block.getBlockId(), 2);
  EXPECT_EQ(block.getBlockType(), BlockType::Regular);
  EXPECT_EQ(block.getOwnedVerts(), Range2D(0, 3, 0, 5));
  EXPECT_EQ(block.getOwnedCells(), Range2D(0, 2, 0, 4));
  EXPECT_EQ(block.getOffsetIntoBlock(), make_array({0, 0}));
  EXPECT_EQ(block.getAllBlockSize(), make_array({2, 4}));

  auto coords = block.getOwnedVertCoords();
  EXPECT_EQ(coords.extent(0), 3);
  EXPECT_EQ(coords.extent(1), 5);
  EXPECT_EQ(coords.extent(2), 2);

  EXPECT_DOUBLE_EQ(coords(0, 0, 0), 1);
  EXPECT_DOUBLE_EQ(coords(0, 0, 1), 4);

  EXPECT_DOUBLE_EQ(coords(2, 4, 0), 0);
  EXPECT_DOUBLE_EQ(coords(2, 4, 1), 0);  
}

TEST(StructuredBlock, SimpleBlockTransform3)
{
  auto coord_func = [](double x, double y) { return std::array<double, 2>{x, 4*y}; };
  MeshBlockSpec spec(2, 4, 3, coord_func);
  StructuredBlock block(spec, 2, BlockType::Regular);

  EXPECT_EQ(block.getBlockId(), 2);
  EXPECT_EQ(block.getBlockType(), BlockType::Regular);
  EXPECT_EQ(block.getOwnedVerts(), Range2D(0, 3, 0, 5));
  EXPECT_EQ(block.getOwnedCells(), Range2D(0, 2, 0, 4));
  EXPECT_EQ(block.getOffsetIntoBlock(), make_array({0, 0}));
  EXPECT_EQ(block.getAllBlockSize(), make_array({2, 4}));

  auto coords = block.getOwnedVertCoords();
  EXPECT_EQ(coords.extent(0), 3);
  EXPECT_EQ(coords.extent(1), 5);
  EXPECT_EQ(coords.extent(2), 2);

  EXPECT_DOUBLE_EQ(coords(0, 0, 0), 0);
  EXPECT_DOUBLE_EQ(coords(0, 0, 1), 4);

  EXPECT_DOUBLE_EQ(coords(2, 4, 0), 1);
  EXPECT_DOUBLE_EQ(coords(2, 4, 1), 0);  
}

