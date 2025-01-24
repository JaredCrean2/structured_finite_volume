#include "mesh/structured_block.h"
#include "mesh/structured_block_interface.h"
#include "mesh/structured_mesh.h"
#include "gtest/gtest.h"

using namespace structured_fv;
using namespace structured_fv::mesh;

namespace
{
void checkCoords(double x0, double y0, double dx, double dy, StructuredBlock::CoordsHostView coords)
{
  for (UInt i=0; i < coords.extent(0); ++i)
    for (UInt j=0; j < coords.extent(1); ++j)
    {
      double x = x0 + i * dx;
      double y = y0 + j * dy;
      EXPECT_NEAR(coords(i, j, 0), x, 1e-13);
      EXPECT_NEAR(coords(i, j, 1), y, 1e-13);
    }
}
}

TEST(StructuredMesh, SingleBlock) {
  MeshSpec meshspec(1, 1, 2);
  meshspec.blocks(0, 0) = MeshBlockSpec(2, 3, 0);

  StructuredMesh mesh(meshspec);

  EXPECT_EQ(mesh.getNumBlocks(), 5);
  EXPECT_EQ(mesh.getNumRegularBlocks(), 1);
  EXPECT_EQ(mesh.getNumGhostBCBlocks(), 4);
  EXPECT_EQ(mesh.getNumBlockInterfaces(), 4);
  EXPECT_EQ(mesh.getNumRegularBlockInterfaces(), 0);
  EXPECT_EQ(mesh.getNumGhostBCBlockInterfaces(), 4);

  {
    const StructuredBlock& block = mesh.getBlock(0);
    EXPECT_EQ(block.getBlockId(), 0);
    EXPECT_EQ(block.getBlockType(), BlockType::Regular);
    EXPECT_EQ(block.getOwnedVerts(), Range2D(0, 3, 0, 4));
    EXPECT_EQ(block.getOwnedCells(), Range2D(0, 2, 0, 3));
    EXPECT_EQ(block.getOffsetIntoBlock()[0], 0);
    EXPECT_EQ(block.getOffsetIntoBlock()[1], 0);
    EXPECT_EQ(block.getAllBlockSize()[0], 2);
    EXPECT_EQ(block.getAllBlockSize()[1], 3);
    auto coords = block.getOwnedVertCoords();
    EXPECT_EQ(coords.extent(0), 3);
    EXPECT_EQ(coords.extent(1), 4);
    checkCoords(0, 0, 0.5, 1.0/3, coords);
  }

  {
    // north ghost block
    const StructuredBlock& block = mesh.getBlock(1);
    EXPECT_EQ(block.getBlockId(), 1);
    EXPECT_EQ(block.getBlockType(), BlockType::GhostBC);
    EXPECT_EQ(block.getOwnedVerts(), Range2D(0, 3, 0, 3));
    EXPECT_EQ(block.getOwnedCells(), Range2D(0, 2, 0, 2));
    EXPECT_EQ(block.getOffsetIntoBlock()[0], 0);
    EXPECT_EQ(block.getOffsetIntoBlock()[1], 0);
    EXPECT_EQ(block.getAllBlockSize()[0], 2);
    EXPECT_EQ(block.getAllBlockSize()[1], 2);
    auto coords = block.getOwnedVertCoords();
    EXPECT_EQ(coords.extent(0), 3);
    EXPECT_EQ(coords.extent(1), 3);
    checkCoords(0, 1, 0.5, 1.0/3, coords);

    EXPECT_EQ(mesh.getBCBlockRange(0), Range(1, 2));
    const StructuredBlockInterface& iface = mesh.getGhostBCBlockInterface(0);
    EXPECT_EQ(iface.getLeftBlockId(), 0);
    EXPECT_EQ(iface.getRightBlockId(), 1);
    EXPECT_EQ(iface.getNeighborDirection(), NeighborDirection::North);
    EXPECT_EQ(iface.getLeftBlockBoundaryCells(), Range2D(0, 2, 2, 3));  
  }

  {
    // East ghost block
    const StructuredBlock& block = mesh.getBlock(2);
    EXPECT_EQ(block.getBlockId(), 2);
    EXPECT_EQ(block.getBlockType(), BlockType::GhostBC);
    EXPECT_EQ(block.getOwnedVerts(), Range2D(0, 3, 0, 4));
    EXPECT_EQ(block.getOwnedCells(), Range2D(0, 2, 0, 3));
    EXPECT_EQ(block.getOffsetIntoBlock()[0], 0);
    EXPECT_EQ(block.getOffsetIntoBlock()[1], 0);
    EXPECT_EQ(block.getAllBlockSize()[0], 2);
    EXPECT_EQ(block.getAllBlockSize()[1], 3);
    auto coords = block.getOwnedVertCoords();
    EXPECT_EQ(coords.extent(0), 3);
    EXPECT_EQ(coords.extent(1), 4);
    checkCoords(1, 0, 0.5, 1.0/3, coords);

    EXPECT_EQ(mesh.getBCBlockRange(1), Range(2, 3));
    const StructuredBlockInterface& iface = mesh.getGhostBCBlockInterface(1);
    EXPECT_EQ(iface.getLeftBlockId(), 0);
    EXPECT_EQ(iface.getRightBlockId(), 2);
    EXPECT_EQ(iface.getNeighborDirection(), NeighborDirection::East);
    EXPECT_EQ(iface.getLeftBlockBoundaryCells(), Range2D(1, 2, 0, 3));    
  }

  {
    // south ghost block
    const StructuredBlock& block = mesh.getBlock(3);
    EXPECT_EQ(block.getBlockId(), 3);
    EXPECT_EQ(block.getBlockType(), BlockType::GhostBC);
    EXPECT_EQ(block.getOwnedVerts(), Range2D(0, 3, 0, 3));
    EXPECT_EQ(block.getOwnedCells(), Range2D(0, 2, 0, 2));
    EXPECT_EQ(block.getOffsetIntoBlock()[0], 0);
    EXPECT_EQ(block.getOffsetIntoBlock()[1], 0);
    EXPECT_EQ(block.getAllBlockSize()[0], 2);
    EXPECT_EQ(block.getAllBlockSize()[1], 2);
    auto coords = block.getOwnedVertCoords();
    EXPECT_EQ(coords.extent(0), 3);
    EXPECT_EQ(coords.extent(1), 3);
    checkCoords(0, -2.0/3, 0.5, 1.0/3, coords);

    EXPECT_EQ(mesh.getBCBlockRange(2), Range(3, 4));
    const StructuredBlockInterface& iface = mesh.getGhostBCBlockInterface(2);
    EXPECT_EQ(iface.getLeftBlockId(), 0);
    EXPECT_EQ(iface.getRightBlockId(), 3);
    EXPECT_EQ(iface.getNeighborDirection(), NeighborDirection::South);
    EXPECT_EQ(iface.getLeftBlockBoundaryCells(), Range2D(0, 2, 0, 1));    
  }    

  {
    std::cout << "\nchecking south block" << std::endl;
    // West ghost block
    const StructuredBlock& block = mesh.getBlock(4);
    EXPECT_EQ(block.getBlockId(), 4);
    EXPECT_EQ(block.getBlockType(), BlockType::GhostBC);
    EXPECT_EQ(block.getOwnedVerts(), Range2D(0, 3, 0, 4));
    EXPECT_EQ(block.getOwnedCells(), Range2D(0, 2, 0, 3));
    EXPECT_EQ(block.getOffsetIntoBlock()[0], 0);
    EXPECT_EQ(block.getOffsetIntoBlock()[1], 0);
    EXPECT_EQ(block.getAllBlockSize()[0], 2);
    EXPECT_EQ(block.getAllBlockSize()[1], 3);
    auto coords = block.getOwnedVertCoords();
    EXPECT_EQ(coords.extent(0), 3);
    EXPECT_EQ(coords.extent(1), 4);
    checkCoords(-1, 0, 0.5, 1.0/3, coords);

    EXPECT_EQ(mesh.getBCBlockRange(3), Range(4, 5));
    const StructuredBlockInterface& iface = mesh.getGhostBCBlockInterface(3);
    EXPECT_EQ(iface.getLeftBlockId(), 0);
    EXPECT_EQ(iface.getRightBlockId(), 4);
    EXPECT_EQ(iface.getNeighborDirection(), NeighborDirection::West);
    EXPECT_EQ(iface.getLeftBlockBoundaryCells(), Range2D(0, 1, 0, 3));    
  }
}