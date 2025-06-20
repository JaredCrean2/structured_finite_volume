#include "mesh/structured_block.h"
#include "mesh/structured_block_interface.h"
#include "mesh/structured_mesh.h"
#include "gtest/gtest.h"

using namespace structured_fv;
using namespace structured_fv::mesh;

namespace
{
void checkCoords(double x0, double y0, const FixedVec<double, 4>& mat, StructuredBlock::CoordsHostView coords)
{
  for (UInt i=0; i < coords.extent(0); ++i)
    for (UInt j=0; j < coords.extent(1); ++j)
    {
      double x = mat[0] * i + mat[1]*j + x0;
      double y = mat[2] * i + mat[3]*j + y0;
      EXPECT_NEAR(coords(i, j, 0), x, 1e-13);
      EXPECT_NEAR(coords(i, j, 1), y, 1e-13);
    }
}

void checkCoords(double x0, double y0, double dx, double dy, StructuredBlock::CoordsHostView coords)
{
  checkCoords(x0, y0, {dx, 0, 0, dy}, coords);
}

FixedVec<Int, 4> make_array(const FixedVec<Int, 4>& vals)
{
  return vals;
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
    EXPECT_EQ(mesh.getBlockInterfaces(0), make_array({0, 1, 2, 3}));
    for (int i=0; i < 4; ++i)
    {
      const StructuredBlockInterface& iface = mesh.getBlockInterface(mesh.getBlockInterfaces(0)[i]);
      EXPECT_EQ(iface.getBlockIdL(), 0);
      EXPECT_EQ(iface.getBlockIdR(), i+1);
    }
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
    EXPECT_EQ(iface.getBlockIdL(), 0);
    EXPECT_EQ(iface.getBlockIdR(), 1);
    EXPECT_EQ(iface.getNeighborDirectionL(), NeighborDirection::North);
    EXPECT_EQ(iface.getNeighborDirectionR(), NeighborDirection::South);
    EXPECT_EQ(iface.getBoundaryCellsL(), Range2D(0, 2, 2, 3));
    EXPECT_EQ(mesh.getBlockInterfaces(1), make_array({-1, -1, 0, -1}));
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
    EXPECT_EQ(iface.getBlockIdL(), 0);
    EXPECT_EQ(iface.getBlockIdR(), 2);
    EXPECT_EQ(iface.getNeighborDirectionL(), NeighborDirection::East);
    EXPECT_EQ(iface.getNeighborDirectionR(), NeighborDirection::West);
    EXPECT_EQ(iface.getBoundaryCellsL(), Range2D(1, 2, 0, 3));
    EXPECT_EQ(mesh.getBlockInterfaces(2), make_array({-1, -1, -1, 1}));

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
    EXPECT_EQ(iface.getBlockIdL(), 0);
    EXPECT_EQ(iface.getBlockIdR(), 3);
    EXPECT_EQ(iface.getNeighborDirectionL(), NeighborDirection::South);
    EXPECT_EQ(iface.getNeighborDirectionR(), NeighborDirection::North);
    EXPECT_EQ(iface.getBoundaryCellsL(), Range2D(0, 2, 0, 1));
    EXPECT_EQ(mesh.getBlockInterfaces(3), make_array({2, -1, -1, -1}));
  }    

  {
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
    EXPECT_EQ(iface.getBlockIdL(), 0);
    EXPECT_EQ(iface.getBlockIdR(), 4);
    EXPECT_EQ(iface.getNeighborDirectionL(), NeighborDirection::West);
    EXPECT_EQ(iface.getNeighborDirectionR(), NeighborDirection::East);
    EXPECT_EQ(iface.getBoundaryCellsL(), Range2D(0, 1, 0, 3));  
    EXPECT_EQ(mesh.getBlockInterfaces(4), make_array({-1, 3, -1, -1}));
  }
}

TEST(StructuredMesh, SingleBlockRotation1)
{
  MeshSpec meshspec(1, 1, 2);
  meshspec.blocks(0, 0) = MeshBlockSpec(2, 3, 1);

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
    checkCoords(1, 0, {0, -1.0/3, 0.5, 0}, coords);
  }

  {
    // north ghost block (in the blocks coordinate system)
    const StructuredBlock& block = mesh.getBlock(4);
    EXPECT_EQ(block.getBlockId(), 4);
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
    checkCoords(0, 0, {0, -1.0/3, 0.5, 0}, coords);

    EXPECT_EQ(mesh.getBCBlockRange(3), Range(4, 5));
    const StructuredBlockInterface& iface = mesh.getGhostBCBlockInterface(3);
    EXPECT_EQ(iface.getBlockIdL(), 0);
    EXPECT_EQ(iface.getBlockIdR(), 4);
    EXPECT_EQ(iface.getNeighborDirectionL(), NeighborDirection::North);
    EXPECT_EQ(iface.getNeighborDirectionR(), NeighborDirection::South);
    EXPECT_EQ(iface.getBoundaryCellsL(), Range2D(0, 2, 2, 3));  
  }


  {
    // East ghost block
    const StructuredBlock& block = mesh.getBlock(1);
    EXPECT_EQ(block.getBlockId(), 1);
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
    checkCoords(1, 1, {0, -1.0/3, 0.5, 0}, coords);

    EXPECT_EQ(mesh.getBCBlockRange(0), Range(1, 2));
    const StructuredBlockInterface& iface = mesh.getGhostBCBlockInterface(0);
    EXPECT_EQ(iface.getBlockIdL(), 0);
    EXPECT_EQ(iface.getBlockIdR(), 1);
    EXPECT_EQ(iface.getNeighborDirectionL(), NeighborDirection::East);
    EXPECT_EQ(iface.getNeighborDirectionR(), NeighborDirection::West);
    EXPECT_EQ(iface.getBoundaryCellsL(), Range2D(1, 2, 0, 3));    
  }

  {
    // south ghost block
    const StructuredBlock& block = mesh.getBlock(2);
    EXPECT_EQ(block.getBlockId(), 2);
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
    checkCoords(5.0/3, 0, {0, -1.0/3, 0.5, 0}, coords);

    EXPECT_EQ(mesh.getBCBlockRange(1), Range(2, 3));
    const StructuredBlockInterface& iface = mesh.getGhostBCBlockInterface(1);
    EXPECT_EQ(iface.getBlockIdL(), 0);
    EXPECT_EQ(iface.getBlockIdR(), 2);
    EXPECT_EQ(iface.getNeighborDirectionL(), NeighborDirection::South);
    EXPECT_EQ(iface.getNeighborDirectionR(), NeighborDirection::North);
    EXPECT_EQ(iface.getBoundaryCellsL(), Range2D(0, 2, 0, 1));    
  }    

  {
    // West ghost block
    const StructuredBlock& block = mesh.getBlock(3);
    EXPECT_EQ(block.getBlockId(), 3);
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
    checkCoords(1, -1, {0, -1.0/3, 0.5, 0}, coords);

    EXPECT_EQ(mesh.getBCBlockRange(2), Range(3, 4));
    const StructuredBlockInterface& iface = mesh.getGhostBCBlockInterface(2);
    EXPECT_EQ(iface.getBlockIdL(), 0);
    EXPECT_EQ(iface.getBlockIdR(), 3);
    EXPECT_EQ(iface.getNeighborDirectionL(), NeighborDirection::West);
    EXPECT_EQ(iface.getNeighborDirectionR(), NeighborDirection::East);
    EXPECT_EQ(iface.getBoundaryCellsL(), Range2D(0, 1, 0, 3));    
  }
}

TEST(StructuredMesh, TwoBlocks)
{
  MeshSpec meshspec(2, 1, 2);
  meshspec.blocks(0, 0) = MeshBlockSpec(2, 3, 0);
  meshspec.blocks(1, 0) = MeshBlockSpec(2, 3, 0);

  StructuredMesh mesh(meshspec);

  EXPECT_EQ(mesh.getNumBlocks(), 8);
  EXPECT_EQ(mesh.getNumRegularBlocks(), 2);
  EXPECT_EQ(mesh.getNumGhostBCBlocks(), 6);
  EXPECT_EQ(mesh.getNumBlockInterfaces(), 7);
  EXPECT_EQ(mesh.getNumRegularBlockInterfaces(), 1);
  EXPECT_EQ(mesh.getNumGhostBCBlockInterfaces(), 6);

  EXPECT_EQ(mesh.getBlockInterfaces(0), make_array({1, 0, 4, 6}));
  EXPECT_EQ(mesh.getBlockInterfaces(1), make_array({2, 3, 5, 0}));

  EXPECT_EQ(mesh.getBlockInterfaces(2), make_array({-1, -1,  1, -1}));
  EXPECT_EQ(mesh.getBlockInterfaces(3), make_array({-1, -1,  2, -1}));
  EXPECT_EQ(mesh.getBlockInterfaces(4), make_array({-1, -1, -1,  3}));
  EXPECT_EQ(mesh.getBlockInterfaces(5), make_array({ 4, -1, -1, -1}));
  EXPECT_EQ(mesh.getBlockInterfaces(6), make_array({ 5, -1, -1, -1}));
  EXPECT_EQ(mesh.getBlockInterfaces(7), make_array({-1,  6, -1, -1}));  
}