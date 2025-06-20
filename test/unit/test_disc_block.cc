#include "gtest/gtest.h"
#include "disc/disc_block.h"
#include "disc/discretization.h"
#include "mesh/structured_mesh.h"
#include "utils/face_iter_per_direction.h"

using namespace structured_fv;

namespace {
class DiscBlockTester : public ::testing::Test
{
  public:
    DiscBlockTester()
    {
      mesh::MeshSpec spec(2, 1, m_num_bc_ghost_cells);
      spec.blocks(0, 0) = mesh::MeshBlockSpec(3, 4);
      spec.blocks(1, 0) = mesh::MeshBlockSpec(5, 4);
      m_mesh = std::make_shared<mesh::StructuredMesh>(spec);
      m_disc = std::make_shared<disc::StructuredDisc>(m_mesh, m_num_bc_ghost_cells, 1);
    }

    int m_num_bc_ghost_cells = 2;
    std::shared_ptr<mesh::StructuredMesh> m_mesh;
    std::shared_ptr<disc::StructuredDisc> m_disc;
};

FixedVec<UInt, 2> make_array2(const FixedVec<UInt, 2>& arr)
{
  return arr;
}

FixedVec<UInt, 4> make_array4(const FixedVec<UInt, 4>& arr)
{
  return arr;
}
}

TEST_F(DiscBlockTester, Counts)
{

  {
    const disc::StructuredBlock& block = m_disc->getBlock(0);
    EXPECT_EQ(block.getBlockId(), 0);
    EXPECT_EQ(block.getBlockType(), mesh::BlockType::Regular);
    EXPECT_EQ(block.getCellDimensions(), make_array2({7, 8}));
    EXPECT_EQ(block.getNumGhostCellsPerDirection(), make_array4({2, 2, 2, 2}));
    EXPECT_EQ(block.getOwnedVerts(), Range2D(2, 6, 2, 7));
    EXPECT_EQ(block.getOwnedCells(), Range2D(2, 5, 2, 6));
    EXPECT_EQ(block.getOwnedAndGhostVerts(), Range2D(0, 8, 0, 9));
    EXPECT_EQ(block.getOwnedAndGhostCells(), Range2D(0, 7, 0, 8));
    EXPECT_EQ(block.getOwnedFaces().getRange(XDirTag()), Range2D(2, 6, 2, 6));
    EXPECT_EQ(block.getOwnedFaces().getRange(YDirTag()), Range2D(2, 5, 2, 7));
    EXPECT_EQ(block.getOwnedAndGhostXFaces().getRange(XDirTag()), Range2D(0, 8, 2, 6));
    EXPECT_EQ(block.getOwnedAndGhostYFaces().getRange(YDirTag()), Range2D(2, 5, 0, 9));
  }

  {
    const disc::StructuredBlock& block = m_disc->getGhostBlock(0);
    EXPECT_EQ(block.getBlockId(), 2);
    EXPECT_EQ(block.getBlockType(), mesh::BlockType::GhostBC);
    EXPECT_EQ(block.getCellDimensions(), make_array2({3, 4}));
    EXPECT_EQ(block.getNumGhostCellsPerDirection(), make_array4({0, 0, 2, 0}));
    EXPECT_EQ(block.getOwnedVerts(), Range2D(0, 4, 2, 5));
    EXPECT_EQ(block.getOwnedCells(), Range2D(0, 3, 2, 4));
    EXPECT_EQ(block.getOwnedAndGhostVerts(), Range2D(0, 4, 0, 5));
    EXPECT_EQ(block.getOwnedAndGhostCells(), Range2D(0, 3, 0, 4));
    EXPECT_EQ(block.getOwnedFaces().getRange(XDirTag()), Range2D(0, 4, 2, 4));
    EXPECT_EQ(block.getOwnedFaces().getRange(YDirTag()), Range2D(0, 3, 2, 5));
    EXPECT_EQ(block.getOwnedAndGhostXFaces().getRange(XDirTag()), Range2D(0, 4, 2, 4));
    EXPECT_EQ(block.getOwnedAndGhostYFaces().getRange(YDirTag()), Range2D(0, 3, 0, 5));   
  }  
}