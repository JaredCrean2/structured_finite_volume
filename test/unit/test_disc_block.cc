#include "gtest/gtest.h"
#include "disc/disc_block.h"
#include "disc/discretization.h"
#include "mesh/structured_mesh.h"
#include "utils/face_iter_per_direction.h"
#include "utils/project_defs.h"

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
    EXPECT_EQ(block.getBlockType(), BlockType::Regular);
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
    EXPECT_EQ(block.getOwnedFacesExceptBC().getRange(XDirTag{}), Range2D(3, 6, 2, 6));
    EXPECT_EQ(block.getOwnedFacesExceptBC().getRange(YDirTag{}), Range2D(2, 5, 3, 6));
  }

  {
    const disc::StructuredBlock& block = m_disc->getGhostBlock(0);
    EXPECT_EQ(block.getBlockId(), 2);
    EXPECT_EQ(block.getBlockType(), BlockType::GhostBC);
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

TEST_F(DiscBlockTester, InteriorFacesPlusBoundaryFaces)
{

  UInt block_id = 0;
  const disc::StructuredBlock& block = m_disc->getBlock(block_id);
  Kokkos::View<int**, HostMemorySpace> xface_counts("xface_counts", block.getCellDimensions()[0]+1, block.getCellDimensions()[1]);
  Kokkos::View<int**, HostMemorySpace> yface_counts("yface_counts", block.getCellDimensions()[0], block.getCellDimensions()[1]+1);

  FaceRangePerDirection face_range = block.getOwnedFacesExceptBC();
  for (UInt i : face_range.getXRange(XDirTag{}))
    for (UInt j : face_range.getYRange(XDirTag{}))
      xface_counts(i, j) += 1;

  for (UInt i : face_range.getXRange(YDirTag{}))
    for (UInt j : face_range.getYRange(YDirTag{}))
      yface_counts(i, j) += 1;

  for (Int iface_id : m_disc->getBlockInterfaces(block_id))
    if (iface_id >= 0)
    {
      const disc::StructuredBlockInterface& iface = m_disc->getBlockInterface(iface_id);
      UInt other_block_id = iface.getBlockIdL() == block_id ? iface.getBlockIdR() : iface.getBlockIdL();
      if (m_disc->getBlock(other_block_id).getBlockType() == BlockType::GhostBC)
      {
        for (UInt i : iface.getOwnedBoundaryFacesL().getXRange())
          for (UInt j : iface.getOwnedBoundaryFacesL().getYRange())
          {
            if (iface.getNeighborDirectionL() == NeighborDirection::East ||
                iface.getNeighborDirectionL() == NeighborDirection::West)
            {
              xface_counts(i, j) += 1;
            } else
            {
              yface_counts(i, j) += 1;
            }
          }
      }
    }

  FaceRangePerDirection all_face_range = block.getOwnedFaces();
  for (UInt i : all_face_range.getXRange(XDirTag{}))
    for (UInt j : all_face_range.getYRange(XDirTag{}))
    {
      EXPECT_EQ(xface_counts(i, j), 1);
    }
   
  for (UInt i : all_face_range.getXRange(YDirTag{}))
    for (UInt j : all_face_range.getYRange(YDirTag{}))
      EXPECT_EQ(yface_counts(i, j), 1);      
}