#include "gtest/gtest.h"
#include "disc/disc_interface.h"
#include "disc/discretization.h"
#include "mesh/adjacent_block_indexer.h"
#include "mesh/structured_mesh.h"

using namespace structured_fv;

namespace {
class DiscIfaceTester : public ::testing::Test
{
  public:
    DiscIfaceTester()
    {
      mesh::MeshSpec spec(2, 1, m_num_bc_ghost_cells);
      spec.blocks(0, 0) = mesh::MeshBlockSpec(3, 4);
      spec.blocks(1, 0) = mesh::MeshBlockSpec(4, 5, 1);
      m_mesh = std::make_shared<mesh::StructuredMesh>(spec);
      m_disc = std::make_shared<disc::StructuredDisc>(m_mesh, m_num_bc_ghost_cells);
    }

    int m_num_bc_ghost_cells = 2;
    std::shared_ptr<mesh::StructuredMesh> m_mesh;
    std::shared_ptr<disc::StructuredDisc> m_disc;
};

std::array<UInt, 2> make_array2(const std::array<UInt, 2>& arr)
{
  return arr;
}

}

TEST_F(DiscIfaceTester, RegularInterface)
{
  const disc::StructuredBlockInterface& iface = m_disc->getBlockInterface(0);
  EXPECT_EQ(iface.getBlockIdL(), 0);
  EXPECT_EQ(iface.getBlockIdR(), 1);
  EXPECT_EQ(iface.getOwnedBoundaryVertsL(), Range2D(5, 6, 2, 7));
  EXPECT_EQ(iface.getOwnedBoundaryVertsR(), Range2D(2, 7, 7, 8));
  EXPECT_EQ(iface.getOwnedBoundaryCellsL(), Range2D(4, 5, 2, 6));
  EXPECT_EQ(iface.getOwnedBoundaryCellsR(), Range2D(2, 6, 6, 7));
  EXPECT_EQ(iface.getNeighborDirectionL(), mesh::NeighborDirection::East);
  EXPECT_EQ(iface.getNeighborDirectionR(), mesh::NeighborDirection::North);

  {
    UInt i_corner = *iface.getOwnedBoundaryVertsL().getXRange().begin();
    UInt j_corner = *iface.getOwnedBoundaryVertsL().getYRange().begin();
    const mesh::AdjacentBlockIndexer& indexerL = iface.getAdjBlockVertIndexerL();
    EXPECT_EQ(indexerL(i_corner+1, j_corner),   make_array2({2, 7}));
    EXPECT_EQ(indexerL(i_corner+2, j_corner),   make_array2({2, 6}));
    EXPECT_EQ(indexerL(i_corner+1, j_corner+1), make_array2({3, 7}));
  }

  {
    UInt i_corner = *iface.getOwnedBoundaryCellsL().getXRange().begin();
    UInt j_corner = *iface.getOwnedBoundaryCellsL().getYRange().begin();
    const mesh::AdjacentBlockIndexer& indexerL = iface.getAdjBlockCellIndexerL();
    EXPECT_EQ(indexerL(i_corner+1, j_corner),   make_array2({2, 6}));
    EXPECT_EQ(indexerL(i_corner+2, j_corner),   make_array2({2, 5}));
    EXPECT_EQ(indexerL(i_corner+1, j_corner+1), make_array2({3, 6}));
  }

  {
    UInt i_corner = *iface.getOwnedBoundaryVertsR().getXRange().begin();
    UInt j_corner = *iface.getOwnedBoundaryVertsR().getYRange().begin();
    const mesh::AdjacentBlockIndexer& indexerR = iface.getAdjBlockVertIndexerR();
    EXPECT_EQ(indexerR(i_corner, j_corner+1),   make_array2({5, 2}));
    EXPECT_EQ(indexerR(i_corner, j_corner+2),   make_array2({4, 2}));
    EXPECT_EQ(indexerR(i_corner+1, j_corner+1), make_array2({5, 3}));
  }

  {
    UInt i_corner = *iface.getOwnedBoundaryCellsR().getXRange().begin();
    UInt j_corner = *iface.getOwnedBoundaryCellsR().getYRange().begin();
    const mesh::AdjacentBlockIndexer& indexerR = iface.getAdjBlockCellIndexerR();
    EXPECT_EQ(indexerR(i_corner, j_corner+1),   make_array2({4, 2}));
    EXPECT_EQ(indexerR(i_corner, j_corner+2),   make_array2({3, 2}));
    EXPECT_EQ(indexerR(i_corner+1, j_corner+1), make_array2({4, 3}));
  }
}

TEST_F(DiscIfaceTester, GhostInterface)
{
  const disc::StructuredBlockInterface& iface = m_disc->getBlockInterface(4);
  EXPECT_EQ(iface.getBlockIdL(), 0);
  EXPECT_EQ(iface.getBlockIdR(), 5);
  EXPECT_EQ(iface.getOwnedBoundaryVertsL(), Range2D(2, 6, 2, 3));
  EXPECT_EQ(iface.getOwnedBoundaryVertsR(), Range2D(0, 4, 2, 3));  
  EXPECT_EQ(iface.getOwnedBoundaryCellsL(), Range2D(2, 5, 2, 3));
  EXPECT_EQ(iface.getOwnedBoundaryCellsR(), Range2D(0, 3, 1, 2));
  EXPECT_EQ(iface.getNeighborDirectionL(), mesh::NeighborDirection::South);
  EXPECT_EQ(iface.getNeighborDirectionR(), mesh::NeighborDirection::North);

  {
    UInt i_corner = *iface.getOwnedBoundaryVertsL().getXRange().begin();
    UInt j_corner = *iface.getOwnedBoundaryVertsL().getYRange().begin();
    const mesh::AdjacentBlockIndexer& indexerL = iface.getAdjBlockVertIndexerL();
    EXPECT_EQ(indexerL(i_corner, j_corner-1),   make_array2({0, 2}));
    EXPECT_EQ(indexerL(i_corner, j_corner-2),   make_array2({0, 1}));
    EXPECT_EQ(indexerL(i_corner+1, j_corner-1), make_array2({1, 2}));
  }

  {
    UInt i_corner = *iface.getOwnedBoundaryCellsL().getXRange().begin();
    UInt j_corner = *iface.getOwnedBoundaryCellsL().getYRange().begin();
    const mesh::AdjacentBlockIndexer& indexerL = iface.getAdjBlockCellIndexerL();
    EXPECT_EQ(indexerL(i_corner, j_corner-1),   make_array2({0, 1}));
    EXPECT_EQ(indexerL(i_corner, j_corner-2),   make_array2({0, 0}));
    EXPECT_EQ(indexerL(i_corner+1, j_corner-1), make_array2({1, 1}));
  }

  {
    UInt i_corner = *iface.getOwnedBoundaryVertsR().getXRange().begin();
    UInt j_corner = *iface.getOwnedBoundaryVertsR().getYRange().begin();
    const mesh::AdjacentBlockIndexer& indexerR = iface.getAdjBlockVertIndexerR();
    EXPECT_EQ(indexerR(i_corner, j_corner+1),   make_array2({2, 2}));
    EXPECT_EQ(indexerR(i_corner, j_corner+2),   make_array2({2, 3}));
    EXPECT_EQ(indexerR(i_corner+1, j_corner+1), make_array2({3, 2}));
  }

  {
    UInt i_corner = *iface.getOwnedBoundaryCellsR().getXRange().begin();
    UInt j_corner = *iface.getOwnedBoundaryCellsR().getYRange().begin();
    const mesh::AdjacentBlockIndexer& indexerR = iface.getAdjBlockCellIndexerR();
    EXPECT_EQ(indexerR(i_corner, j_corner+1),   make_array2({2, 2}));
    EXPECT_EQ(indexerR(i_corner, j_corner+2),   make_array2({2, 3}));
    EXPECT_EQ(indexerR(i_corner+1, j_corner+1), make_array2({3, 2}));
  }
}