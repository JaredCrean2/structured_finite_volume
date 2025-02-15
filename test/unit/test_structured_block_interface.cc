#include "gtest/gtest.h"
#include "mesh/adjacent_block_indexer.h"
#include "mesh/block_spec.h"
#include "mesh/structured_block.h"
#include "utils/project_defs.h"
#include "mesh/structured_block_interface.h"

using namespace structured_fv;
using namespace structured_fv::mesh;

namespace {

std::array<UInt, 2> make_array(const std::array<UInt, 2>& arr)
{
  return arr;
}

class StructuredBlockInterfaceTester : public ::testing::Test
{
  public:
    StructuredBlockInterfaceTester() :
      m_center_block(MeshBlockSpec(3, 4, 0), 0, BlockType::Regular),
      m_east_block(  MeshBlockSpec(6, 4, 0), 1, BlockType::Regular),
      m_north_block( MeshBlockSpec(3, 6, 0), 2, BlockType::Regular),
      m_east_iface(m_center_block, NeighborDirection::East, m_center_block.getOwnedCells().getYRange(), {1, 2}, 
                   m_east_block, m_east_block.getOwnedCells().getYRange()),
      m_north_iface(m_center_block, NeighborDirection::North, m_center_block.getOwnedCells().getXRange(), {1, 2}, 
                   m_north_block, m_north_block.getOwnedCells().getXRange())                   
    {}

  protected:
    StructuredBlock m_center_block;
    StructuredBlock m_east_block;
    StructuredBlock m_north_block;
    StructuredBlockInterface m_east_iface;
    StructuredBlockInterface m_north_iface;
};

class StructuredBlockInterfaceRotationTester : public ::testing::Test
{
  public:
    StructuredBlockInterfaceRotationTester() :
      m_center_block(MeshBlockSpec(3, 4, 0), 0, BlockType::Regular),
      m_east_block(  MeshBlockSpec(4, 6, 1), 1, BlockType::Regular),
      m_north_block( MeshBlockSpec(3, 6, 2), 2, BlockType::Regular),
      m_east_iface(m_center_block, NeighborDirection::East, m_center_block.getOwnedCells().getYRange(), {-2, 1}, 
                   m_east_block, m_east_block.getOwnedCells().getXRange()),
      m_north_iface(m_center_block, NeighborDirection::North, m_center_block.getOwnedCells().getXRange(), {-1, -2}, 
                   m_north_block, m_north_block.getOwnedCells().getXRange())                   
    {}

  protected:
    StructuredBlock m_center_block;
    StructuredBlock m_east_block;
    StructuredBlock m_north_block;
    StructuredBlockInterface m_east_iface;
    StructuredBlockInterface m_north_iface;
};

}


TEST_F(StructuredBlockInterfaceTester, East)
{
  EXPECT_EQ(m_east_iface.getBlockIdL(), 0);
  EXPECT_EQ(m_east_iface.getBlockIdR(), 1);
  EXPECT_EQ(m_east_iface.getNeighborDirectionL(), NeighborDirection::East);
  EXPECT_EQ(m_east_iface.getNeighborDirectionR(), NeighborDirection::West);
  EXPECT_EQ(m_east_iface.getBoundaryCellsL(), Range2D(2, 3, 0, 4));
  EXPECT_EQ(m_east_iface.getBoundaryCellsR(), Range2D(0, 1, 0, 4));

  const auto& indexerL = m_east_iface.getAdjacentBlockIndexerL();
  EXPECT_EQ(indexerL(3, 0), make_array({0, 0}));
  EXPECT_EQ(indexerL(4, 0), make_array({1, 0}));
  EXPECT_EQ(indexerL(3, 1), make_array({0, 1}));

  const auto& indexerR = m_east_iface.getAdjacentBlockIndexerR();
  EXPECT_EQ(indexerR(-1, 0), make_array({2, 0}));
  EXPECT_EQ(indexerR(-2, 0), make_array({1, 0}));
  EXPECT_EQ(indexerR(-1, 1), make_array({2, 1}));  
}

TEST_F(StructuredBlockInterfaceTester, North)
{
  EXPECT_EQ(m_north_iface.getBlockIdL(), 0);
  EXPECT_EQ(m_north_iface.getBlockIdR(), 2);
  EXPECT_EQ(m_north_iface.getNeighborDirectionL(), NeighborDirection::North);
  EXPECT_EQ(m_north_iface.getNeighborDirectionR(), NeighborDirection::South);
  EXPECT_EQ(m_north_iface.getBoundaryCellsL(), Range2D(0, 3, 3, 4));
  EXPECT_EQ(m_north_iface.getBoundaryCellsR(), Range2D(0, 3, 0, 1));

  const auto& indexerL = m_north_iface.getAdjacentBlockIndexerL();
  EXPECT_EQ(indexerL(0, 4), make_array({0, 0}));
  EXPECT_EQ(indexerL(0, 5), make_array({0, 1}));
  EXPECT_EQ(indexerL(1, 4), make_array({1, 0}));

  const auto& indexerR = m_north_iface.getAdjacentBlockIndexerR();
  EXPECT_EQ(indexerR(0, -1), make_array({0, 3}));
  EXPECT_EQ(indexerR(0, -2), make_array({0, 2}));
  EXPECT_EQ(indexerR(1, -1), make_array({1, 3}));  
}


TEST_F(StructuredBlockInterfaceRotationTester, East)
{
  EXPECT_EQ(m_east_iface.getBlockIdL(), 0);
  EXPECT_EQ(m_east_iface.getBlockIdR(), 1);
  EXPECT_EQ(m_east_iface.getNeighborDirectionL(), NeighborDirection::East);
  EXPECT_EQ(m_east_iface.getNeighborDirectionR(), NeighborDirection::North);
  EXPECT_EQ(m_east_iface.getBoundaryCellsL(), Range2D(2, 3, 0, 4));
  EXPECT_EQ(m_east_iface.getBoundaryCellsR(), Range2D(0, 4, 5, 6));

  const auto& indexerL = m_east_iface.getAdjacentBlockIndexerL();
  EXPECT_EQ(indexerL(3, 0), make_array({0, 5}));
  EXPECT_EQ(indexerL(4, 0), make_array({0, 4}));
  EXPECT_EQ(indexerL(3, 1), make_array({1, 5}));

  const auto& indexerR = m_east_iface.getAdjacentBlockIndexerR();
  EXPECT_EQ(indexerR(0, 6), make_array({2, 0}));
  EXPECT_EQ(indexerR(0, 7), make_array({1, 0}));
  EXPECT_EQ(indexerR(1, 6), make_array({2, 1}));  
}

TEST_F(StructuredBlockInterfaceRotationTester, North)
{
  EXPECT_EQ(m_north_iface.getBlockIdL(), 0);
  EXPECT_EQ(m_north_iface.getBlockIdR(), 2);
  EXPECT_EQ(m_north_iface.getNeighborDirectionL(), NeighborDirection::North);
  EXPECT_EQ(m_north_iface.getNeighborDirectionR(), NeighborDirection::North);
  EXPECT_EQ(m_north_iface.getBoundaryCellsL(), Range2D(0, 3, 3, 4));
  EXPECT_EQ(m_north_iface.getBoundaryCellsR(), Range2D(0, 3, 5, 6));

  const auto& indexerL = m_north_iface.getAdjacentBlockIndexerL();
  EXPECT_EQ(indexerL(0, 4), make_array({2, 5}));
  EXPECT_EQ(indexerL(0, 5), make_array({2, 4}));
  EXPECT_EQ(indexerL(1, 4), make_array({1, 5}));

  const auto& indexerR = m_north_iface.getAdjacentBlockIndexerR();
  EXPECT_EQ(indexerR(0, 6), make_array({2, 3}));
  EXPECT_EQ(indexerR(0, 7), make_array({2, 2}));
  EXPECT_EQ(indexerR(1, 6), make_array({1, 3}));  
}