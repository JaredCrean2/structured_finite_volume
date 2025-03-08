#include "gtest/gtest.h"
#include "disc/disc_block.h"
#include "disc/disc_interface.h"
#include "disc/discretization.h"
#include "mesh/adjacent_block_indexer.h"

using namespace structured_fv;

namespace {
class DiscTester : public ::testing::Test
{
  public:
    DiscTester()
    {
      mesh::MeshSpec spec(2, 1, m_num_bc_ghost_cells);
      spec.blocks(0, 0) = mesh::MeshBlockSpec(3, 4);
      spec.blocks(1, 0) = mesh::MeshBlockSpec(5, 4);
      m_mesh = std::make_shared<mesh::StructuredMesh>(spec);
      m_disc = std::make_shared<disc::StructuredDisc>(m_mesh, m_num_bc_ghost_cells);
    }

    int m_num_bc_ghost_cells = 2;
    std::shared_ptr<mesh::StructuredMesh> m_mesh;
    std::shared_ptr<disc::StructuredDisc> m_disc;
};

}

TEST_F(DiscTester, Counts)
{
  EXPECT_EQ(m_disc->getNumBlocks(), 8);
  EXPECT_EQ(m_disc->getNumRegularBlocks(), 2);
  EXPECT_EQ(m_disc->getNumGhostBCBlocks(), 6);

  EXPECT_EQ(m_disc->getNumBlockInterfaces(), 7);
  EXPECT_EQ(m_disc->getNumRegularBlockInterfaces(), 1);
  EXPECT_EQ(m_disc->getNumGhostBCBlockInterfaces(), 6);
}

TEST_F(DiscTester, ElementFieldUpdator)
{
  int nvals_per_element = 2;
  disc::ElementField<int> field(*m_disc, nvals_per_element);

  auto vals = [](UInt block, UInt i, UInt j, UInt k)
  {
    return block + 2*i + 3*j + 4*k;
  };

  for (UInt b=0; b < field.getNumBlocks(); ++b)
  {
    auto block = m_disc->getBlock(b);
    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        for (int k=0; k < nvals_per_element; ++k)
          field(b, i, j, k) = vals(b, i, j, k);
  }

  field.updateGhostValues();

  {
    UInt block_id = 0;
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    const disc::StructuredBlockInterface& north_iface = m_disc->getBlockInterface(1);
    const disc::StructuredBlockInterface& south_iface = m_disc->getBlockInterface(4);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        for (int k=0; k < nvals_per_element; ++k)
          EXPECT_EQ(field(block_id, i, j, k), vals(block_id, i, j, k));

    for (UInt i : block.getOwnedCells().getXRange())
    {
      for (UInt k=0; k < nvals_per_element; ++k)
      {
        UInt j1 = *block.getOwnedCells().getYRange().end();
        UInt j2 = j1 + 1;        
        EXPECT_EQ(field(block_id, i, j1, k), field(north_iface.getBlockIdR(), i-2, 2, k));
        EXPECT_EQ(field(block_id, i, j1, k), vals(north_iface.getBlockIdR(), i-2, 2, k));

        EXPECT_EQ(field(block_id, i, j2, k), field(north_iface.getBlockIdR(), i-2, 3, k));
        EXPECT_EQ(field(block_id, i, j2, k), vals(north_iface.getBlockIdR(), i-2, 3, k));
      
        j1 = 0;
        j2 = 1;
        EXPECT_EQ(field(block_id, i, j1, k), field(south_iface.getBlockIdR(), i-2, 0, k));
        EXPECT_EQ(field(block_id, i, j1, k), vals(south_iface.getBlockIdR(), i-2, 0, k));

        EXPECT_EQ(field(block_id, i, j2, k), field(south_iface.getBlockIdR(), i-2, 1, k));
        EXPECT_EQ(field(block_id, i, j2, k), vals(south_iface.getBlockIdR(), i-2, 1, k));
      }
    }

    const disc::StructuredBlockInterface& east_iface = m_disc->getBlockInterface(0);
    const disc::StructuredBlockInterface& west_iface = m_disc->getBlockInterface(6);

    for (UInt j : block.getOwnedCells().getYRange())
    {
      for (UInt k=0; k < nvals_per_element; ++k)
      {
        UInt i1 = *block.getOwnedCells().getXRange().end();
        UInt i2 = i1 + 1;
        EXPECT_EQ(field(block_id, i1, j, 0), field(east_iface.getBlockIdR(), 2, j, 0));
        EXPECT_EQ(field(block_id, i1, j, 0), vals(east_iface.getBlockIdR(), 2, j, 0));

        EXPECT_EQ(field(block_id, i2, j, k), field(east_iface.getBlockIdR(), 3, j, k));
        EXPECT_EQ(field(block_id, i2, j, k), field(east_iface.getBlockIdR(), 3, j, k));


        i1 = 0;
        i2 = 1;
        EXPECT_EQ(field(block_id, i1, j, k), field(west_iface.getBlockIdR(), 0, j-2, k));
        EXPECT_EQ(field(block_id, i1, j, k), vals(west_iface.getBlockIdR(), 0, j-2, k));

        EXPECT_EQ(field(block_id, i2, j, k), field(west_iface.getBlockIdR(), 1, j-2, k));
        EXPECT_EQ(field(block_id, i2, j, k), vals(west_iface.getBlockIdR(), 1, j-2, k));
      }
    }
  }

  {
    UInt block_id = 1;
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    const disc::StructuredBlockInterface& north_iface = m_disc->getBlockInterface(2);
    const disc::StructuredBlockInterface& south_iface = m_disc->getBlockInterface(5);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        for (int k=0; k < nvals_per_element; ++k)
          EXPECT_EQ(field(block_id, i, j, k), vals(block_id, i, j, k));

    for (UInt i : block.getOwnedCells().getXRange())
    {
      for (UInt k=0; k < nvals_per_element; ++k)
      {
        UInt j1 = *block.getOwnedCells().getYRange().end();
        UInt j2 = j1 + 1;        
        EXPECT_EQ(field(block_id, i, j1, k), field(north_iface.getBlockIdR(), i-2, 2, k));
        EXPECT_EQ(field(block_id, i, j1, k), vals(north_iface.getBlockIdR(), i-2, 2, k));

        EXPECT_EQ(field(block_id, i, j2, k), field(north_iface.getBlockIdR(), i-2, 3, k));
        EXPECT_EQ(field(block_id, i, j2, k), vals(north_iface.getBlockIdR(), i-2, 3, k));
      
        j1 = 0;
        j2 = 1;
        EXPECT_EQ(field(block_id, i, j1, k), field(south_iface.getBlockIdR(), i-2, 0, k));
        EXPECT_EQ(field(block_id, i, j1, k), vals(south_iface.getBlockIdR(), i-2, 0, k));

        EXPECT_EQ(field(block_id, i, j2, k), field(south_iface.getBlockIdR(), i-2, 1, k));
        EXPECT_EQ(field(block_id, i, j2, k), vals(south_iface.getBlockIdR(), i-2, 1, k));
      }
    }

    const disc::StructuredBlockInterface& east_iface = m_disc->getBlockInterface(3);
    const disc::StructuredBlockInterface& west_iface = m_disc->getBlockInterface(0);

    for (UInt j : block.getOwnedCells().getYRange())
    {
      for (UInt k=0; k < nvals_per_element; ++k)
      {
        UInt i1 = *block.getOwnedCells().getXRange().end();
        UInt i2 = i1 + 1;
        EXPECT_EQ(field(block_id, i1, j, k), field(east_iface.getBlockIdR(), 2, j-2, k));
        EXPECT_EQ(field(block_id, i1, j, k), vals(east_iface.getBlockIdR(), 2, j-2, k));

        EXPECT_EQ(field(block_id, i2, j, k), field(east_iface.getBlockIdR(), 3, j-2, k));
        EXPECT_EQ(field(block_id, i2, j, k), field(east_iface.getBlockIdR(), 3, j-2, k));


        i1 = 0;
        i2 = 1;
        EXPECT_EQ(field(block_id, i1, j, k), field(west_iface.getBlockIdL(), 3, j, k));
        EXPECT_EQ(field(block_id, i1, j, k), vals(west_iface.getBlockIdL(), 3, j, k));

        EXPECT_EQ(field(block_id, i2, j, k), field(west_iface.getBlockIdL(), 4, j, k));
        EXPECT_EQ(field(block_id, i2, j, k), vals(west_iface.getBlockIdL(), 4, j, k));
      }
    }
  }  


  {
    UInt block_id = 2;
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    const disc::StructuredBlockInterface& south_iface = m_disc->getBlockInterface(1);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        for (int k=0; k < nvals_per_element; ++k)
          EXPECT_EQ(field(block_id, i, j, k), vals(block_id, i, j, k));

    for (UInt i : block.getOwnedCells().getXRange())
    {
      for (UInt k=0; k < nvals_per_element; ++k)
      {
        UInt j1 = 0;
        UInt j2 = 1;
        EXPECT_EQ(field(block_id, i, j1, k), field(south_iface.getBlockIdL(), i+2, 4, k));
        EXPECT_EQ(field(block_id, i, j1, k), vals(south_iface.getBlockIdL(), i+2, 4, k));

        EXPECT_EQ(field(block_id, i, j2, k), field(south_iface.getBlockIdL(), i+2, 5, k));
        EXPECT_EQ(field(block_id, i, j2, k), vals(south_iface.getBlockIdL(), i+2, 5, k));
      }
    }
  }
}

TEST_F(DiscTester, VertFieldUpdator)
{
  int nvals_per_element = 2;
  disc::VertexField<int> field(*m_disc, nvals_per_element);

  auto vals = [](UInt block, UInt i, UInt j, UInt k)
  {
    return block + 2*i + 3*j + 4*k;
  };

  for (UInt b=0; b < field.getNumBlocks(); ++b)
  {
    auto block = m_disc->getBlock(b);
    for (UInt i : block.getOwnedVerts().getXRange())
      for (UInt j : block.getOwnedVerts().getYRange())
        for (int k=0; k < nvals_per_element; ++k)
          field(b, i, j, k) = vals(b, i, j, k);
  }

  field.updateGhostValues();

  {
    UInt block_id = 0;
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    const disc::StructuredBlockInterface& north_iface = m_disc->getBlockInterface(1);
    const disc::StructuredBlockInterface& south_iface = m_disc->getBlockInterface(4);

    for (UInt i : block.getOwnedVerts().getXRange())
      for (UInt j : block.getOwnedVerts().getYRange())
        for (int k=0; k < nvals_per_element; ++k)
          EXPECT_EQ(field(block_id, i, j, k), vals(block_id, i, j, k));

    for (UInt i : block.getOwnedVerts().getXRange())
    {
      for (UInt k=0; k < nvals_per_element; ++k)
      {
        UInt j1 = *block.getOwnedVerts().getYRange().end();
        UInt j2 = j1 + 1;        
        EXPECT_EQ(field(block_id, i, j1, k), field(north_iface.getBlockIdR(), i-2, 3, k));
        EXPECT_EQ(field(block_id, i, j1, k), vals(north_iface.getBlockIdR(), i-2, 3, k));

        EXPECT_EQ(field(block_id, i, j2, k), field(north_iface.getBlockIdR(), i-2, 4, k));
        EXPECT_EQ(field(block_id, i, j2, k), vals(north_iface.getBlockIdR(), i-2, 4, k));
      
        j1 = 0;
        j2 = 1;
        EXPECT_EQ(field(block_id, i, j1, k), field(south_iface.getBlockIdR(), i-2, 0, k));
        EXPECT_EQ(field(block_id, i, j1, k), vals(south_iface.getBlockIdR(), i-2, 0, k));

        EXPECT_EQ(field(block_id, i, j2, k), field(south_iface.getBlockIdR(), i-2, 1, k));
        EXPECT_EQ(field(block_id, i, j2, k), vals(south_iface.getBlockIdR(), i-2, 1, k));
      }
    }

    const disc::StructuredBlockInterface& east_iface = m_disc->getBlockInterface(0);
    const disc::StructuredBlockInterface& west_iface = m_disc->getBlockInterface(6);

    for (UInt j : block.getOwnedVerts().getYRange())
    {
      for (UInt k=0; k < nvals_per_element; ++k)
      {
        UInt i1 = *block.getOwnedVerts().getXRange().end();
        UInt i2 = i1 + 1;
        EXPECT_EQ(field(block_id, i1, j, 0), field(east_iface.getBlockIdR(), 3, j, 0));  //TODO: why 0
        EXPECT_EQ(field(block_id, i1, j, 0), vals(east_iface.getBlockIdR(), 3, j, 0));

        EXPECT_EQ(field(block_id, i2, j, k), field(east_iface.getBlockIdR(), 4, j, k));
        EXPECT_EQ(field(block_id, i2, j, k), field(east_iface.getBlockIdR(), 4, j, k));


        i1 = 0;
        i2 = 1;
        EXPECT_EQ(field(block_id, i1, j, k), field(west_iface.getBlockIdR(), 0, j-2, k));
        EXPECT_EQ(field(block_id, i1, j, k), vals(west_iface.getBlockIdR(), 0, j-2, k));

        EXPECT_EQ(field(block_id, i2, j, k), field(west_iface.getBlockIdR(), 1, j-2, k));
        EXPECT_EQ(field(block_id, i2, j, k), vals(west_iface.getBlockIdR(), 1, j-2, k));
      }
    }
  }

  {
    UInt block_id = 1;
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    const disc::StructuredBlockInterface& north_iface = m_disc->getBlockInterface(2);
    const disc::StructuredBlockInterface& south_iface = m_disc->getBlockInterface(5);

    for (UInt i : block.getOwnedVerts().getXRange())
      for (UInt j : block.getOwnedVerts().getYRange())
        for (int k=0; k < nvals_per_element; ++k)
          EXPECT_EQ(field(block_id, i, j, k), vals(block_id, i, j, k));

    for (UInt i : block.getOwnedVerts().getXRange())
    {
      for (UInt k=0; k < nvals_per_element; ++k)
      {
        UInt j1 = *block.getOwnedVerts().getYRange().end();
        UInt j2 = j1 + 1;        
        EXPECT_EQ(field(block_id, i, j1, k), field(north_iface.getBlockIdR(), i-2, 3, k));
        EXPECT_EQ(field(block_id, i, j1, k), vals(north_iface.getBlockIdR(), i-2, 3, k));

        EXPECT_EQ(field(block_id, i, j2, k), field(north_iface.getBlockIdR(), i-2, 4, k));
        EXPECT_EQ(field(block_id, i, j2, k), vals(north_iface.getBlockIdR(), i-2, 4, k));
      
        j1 = 0;
        j2 = 1;
        EXPECT_EQ(field(block_id, i, j1, k), field(south_iface.getBlockIdR(), i-2, 0, k));
        EXPECT_EQ(field(block_id, i, j1, k), vals(south_iface.getBlockIdR(), i-2, 0, k));

        EXPECT_EQ(field(block_id, i, j2, k), field(south_iface.getBlockIdR(), i-2, 1, k));
        EXPECT_EQ(field(block_id, i, j2, k), vals(south_iface.getBlockIdR(), i-2, 1, k));
      }
    }

    const disc::StructuredBlockInterface& east_iface = m_disc->getBlockInterface(3);
    const disc::StructuredBlockInterface& west_iface = m_disc->getBlockInterface(0);

    for (UInt j : block.getOwnedVerts().getYRange())
    {
      for (UInt k=0; k < nvals_per_element; ++k)
      {
        UInt i1 = *block.getOwnedVerts().getXRange().end();
        UInt i2 = i1 + 1;
        EXPECT_EQ(field(block_id, i1, j, k), field(east_iface.getBlockIdR(), 3, j-2, k));
        EXPECT_EQ(field(block_id, i1, j, k), vals(east_iface.getBlockIdR(), 3, j-2, k));

        EXPECT_EQ(field(block_id, i2, j, k), field(east_iface.getBlockIdR(), 4, j-2, k));
        EXPECT_EQ(field(block_id, i2, j, k), field(east_iface.getBlockIdR(), 4, j-2, k));


        i1 = 0;
        i2 = 1;
        EXPECT_EQ(field(block_id, i1, j, k), field(west_iface.getBlockIdL(), 3, j, k));
        EXPECT_EQ(field(block_id, i1, j, k), vals(west_iface.getBlockIdL(), 3, j, k));

        EXPECT_EQ(field(block_id, i2, j, k), field(west_iface.getBlockIdL(), 4, j, k));
        EXPECT_EQ(field(block_id, i2, j, k), vals(west_iface.getBlockIdL(), 4, j, k));
      }
    }
  }  


  {
    UInt block_id = 2;
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    const disc::StructuredBlockInterface& south_iface = m_disc->getBlockInterface(1);

    for (UInt i : block.getOwnedVerts().getXRange())
      for (UInt j : block.getOwnedVerts().getYRange())
        for (int k=0; k < nvals_per_element; ++k)
          EXPECT_EQ(field(block_id, i, j, k), vals(block_id, i, j, k));

    for (UInt i : block.getOwnedVerts().getXRange())
    {
      for (UInt k=0; k < nvals_per_element; ++k)
      {
        UInt j1 = 0;
        UInt j2 = 1;
        EXPECT_EQ(field(block_id, i, j1, k), field(south_iface.getBlockIdL(), i+2, 4, k));
        EXPECT_EQ(field(block_id, i, j1, k), vals(south_iface.getBlockIdL(), i+2, 4, k));

        EXPECT_EQ(field(block_id, i, j2, k), field(south_iface.getBlockIdL(), i+2, 5, k));
        EXPECT_EQ(field(block_id, i, j2, k), vals(south_iface.getBlockIdL(), i+2, 5, k));
      }
    }
  }
}