#include "gtest/gtest.h"
#include "disc/disc_block.h"
#include "disc/disc_interface.h"
#include "disc/discretization.h"
#include "disc/elem_field.h"
#include "disc/face_field.h"
#include "disc/vert_field.h"
#include "mesh/adjacent_block_indexer.h"
#include "utils/face_iter_per_direction.h"

using namespace structured_fv;

namespace {
class DiscTester : public ::testing::Test
{
  public:
    DiscTester() :
      spec(2, 1, m_num_bc_ghost_cells)
    {
      //spec = mesh::MeshSpec(2, 1, m_num_bc_ghost_cells);
      spec.blocks(0, 0) = mesh::MeshBlockSpec(3, 4, 0, [](Real x, Real y) { return std::array<Real, 2>{x, y}; });
      spec.blocks(1, 0) = mesh::MeshBlockSpec(5, 4, 0, [](Real x, Real y) { return std::array<Real, 2>{x+1, y}; });
      m_mesh = std::make_shared<mesh::StructuredMesh>(spec);
      m_disc = std::make_shared<disc::StructuredDisc>(m_mesh, m_num_bc_ghost_cells);
    }

    int m_num_bc_ghost_cells = 2;
    mesh::MeshSpec spec;
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

TEST_F(DiscTester, CoordField)
{
  const disc::VertexField<Real>& coord_field = *m_disc->getCoordField();
  {
    UInt block_id = 0;
    Real dx = 1.0/3;
    Real dy = 1.0/4;
    Real x0 = -2*dx;
    Real y0 = -2*dy;
    Range2D verts(0, 5, 0, 9);
    for (UInt i : verts.getXRange())
      for (UInt j : verts.getYRange())
      {
        if ((i < 2 && j < 2) || (i < 2 && j > 6))
          continue;

        EXPECT_NEAR(coord_field(block_id, i, j, 0), x0 + i*dx, 1e-13);
        EXPECT_NEAR(coord_field(block_id, i, j, 1), y0 + j*dy, 1e-13);
      }

    x0 = 1;
    y0 = 0;
    dx = 1.0/5;
    dy = 1.0/4;
    verts = Range2D(5, 8, 2, 7);

    for (UInt i : verts.getXRange())
      for (UInt j : verts.getYRange())
      {
        EXPECT_NEAR(coord_field(block_id, i, j, 0), x0 + (i - 5)*dx, 1e-13);
        EXPECT_NEAR(coord_field(block_id, i, j, 1), y0 + (j - 2)*dy, 1e-13);
      }    
  }

  {
    UInt block_id = 1;
    Real dx = 1.0/3;
    Real dy = 1.0/4;
    Real x0 = 1 - 2*dx;
    Real y0 = 0;
    Range2D verts(0, 3, 2, 7);
    for (UInt i : verts.getXRange())
      for (UInt j : verts.getYRange())
      {
        EXPECT_NEAR(coord_field(block_id, i, j, 0), x0 + i*dx, 1e-13);
        EXPECT_NEAR(coord_field(block_id, i, j, 1), y0 + (j-2)*dy, 1e-13);
      }

    dx = 1.0/5;
    dy = 1.0/4;
    x0 = 1;
    y0 = -2*dy;
    verts = Range2D(2, 9, 0, 9);
    for (UInt i : verts.getXRange())
      for (UInt j : verts.getYRange())
      {
        if (i > 7 && (j < 2 || j > 6))
          continue;
        
        EXPECT_NEAR(coord_field(block_id, i, j, 0), x0 + (i - 2)*dx, 1e-13);
        EXPECT_NEAR(coord_field(block_id, i, j, 1), y0 + j*dy, 1e-13);
      }

  }
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

TEST_F(DiscTester, FaceFieldSetConstant)
{
  int val = 42;
  disc::FaceField<int> field(*m_disc, 1);
  field.set(val);

  for (UInt block_id=0; block_id < m_disc->getNumBlocks(); ++block_id)
  {
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    FaceRangePerDirection faces = block.getOwnedAndGhostFacesWithCorners();
    auto east_data = field.getData(block_id, mesh::NeighborDirection::East);
    auto north_data = field.getData(block_id, mesh::NeighborDirection::North);

    for (UInt i : faces.getXRange(XDirTag()))
      for (UInt j : faces.getYRange(XDirTag()))
        EXPECT_EQ(east_data(i, j, 0), val);

    for (UInt i : faces.getXRange(YDirTag()))
      for (UInt j : faces.getYRange(YDirTag()))
        EXPECT_EQ(north_data(i, j, 0), val);      
  }
}

TEST_F(DiscTester, FaceFieldSetField)
{
  disc::FaceField<double> field(*m_disc, 2);
  auto f = [](Real x, Real y) { return std::array<double, 2>{x, y}; };
  field.set(f);

  std::array<double, 2> dx = {1.0/spec.blocks(0, 0).num_cells_x, 1.0/spec.blocks(1, 0).num_cells_x};
  std::array<double, 2> dy = {1.0/spec.blocks(0, 0).num_cells_y, 1.0/spec.blocks(1, 0).num_cells_y};
  std::array<double, 2> x0 = {0, 1};
  std::array<double, 2> y0 = {0, 0};

  for (UInt block_id=0; block_id < m_disc->getNumRegularBlocks(); ++block_id)
  {
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    FaceRangePerDirection faces = block.getOwnedFaces();
    auto east_data = field.getData(block_id, mesh::NeighborDirection::East);
    auto north_data = field.getData(block_id, mesh::NeighborDirection::North);

    for (UInt i : faces.getXRange(XDirTag()))
      for (UInt j : faces.getYRange(XDirTag()))
      {
        EXPECT_NEAR(east_data(i, j, 0), x0[block_id] + (i-2)*dx[block_id], 1e-13);
        EXPECT_NEAR(east_data(i, j, 1), y0[block_id] + (j-2 + 0.5)*dy[block_id], 1e-13);
      }

    for (UInt i : faces.getXRange(YDirTag()))
      for (UInt j : faces.getYRange(YDirTag()))
      {
        EXPECT_NEAR(north_data(i, j, 0), x0[block_id] + (i-2 + 0.5)*dx[block_id], 1e-13);
        EXPECT_NEAR(north_data(i, j, 1), y0[block_id] + (j-2)*dy[block_id], 1e-13);
      }      
  }
}