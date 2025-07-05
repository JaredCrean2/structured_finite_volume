#include "gtest/gtest.h"
#include "disc/disc_vector.h"
#include "mesh/adjacent_block_indexer.h"
#include "physics/physics_model.h"

using namespace structured_fv;

namespace {
class PhysicsModelTester : public ::testing::Test
{
  public:
    PhysicsModelTester() :
      spec(2, 1, m_num_bc_ghost_cells)
    {
      //spec = mesh::MeshSpec(2, 1, m_num_bc_ghost_cells);
      spec.blocks(0, 0) = mesh::MeshBlockSpec(3, 4, 0, [](Real x, Real y) { return FixedVec<Real, 2>{x, y}; });
      spec.blocks(1, 0) = mesh::MeshBlockSpec(5, 4, 0, [](Real x, Real y) { return FixedVec<Real, 2>{x+1, y}; });
      m_mesh = std::make_shared<mesh::StructuredMesh>(spec);
      m_disc = std::make_shared<disc::StructuredDisc>(m_mesh, m_num_bc_ghost_cells, m_dofs_per_cell);
    }

    UInt m_dofs_per_cell = 4;
    int m_num_bc_ghost_cells = 2;
    mesh::MeshSpec spec;
    std::shared_ptr<mesh::StructuredMesh> m_mesh;
    std::shared_ptr<disc::StructuredDisc> m_disc;
};

}

TEST_F(PhysicsModelTester, VecToField)
{
  auto vec = std::make_shared<disc::DiscVector<Real>>(m_disc, "vec");
  auto field = std::make_shared<disc::ElementField<Real>>(*m_disc, m_dofs_per_cell);
  auto func =  [](Real x, Real y, UInt v) { return x + 2*y + v; };

  auto coord_field    = m_disc->getCoordField();
  auto dof_nums_field = m_disc->getDofNumbering();
  for (UInt block_id : m_disc->getRegularBlocksIds())
  {
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    const auto& coords   = coord_field->getData(block_id);
    const auto& dof_nums = dof_nums_field->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        for (UInt k=0; k < m_dofs_per_cell; ++k)
        {
          const auto [x, y] = disc::computeCellCentroid(coords, i, j);
          (*vec)(dof_nums(i, j, k)) = func(x, y, k);
        }
  }


  vecToField(m_disc, vec, field);

  for (UInt block_id : m_disc->getRegularBlocksIds())
  {
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    const auto& field_vals = field->getData(block_id);
    const auto& coords   = coord_field->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        for (UInt k=0; k < m_dofs_per_cell; ++k)
        {
          const auto [x, y] = disc::computeCellCentroid(coords, i, j);
          EXPECT_EQ(field_vals(i, j, k), func(x, y, k));
        }
  }

  const disc::StructuredBlockInterface& iface = m_disc->getBlockInterface(0);
  const mesh::AdjacentBlockIndexer& indexerL = iface.getAdjBlockCellIndexerL();
  const mesh::AdjacentBlockIndexer& indexerR = iface.getAdjBlockCellIndexerR();
  const auto& field_valsL = field->getData(iface.getBlockIdL());
  const auto& field_valsR = field->getData(iface.getBlockIdR());
  const auto& coordsL = coord_field->getData(iface.getBlockIdL());
  const auto& coordsR = coord_field->getData(iface.getBlockIdR());

  for (UInt i : iface.getOwnedBoundaryCellsL().getXRange())
    for (UInt j : iface.getOwnedBoundaryCellsL().getYRange())
      for (int k=0; k < m_num_bc_ghost_cells; ++k)
      {
        const auto [iprime, jprime] = computeIndices(iface.getNeighborDirectionL(), k+1, i, j);
        const auto [ineighbor, jneighbor] = indexerL(iprime, jprime);
        const auto [x, y] = disc::computeCellCentroid(coordsL, iprime, jprime);
        const auto [x2, y2] = disc::computeCellCentroid(coordsR, ineighbor, jneighbor);

        EXPECT_EQ(x, x2);
        EXPECT_EQ(y, y2);

        for (UInt v=0; v < m_dofs_per_cell; ++v)
        {
          EXPECT_EQ(field_valsL(iprime, jprime, v), func(x, y, v));
          EXPECT_EQ(field_valsR(ineighbor, jneighbor, v), func(x, y, v));
        }
      }

  for (UInt i : iface.getOwnedBoundaryCellsR().getXRange())
    for (UInt j : iface.getOwnedBoundaryCellsR().getYRange())
      for (int k=0; k < m_num_bc_ghost_cells; ++k)
      {
        const auto [iprime, jprime] = computeIndices(iface.getNeighborDirectionR(), k+1, i, j);
        const auto [ineighbor, jneighbor] = indexerR(iprime, jprime);
        const auto [x, y] = disc::computeCellCentroid(coordsR, iprime, jprime);
        const auto [x2, y2] = disc::computeCellCentroid(coordsL, ineighbor, jneighbor);

        EXPECT_EQ(x, x2);
        EXPECT_EQ(y, y2);

        for (UInt v=0; v < m_dofs_per_cell; ++v)
        {
          EXPECT_EQ(field_valsR(iprime, jprime, v), func(x, y, v));
          EXPECT_EQ(field_valsL(ineighbor, jneighbor, v), func(x, y, v));
        }
      }
}

TEST_F(PhysicsModelTester, VecToFieldDot)
{
  Real h = 1e-40;
  auto vec = std::make_shared<disc::DiscVector<Real>>(m_disc, "vec");
  auto vec_dot = std::make_shared<disc::DiscVector<Real>>(m_disc, "vec_dot");

  auto field = std::make_shared<disc::ElementField<Complex>>(*m_disc, m_dofs_per_cell);
  auto func_real = [](Real x, Real y, UInt v) { return x + 2*y + v; };
  auto func_imag = [](Real x, Real y, UInt v) { return 2*(x + 2*y + v); };

  auto coord_field    = m_disc->getCoordField();
  auto dof_nums_field = m_disc->getDofNumbering();
  for (UInt block_id : m_disc->getRegularBlocksIds())
  {
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    const auto& coords   = coord_field->getData(block_id);
    const auto& dof_nums = dof_nums_field->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        for (UInt k=0; k < m_dofs_per_cell; ++k)
        {
          const auto [x, y] = disc::computeCellCentroid(coords, i, j);
          (*vec)(dof_nums(i, j, k)) = func_real(x, y, k);
          (*vec_dot)(dof_nums(i, j, k)) = func_imag(x, y, k);
        }
  }


  vecToFieldDot(m_disc, vec, vec_dot, field, h);

  for (UInt block_id : m_disc->getRegularBlocksIds())
  {
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    const auto& field_vals = field->getData(block_id);
    const auto& coords   = coord_field->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        for (UInt k=0; k < m_dofs_per_cell; ++k)
        {
          const auto [x, y] = disc::computeCellCentroid(coords, i, j);
          EXPECT_EQ(field_vals(i, j, k), Complex(func_real(x, y, k), h*func_imag(x, y, k)));
        }
  }

  const disc::StructuredBlockInterface& iface = m_disc->getBlockInterface(0);
  const mesh::AdjacentBlockIndexer& indexerL = iface.getAdjBlockCellIndexerL();
  const mesh::AdjacentBlockIndexer& indexerR = iface.getAdjBlockCellIndexerR();
  const auto& field_valsL = field->getData(iface.getBlockIdL());
  const auto& field_valsR = field->getData(iface.getBlockIdR());
  const auto& coordsL = coord_field->getData(iface.getBlockIdL());
  const auto& coordsR = coord_field->getData(iface.getBlockIdR());

  for (UInt i : iface.getOwnedBoundaryCellsL().getXRange())
    for (UInt j : iface.getOwnedBoundaryCellsL().getYRange())
      for (int k=0; k < m_num_bc_ghost_cells; ++k)
      {
        const auto [iprime, jprime] = computeIndices(iface.getNeighborDirectionL(), k+1, i, j);
        const auto [ineighbor, jneighbor] = indexerL(iprime, jprime);
        const auto [x, y] = disc::computeCellCentroid(coordsL, iprime, jprime);
        const auto [x2, y2] = disc::computeCellCentroid(coordsR, ineighbor, jneighbor);

        EXPECT_EQ(x, x2);
        EXPECT_EQ(y, y2);

        for (UInt v=0; v < m_dofs_per_cell; ++v)
        {
          EXPECT_EQ(field_valsL(iprime, jprime, v), Complex(func_real(x, y, v), h*func_imag(x, y, v)));
          EXPECT_EQ(field_valsR(ineighbor, jneighbor, v), Complex(func_real(x, y, v), h*func_imag(x, y, v)));
        }
      }

  for (UInt i : iface.getOwnedBoundaryCellsR().getXRange())
    for (UInt j : iface.getOwnedBoundaryCellsR().getYRange())
      for (int k=0; k < m_num_bc_ghost_cells; ++k)
      {
        const auto [iprime, jprime] = computeIndices(iface.getNeighborDirectionR(), k+1, i, j);
        const auto [ineighbor, jneighbor] = indexerR(iprime, jprime);
        const auto [x, y] = disc::computeCellCentroid(coordsR, iprime, jprime);
        const auto [x2, y2] = disc::computeCellCentroid(coordsL, ineighbor, jneighbor);

        EXPECT_EQ(x, x2);
        EXPECT_EQ(y, y2);

        for (UInt v=0; v < m_dofs_per_cell; ++v)
        {
          EXPECT_EQ(field_valsR(iprime, jprime, v), Complex(func_real(x, y, v), h*func_imag(x, y, v)));
          EXPECT_EQ(field_valsL(ineighbor, jneighbor, v), Complex(func_real(x, y, v), h*func_imag(x, y, v)));
        }
      }
}


TEST_F(PhysicsModelTester, FieldToVec)
{
  auto vec = std::make_shared<disc::DiscVector<Real>>(m_disc, "vec");
  auto field = std::make_shared<disc::ElementField<Real>>(*m_disc, m_dofs_per_cell);
  auto func =  [](Real x, Real y, UInt v) { return x + 2*y + v; };

  auto coord_field    = m_disc->getCoordField();
  auto dof_nums_field = m_disc->getDofNumbering();
  for (UInt block_id : m_disc->getRegularBlocksIds())
  {
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    const auto& coords   = coord_field->getData(block_id);
    const auto& field_vals = field->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        for (UInt k=0; k < m_dofs_per_cell; ++k)
        {
          const auto [x, y] = disc::computeCellCentroid(coords, i, j);
          field_vals(i, j, k) =  func(x, y, k);;
        }
  }


  fieldToVec(m_disc, field, vec);

  for (UInt block_id : m_disc->getRegularBlocksIds())
  {
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    const auto& coords   = coord_field->getData(block_id);
    const auto& dof_nums = dof_nums_field->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        for (UInt k=0; k < m_dofs_per_cell; ++k)
        {
          const auto [x, y] = disc::computeCellCentroid(coords, i, j);
          EXPECT_EQ((*vec)(dof_nums(i, j, k)), func(x, y, k));
        }
  }
}

TEST_F(PhysicsModelTester, FieldToVecDot)
{
  Real h = 1e-40;
  auto vec = std::make_shared<disc::DiscVector<Real>>(m_disc, "vec");
  auto vec_dot = std::make_shared<disc::DiscVector<Real>>(m_disc, "vec_dot");
  auto field = std::make_shared<disc::ElementField<Complex>>(*m_disc, m_dofs_per_cell);
  auto func_real =  [](Real x, Real y, UInt v) { return x + 2*y + v; };
  auto func_imag =  [](Real x, Real y, UInt v) { return 2*(x + 2*y + v); };


  auto coord_field    = m_disc->getCoordField();
  auto dof_nums_field = m_disc->getDofNumbering();
  for (UInt block_id : m_disc->getRegularBlocksIds())
  {
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    const auto& coords   = coord_field->getData(block_id);
    const auto& field_vals = field->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        for (UInt k=0; k < m_dofs_per_cell; ++k)
        {
          const auto [x, y] = disc::computeCellCentroid(coords, i, j);
          field_vals(i, j, k) =  Complex(func_real(x, y, k), h*func_imag(x, y, k));
        }
  }


  fieldToVecDot(m_disc, field, vec_dot, h);

  for (UInt block_id : m_disc->getRegularBlocksIds())
  {
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    const auto& coords   = coord_field->getData(block_id);
    const auto& dof_nums = dof_nums_field->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        for (UInt k=0; k < m_dofs_per_cell; ++k)
        {
          const auto [x, y] = disc::computeCellCentroid(coords, i, j);
          EXPECT_DOUBLE_EQ((*vec_dot)(dof_nums(i, j, k)), func_imag(x, y, k));
        }
  }
}