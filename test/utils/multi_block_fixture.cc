#include "multi_block_fixture.h"
#include "rotation_utils.h"
#include "gtest/gtest.h"

namespace test_utils {

MultiBlockFixture::MultiBlockFixture()
    : spec1(1, 1, m_num_bc_ghost_cells), spec2(2, 2, m_num_bc_ghost_cells) {
  spec1.blocks(0, 0) = mesh::MeshBlockSpec(6, 8, 0, [](Real x, Real y) {
    return FixedVec<Real, 2>{2 * x, 2 * y};
  });
  m_mesh1 = std::make_shared<mesh::StructuredMesh>(spec1);
  m_disc1 = std::make_shared<disc::StructuredDisc>(
      m_mesh1, m_num_bc_ghost_cells, m_dofs_per_cell);

  auto f1 = [](Real x, Real y) { return FixedVec<Real, 2>{x, y}; };
  auto f2 = [](Real x, Real y) { return FixedVec<Real, 2>{x + 1, y}; };
  auto f3 = [](Real x, Real y) { return FixedVec<Real, 2>{x, y + 1}; };
  auto f4 = [](Real x, Real y) { return FixedVec<Real, 2>{x + 1, y + 1}; };
  spec2.blocks(0, 0) = test_utils::rotate_block(mesh::MeshBlockSpec(3, 4, 0, f1), 0);
  spec2.blocks(1, 0) = test_utils::rotate_block(mesh::MeshBlockSpec(3, 4, 0, f2), 1);
  spec2.blocks(0, 1) = test_utils::rotate_block(mesh::MeshBlockSpec(3, 4, 0, f3), 2);
  spec2.blocks(1, 1) = test_utils::rotate_block(mesh::MeshBlockSpec(3, 4, 0, f4), 3);

  m_mesh1 = std::make_shared<mesh::StructuredMesh>(spec1);
  m_disc1 = std::make_shared<disc::StructuredDisc>(
      m_mesh1, m_num_bc_ghost_cells, m_dofs_per_cell);

  m_mesh2 = std::make_shared<mesh::StructuredMesh>(spec2);
  m_disc2 = std::make_shared<disc::StructuredDisc>(
      m_mesh2, m_num_bc_ghost_cells, m_dofs_per_cell); 
}

void MultiBlockFixture::setup_disc(UInt dofs_per_cell)
{
  m_dofs_per_cell = dofs_per_cell;

  m_disc1 = std::make_shared<disc::StructuredDisc>(
      m_mesh1, m_num_bc_ghost_cells, m_dofs_per_cell);
  m_disc2 = std::make_shared<disc::StructuredDisc>(
      m_mesh2, m_num_bc_ghost_cells, m_dofs_per_cell);

  m_solution1 = std::make_shared<disc::DiscVector<Real>>(m_disc1, "solution");
  m_residual1 = std::make_shared<disc::DiscVector<Real>>(m_disc1, "residual");
  m_solution1->set(0);
  m_residual1->set(0);

  m_solution2 = std::make_shared<disc::DiscVector<Real>>(m_disc2, "solution");
  m_residual2 = std::make_shared<disc::DiscVector<Real>>(m_disc2, "residual");
  m_solution2->set(0);
  m_residual2->set(0);      
} 



// given the indices of a cell on the single block of mesh1, compute the
// (block, i, j) for the same cell on mesh 2
std::tuple<UInt, UInt, UInt> MultiBlockFixture::get_idx(UInt i, UInt j) const
{
  assert(i >= m_num_bc_ghost_cells && j >= m_num_bc_ghost_cells);
  UInt block, i2, j2, rotation;
  UInt nx = spec1.blocks(0, 0).num_cells_x;
  UInt ny = spec1.blocks(0, 0).num_cells_y;

  UInt ioffset = i - m_num_bc_ghost_cells;
  UInt joffset = j - m_num_bc_ghost_cells;
  bool first_col = ioffset < nx/2;
  bool first_row = joffset < ny/2;

  if (first_col && first_row)
  {
    block = 0;
    i2 = ioffset;
    j2 = joffset;
    rotation = spec2.blocks(0, 0).rotation;
  } else if (!first_col && first_row)
  {
    block = 1;
    i2 = ioffset - nx/2;
    j2 = joffset;
    rotation = spec2.blocks(1, 0).rotation;
  } else if (first_col && !first_row)
  {
    block = 2;
    i2 = ioffset;
    j2 = joffset - ny/2;
    rotation = spec2.blocks(0, 1).rotation;
  } else if (!first_col && !first_row)
  {
    block = 3;
    i2 = ioffset - nx/2;
    j2 = joffset - ny/2;
    rotation = spec2.blocks(1, 1).rotation;
  } else
  {
    throw std::runtime_error("unreachable");  // appease the compiler
  }

  auto [i3, j3] = test_utils::compute_indices(i2, j2, nx/2, ny/2, rotation);
  return {block, i3 + m_num_bc_ghost_cells, j3 + m_num_bc_ghost_cells};
}

void MultiBlockFixture::test_residuals_equal(
                          std::shared_ptr<disc::DiscVector<Real>> residual1,
                          std::shared_ptr<disc::DiscVector<Real>> residual2) const
{
  const disc::StructuredBlock &block1 = m_disc1->getBlock(0);
  structured_fv::disc::ElementFieldPtr<GlobalDof> dof_numbering1 = m_disc1->getDofNumbering();
  structured_fv::disc::ElementFieldPtr<GlobalDof> dof_numbering2 = m_disc2->getDofNumbering();
  for (UInt i : block1.getOwnedCells().getXRange())
    for (UInt j : block1.getOwnedCells().getYRange())
    {
      auto [block_id2, i2, j2] = get_idx(i, j);
      GlobalDof dof1 = (*dof_numbering1)(0, i, j, 0);
      GlobalDof dof2 = (*dof_numbering2)(block_id2, i2, j2, 0);

      const auto& coords1 = m_disc1->getCoordField()->getData(0);
      const auto& coords2 = m_disc2->getCoordField()->getData(block_id2);
      auto centroid1 = disc::computeCellCentroid(coords1, i, j);
      auto centroid2 = disc::computeCellCentroid(coords2, i2, j2);

      EXPECT_NEAR(centroid1[0], centroid2[0], 1e-13);
      EXPECT_NEAR(centroid1[1], centroid2[1], 1e-13);
      EXPECT_NEAR((*residual1)(dof1), (*residual2)(dof2), 1e-11);
    }
}

}