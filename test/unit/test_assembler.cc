#include "gtest/gtest.h"
#include "linear_system/assembler.h"
#include "linear_system/large_matrix.h"
#include "linear_system/large_matrix_dense.h"

using namespace structured_fv;

namespace {
using MatrixType = linear_system::LargeMatrixDense;
using AssemblerType = linear_system::Assembler<MatrixType>;

class AssemblerTester : public ::testing::Test
{
  public:
    AssemblerTester() :
      spec(1, 1, m_num_bc_ghost_cells)
    {
      spec.blocks(0, 0) = mesh::MeshBlockSpec(3, 4, 0, [](Real x, Real y) { return FixedVec<Real, 2>{x, y}; });
      m_mesh = std::make_shared<mesh::StructuredMesh>(spec);
      m_disc = std::make_shared<disc::StructuredDisc>(m_mesh, m_num_bc_ghost_cells, m_dofs_per_cell);
      linear_system::LargeMatrixOpts opts;
      m_matrix = std::make_shared<MatrixType>("A", m_disc->getNumDofs(), m_disc->getNumDofs(), opts);
      m_assembler = std::make_shared<AssemblerType>(m_disc, *m_matrix);
    }

    static constexpr UInt m_dofs_per_cell = 3;
    int m_num_bc_ghost_cells = 2;
    mesh::MeshSpec spec;
    std::shared_ptr<mesh::StructuredMesh> m_mesh;
    std::shared_ptr<disc::StructuredDisc> m_disc;
    std::shared_ptr<MatrixType> m_matrix;
    std::shared_ptr<AssemblerType> m_assembler;
};

}

TEST_F(AssemblerTester, AssembleValues)
{
  m_assembler->setBlock(0);
  FixedVec<linear_system::Indices, 1> row_indices = {{{2, 2}}};
  FixedVec<linear_system::Indices, 2> col_indices = {{{2, 2}, {2, 3}}};
  Matrix<Real, m_dofs_per_cell, m_dofs_per_cell*2> mat({1, 2, 3,
                                                        4, 5, 6,
                                                        7, 8, 9,
                                                        10, 11, 12,
                                                        13, 14, 15,
                                                        16, 17, 18});
  m_assembler->assembleValues(row_indices, col_indices, mat);

  const auto& dof_nums = m_disc->getDofNumbering()->getData(0);
  UInt row = 0;
  for (linear_system::Indices row_idx : row_indices)
    for (UInt k1=0; k1 < m_dofs_per_cell; ++k1)
    {
      UInt col = 0;
      UInt dof_row = dof_nums(row_idx.i, row_idx.j, k1);
      for (linear_system::Indices col_idx : col_indices)
        for (UInt k2=0; k2 < m_dofs_per_cell; ++k2)
        {
          UInt dof_col = dof_nums(col_idx.i, col_idx.j, k2);
          EXPECT_EQ((*m_matrix)(dof_row, dof_col), mat(row, col));
          col++;
        }

      row++;
    }
}

TEST_F(AssemblerTester, AssembleValuesWithAlpha)
{
  Real alpha = 2;
  m_assembler->setBlock(0);
  m_assembler->setAlpha(alpha);
  FixedVec<linear_system::Indices, 1> row_indices = {{{2, 2}}};
  FixedVec<linear_system::Indices, 2> col_indices = {{{2, 2}, {2, 3}}};
  Matrix<Real, m_dofs_per_cell, m_dofs_per_cell*2> mat({1, 2, 3,
                                                        4, 5, 6,
                                                        7, 8, 9,
                                                        10, 11, 12,
                                                        13, 14, 15,
                                                        16, 17, 18});

  Matrix<Real, m_dofs_per_cell, m_dofs_per_cell*2>  mat_orig = mat;                                                     
  m_assembler->assembleValues(row_indices, col_indices, mat);

  const auto& dof_nums = m_disc->getDofNumbering()->getData(0);
  UInt row = 0;
  for (linear_system::Indices row_idx : row_indices)
    for (UInt k1=0; k1 < m_dofs_per_cell; ++k1)
    {
      UInt col = 0;
      UInt dof_row = dof_nums(row_idx.i, row_idx.j, k1);
      for (linear_system::Indices col_idx : col_indices)
        for (UInt k2=0; k2 < m_dofs_per_cell; ++k2)
        {
          UInt dof_col = dof_nums(col_idx.i, col_idx.j, k2);
          EXPECT_EQ((*m_matrix)(dof_row, dof_col), alpha*mat_orig(row, col));
          col++;
        }

      row++;
    }
}

TEST_F(AssemblerTester, CopyAssembler)
{
  m_assembler->setBlock(0);
  FixedVec<linear_system::Indices, 1> row_indices = {{{2, 2}}};
  FixedVec<linear_system::Indices, 2> col_indices = {{{2, 2}, {2, 3}}};
  Matrix<Real, m_dofs_per_cell, m_dofs_per_cell*2> mat({1, 2, 3,
                                                        4, 5, 6,
                                                        7, 8, 9,
                                                        10, 11, 12,
                                                        13, 14, 15,
                                                        16, 17, 18});
  m_assembler->assembleValues(row_indices, col_indices, mat);

  AssemblerType assembler2(*m_assembler);
  assembler2.assembleValues(row_indices, col_indices, mat);

  const auto& dof_nums = m_disc->getDofNumbering()->getData(0);
  UInt row = 0;
  for (linear_system::Indices row_idx : row_indices)
    for (UInt k1=0; k1 < m_dofs_per_cell; ++k1)
    {
      UInt col = 0;
      UInt dof_row = dof_nums(row_idx.i, row_idx.j, k1);
      for (linear_system::Indices col_idx : col_indices)
        for (UInt k2=0; k2 < m_dofs_per_cell; ++k2)
        {
          UInt dof_col = dof_nums(col_idx.i, col_idx.j, k2);
          EXPECT_EQ((*m_matrix)(dof_row, dof_col), 2*mat(row, col));
          col++;
        }

      row++;
    }
}