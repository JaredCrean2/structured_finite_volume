#ifndef STRUCTURED_FINITE_VOLUME_LINEAR_SYSTEM_ASSEMBLER_H
#define STRUCTURED_FINITE_VOLUME_LINEAR_SYSTEM_ASSEMBLER_H

#include "assembler_base.h"
#include "utils/project_defs.h"
#include "utils/matrix.h"
#include "disc/discretization.h"
#include "disc/elem_field.h"
#include "disc/discretization.h"
#include <View/Kokkos_ViewCtor.hpp>

namespace structured_fv {
namespace linear_system {

static_assert(std::is_signed<GlobalDof>::value, "DofInt must be a signed integral type");

struct Indices
{
  UInt i;
  UInt j;
};

template <typename LargeMatrixType>
class Assembler : public AssemblerBase
{
  public:
    explicit Assembler(disc::StructuredDiscPtr disc, LargeMatrixType matrix) :
      m_disc(disc),
      m_matrix(matrix),
      m_dof_nums()
    {}

    virtual ~Assembler() {}

    void setAlpha(Real alpha) override { m_alpha = alpha; }

    Real getAlpha() const override { return m_alpha; }

    void setBlock(UInt block_id)
    {
      m_dof_nums = m_disc->getDofNumbering()->getData(block_id);
    }

    template <UInt M, UInt N>
    void assembleValues(const std::array<Indices, M>& row_indices, const std::array<Indices, N>& col_indices, Matrix<Real, M, N>& jac)
    {
      for (UInt i=0; i < M; ++i)
        for (UInt j=0; j < N; ++j)
          jac[i][j] *= m_alpha;

      std::array<GlobalDof, M> row_dofs;
      std::array<GlobalDof, N> col_dofs;
      for (UInt i=0; i < M; ++i)
        row_dofs[i] = m_dof_nums(row_indices[i].i, row_indices[i].j);

      for (UInt j=0; j < N; ++j)
        col_dofs[j] = m_dof_nums(col_dofs[j].i, col_dofs[j].j);

      m_matrix.assembleValues(row_dofs, col_dofs, jac);  
    }

  private:
    disc::StructuredDiscPtr m_disc;
    LargeMatrixType m_matrix;
    disc::ElementField<GlobalDof>::ConstFieldData m_dof_nums;
    Real m_alpha = 1;
};

template <typename LargeMatrixType>
using AssemblerPtr = std::shared_ptr<Assembler<LargeMatrixType>>;

}
}

#endif