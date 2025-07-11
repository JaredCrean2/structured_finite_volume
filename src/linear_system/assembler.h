#ifndef STRUCTURED_FINITE_VOLUME_LINEAR_SYSTEM_ASSEMBLER_H
#define STRUCTURED_FINITE_VOLUME_LINEAR_SYSTEM_ASSEMBLER_H

#include "assembler_base.h"
#include "utils/project_defs.h"
#include "utils/matrix.h"
#include "disc/discretization.h"
#include "disc/elem_field.h"
#include "disc/discretization.h"
#include <View/Kokkos_ViewCtor.hpp>
#include <cassert>
#include <iosfwd>

namespace structured_fv {
namespace linear_system {

static_assert(std::is_signed<GlobalDof>::value, "DofInt must be a signed integral type");

struct Indices
{
  UInt i;
  UInt j;
};

std::ostream& operator<<(std::ostream& os, const Indices& ind);

template <typename LargeMatrixType>
class Assembler final : public AssemblerBase
{
  public:
    explicit Assembler(disc::StructuredDiscPtr disc, LargeMatrixType matrix) :
      m_disc(*disc),
      m_matrix(matrix),
      m_dof_nums()
    {}

    ~Assembler() = default;

    Assembler operator=(const Assembler& other) = delete;

    void setAlpha(Real alpha) override { m_alpha = alpha; }

    Real getAlpha() const override { return m_alpha; }

    void setBlock(UInt block_id)
    {
      m_dof_nums = m_disc.getDofNumbering()->getData(block_id);
    }

    template <UInt M, UInt N, UInt K, UInt P>
    void assembleValues(const FixedVec<Indices, M>& row_indices, const FixedVec<Indices, N>& col_indices, Matrix<Real, K, P>& jac)
    {
      static_assert(K % M == 0, "size of jacobian contribution should be num_dof_per_node times number of cells");
      static_assert(P % M == 0, "size of jacobian contribution should be num_dof_per_node times number of cells");
      static_assert(K/M == P/N, "number of dofs per node must match for rows and columns");
      assert(K/M == m_dof_nums.extent(2));
      assert(P/N == m_dof_nums.extent(2));
      constexpr UInt num_dof_per_node = K/M;

      for (UInt i=0; i < K; ++i)
        for (UInt j=0; j < P; ++j)
          jac(i, j) *= m_alpha;

      FixedVec<GlobalDof, K> row_dofs;
      FixedVec<GlobalDof, P> col_dofs;
      for (UInt i=0; i < M; ++i)
        for (UInt k=0; k < num_dof_per_node; ++k)
          row_dofs[i*num_dof_per_node + k] = m_dof_nums(row_indices[i].i, row_indices[i].j, k);

      for (UInt j=0; j < N; ++j)
        for (UInt k=0; k < num_dof_per_node; ++k)
          col_dofs[j*num_dof_per_node + k] = m_dof_nums(col_indices[j].i, col_indices[j].j, k);

      m_matrix.assembleValues(row_dofs, col_dofs, jac);  
    }

  private:
    disc::StructuredDisc& m_disc;
    LargeMatrixType m_matrix;
    disc::ElementField<GlobalDof>::ConstFieldData m_dof_nums;
    Real m_alpha = 1;
};

template <typename LargeMatrixType>
using AssemblerPtr = std::shared_ptr<Assembler<LargeMatrixType>>;

}
}

#endif