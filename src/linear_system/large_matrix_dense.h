#ifndef STRUCTURED_FINITE_VOLUME_LINEAR_SYSTEM_LARGE_MATRIX_DENSE_H
#define STRUCTURED_FINITE_VOLUME_LINEAR_SYSTEM_LARGE_MATRIX_DENSE_H

#include "large_matrix.h"
#include "linear_system/sparsity_pattern_dense.h"
#include "assembler.h"
#include "bla_wrapper.h"
//#include <decl/Kokkos_Declare_SERIAL.hpp>
#include <memory>
#include <iosfwd>

namespace structured_fv {
namespace linear_system {

class LargeMatrixDense final : public LargeMatrix
{
  public:
    LargeMatrixDense(const std::string& name, GlobalDof mlocal, GlobalDof nlocal, LargeMatrixOpts opts) :
      LargeMatrix(mlocal, nlocal, std::make_shared<SparsityPatternDense>(mlocal)),
      m_name(name),
      m_opts(opts),
      m_matrix(name, mlocal, nlocal),
      m_matrix_factorization(name + "_factorization", 0, 0),
      m_ipiv(name + "_ipiv", 0)
    {
      if (!opts.factor_in_place)
        Kokkos::resize(m_matrix_factorization, mlocal, nlocal);
    }

    const Real& operator()(int i, int j) const { return m_matrix(i, j); }

    std::shared_ptr<LargeMatrix> clone() override
    {
      return std::make_shared<LargeMatrixDense>(m_name, getMLocal(), getNLocal(), m_opts);
    }

    void printToStdout();

    template <UInt M, UInt N, size_t M2, size_t N2>
    void assembleValues(const std::array<GlobalDof, M2>& dofs_row, const std::array<GlobalDof, N2>& dofs_col,
                        const Matrix<Real, M, N>& jac)
    {
      static_assert(M == M2);
      static_assert(N == N2);
      checkDofsForAssembly(dofs_row, dofs_col);

      for (unsigned int j=0; j < dofs_col.size(); ++j)
      {
        if (dofs_col[j] < 0)
          continue;

        for (unsigned int i=0; i < dofs_row.size(); ++i)
          if (dofs_row[i] >= 0)
            m_matrix(dofs_row[i], dofs_col[j]) += jac(i, j);
            //m_matrix[getIdx(dofs_rows[i], dofs_cols[j])] += jac[i][j];
      }
    }

    AssemblerPtr<LargeMatrixDense> getAssembler(disc::StructuredDiscPtr disc) const;


  protected:
    void zeroMatrix_impl() override;

    void factor_impl() override;

    void solve_impl(ConstVectorType& b, VectorType& x) override;

    void matVec_impl(ConstVectorType& x, VectorType& b) override;

    void axpy_impl(Real alpha, std::shared_ptr<LargeMatrix> x) override;

    AssemblerBasePtr getAssembler_impl(disc::StructuredDiscPtr disc) const override;

  private:
    using Matrix = Kokkos::View<Real**, Kokkos::LayoutLeft>;

    std::string m_name;
    LargeMatrixOpts m_opts;
    Matrix m_matrix;
    Matrix m_matrix_factorization;
    Kokkos::View<lapack_int*, Kokkos::LayoutLeft> m_ipiv;
    //std::vector<Real> m_matrix;
    //std::vector<Real> m_matrix_factorization;
    //std::vector<lapack_int> m_ipiv;

    friend std::ostream& operator<<(std::ostream& os, const LargeMatrixDense& mat);

};

std::ostream& operator<<(std::ostream& os, const LargeMatrixDense& mat);


} // namespace
}

#endif