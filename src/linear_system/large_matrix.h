#ifndef STRUCTURED_FINITE_VOLUME_LINEAR_SYSTEM_LARGE_MATRIX_H
#define STRUCTURED_FINITE_VOLUME_LINEAR_SYSTEM_LARGE_MATRIX_H

#include "disc/discretization.h"
#include "utils/project_defs.h"
#include "utils/matrix.h"
#include "disc/disc_vector.h"
#include "linear_system/sparsity_pattern.h"
#include "assembler_base.h"
#include "utils/error_handling.h"

namespace structured_fv {
namespace linear_system {

struct LargeMatrixOpts
{
  virtual ~LargeMatrixOpts() = default;

  bool is_structurally_symmetric = false;
  bool is_value_symmetric        = false;
  bool factor_in_place           = false;
};

using VectorType = disc::DiscVector<Real>::ViewType;
using ConstVectorType = const disc::DiscVector<Real>::ConstViewType;


//TODO: should m_local and n_local be specified by SparsityPattern?
// The implementations of this class need to have shallow-copy semantics
class LargeMatrix
{
  public:
    LargeMatrix(GlobalDof m_local, GlobalDof n_local, std::shared_ptr<SparsityPattern> sparsity) :
      m_mlocal(m_local),
      m_nlocal(n_local),
      m_sparsity_pattern(sparsity),
      m_is_factored(false)
    {}

    virtual ~LargeMatrix() = default;

    LargeMatrix(const LargeMatrix&) = default;

    LargeMatrix& operator=(const LargeMatrix&) = default;

    virtual std::shared_ptr<LargeMatrix> clone() = 0;

    GlobalDof getMLocal() const { return m_mlocal; }

    GlobalDof getNLocal() const { return m_nlocal; }

    FixedVec<GlobalDof, 2> getSize() const { return {m_mlocal, m_nlocal}; }

    std::shared_ptr<SparsityPattern> getSparsityPattern() const { return m_sparsity_pattern; }

    void zeroMatrix()
    {
      m_is_factored = false;
      zeroMatrix_impl();
    }
/*
    // if dof < 0, the corresponding entries are ignored
    template <UInt M, UInt N>
    void assembleValues(const FixedVec<GlobalDof, M>& dofs_row, const FixedVec<GlobalDof, N>& dofs_col,  const Matrix<Real, M, N>& jac)
    {
#ifndef NDEBUG
      assert(!m_is_factored);
      checkDofUniqueness(dofs_row);
      checkDofUniqueness(dofs_col);
#endif

      assembleValues_impl(dofs_row, dofs_col, jac);
    }
*/       
    void finishMatrixAssembly() { finishMatrixAssembly_impl(); };

    // factor the matrix/update the preconditioner.
    void factor()
    {
      factor_impl();
      m_is_factored = true;
    }

    // solve Ax = b for x.  The initial value of x might be used for an initial guess.
    // factor() must have been called first
    void solve(ConstVectorType& b, VectorType& x)
    {
      assertAlways(m_is_factored, "must factor matrix before solving a linear system");
      solve_impl(b, x);
    }

    // compute b = A * x, overwriting b
    void matVec(ConstVectorType & x, VectorType& b)
    {
      assert(x.size() >= getNLocal());
      assert(b.size() >= getMLocal());
      matVec_impl(x, b);
    }

    void axpy(Real alpha, std::shared_ptr<LargeMatrix> x)
    {
      m_is_factored = false;
      axpy_impl(alpha, x);
    }

    // Assembler is the base class of an Assembler<LargeMatrixType>, which can access the
    // temlated assembleValues() method
    AssemblerBasePtr getAssembler(disc::StructuredDiscPtr disc) const { return getAssembler_impl(disc); }


  protected:

    bool getIsFactored() const { return m_is_factored; }

    virtual void zeroMatrix_impl() = 0;

    virtual void finishMatrixAssembly_impl() {};

    // factor the matrix/update the preconditioner.
    virtual void factor_impl() = 0;

    // solve Ax = b for x.  The initial value of x might be used for an initial guess.
    // factor() must have been called first
    virtual void solve_impl(ConstVectorType & b, VectorType& x) = 0;

    virtual void matVec_impl(ConstVectorType& x, VectorType& b) = 0;

    virtual void axpy_impl(Real alpha, std::shared_ptr<LargeMatrix> x) = 0;

    virtual AssemblerBasePtr getAssembler_impl(disc::StructuredDiscPtr disc) const = 0;

    template <UInt M, UInt N>
    void checkDofsForAssembly(const FixedVec<GlobalDof, M>& dofs_row, const FixedVec<GlobalDof, N>& dofs_col) const
    {
#ifndef NDEBUG
      assert(!m_is_factored);
      checkDofUniqueness(dofs_row);
      checkDofUniqueness(dofs_col);
#endif
    }

    template <UInt M>
    void checkDofUniqueness(const FixedVec<GlobalDof, M>& dofs) const
    {
      std::vector<GlobalDof> dofs_copy;
      for (unsigned int i=0; i < dofs.size(); ++i)
        if (dofs[i] >= 0)
          dofs_copy.push_back(i);

      auto it = std::unique(dofs_copy.begin(), dofs_copy.end());
      assertAlways(it == dofs_copy.end(), "dofs passed to LargeMatrix are not unique");
    }

  private:
    GlobalDof m_mlocal;
    GlobalDof m_nlocal;
    std::shared_ptr<SparsityPattern> m_sparsity_pattern;
    bool m_is_factored;
};

using LargeMatrixPtr = std::shared_ptr<LargeMatrix>;

void solve(LargeMatrixPtr mat, const disc::DiscVectorPtr<Real> b, disc::DiscVectorPtr<Real> x);

} // namespace
}

#endif
