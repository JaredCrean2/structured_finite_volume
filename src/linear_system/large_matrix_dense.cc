#include "linear_system/large_matrix_dense.h"
#include "linear_system/bla_wrapper.h"

namespace structured_fv {
namespace linear_system {

void LargeMatrixDense::printToStdout()
{
  auto& self = *this;
  for (int i=0; i < getMLocal(); ++i)
  {
    for (int j=0; j < getNLocal(); ++j)
      std::cout << self(i, j) << ", ";

    std::cout << std::endl;
  }
}

AssemblerPtr<LargeMatrixDense> LargeMatrixDense::getAssembler(disc::StructuredDiscPtr disc) const
{
  return std::make_shared<Assembler<LargeMatrixDense>>(disc, *this);
}


void LargeMatrixDense::zeroMatrix_impl()
{
  Kokkos::deep_copy(m_matrix, 0);
}

void LargeMatrixDense::factor_impl()
{
  
  Matrix* matrix_factorization;
  if (m_opts.factor_in_place)
  {
    matrix_factorization = &m_matrix;
  } else
  {
    matrix_factorization = &m_matrix_factorization;
    Kokkos::deep_copy(*matrix_factorization, m_matrix);
    //*matrix_factorization = m_matrix;
  }


  if (m_opts.is_value_symmetric)
  {
    potrf('L', *matrix_factorization);
  } else
  {
    //m_ipiv.resize(getMLocal());
    if (m_ipiv.extent(0) == 0)
      Kokkos::resize(m_ipiv, getMLocal());
    getrf(*matrix_factorization, m_ipiv);
  }

}



void LargeMatrixDense::solve_impl(ConstVectorType& b, VectorType& x)
{
  assert(b.extent(0) == getMLocal());
  assert(x.extent(0) == getMLocal());

  //std::copy(&(b[0]), (&b[0] + getMLocal()), &(x[0]));
  Kokkos::deep_copy(x, b);

  auto& matrix_factorization = m_opts.factor_in_place ? m_matrix : m_matrix_factorization;

  Kokkos::View<const double**, Kokkos::LayoutLeft> mat_copy(matrix_factorization);
  if (m_opts.is_value_symmetric)
  {
    //potrs('L', getMLocal(), 1, matrix_factorization.data(), getMLocal(), &(x[0]), getMLocal());
    potrs('L', matrix_factorization, x);
  } else
  {
    getrs('N', matrix_factorization, m_ipiv, x);
  }
}


void LargeMatrixDense::matVec_impl(ConstVectorType& x, VectorType& b)
{
  assertAlways(!getIsFactored(), "cannot do mat-vec when matrix is factored");

  if (m_opts.is_value_symmetric)
    symv('L', 1.0, m_matrix, x, 0.0, b);
  else
    gemv('N', 1, m_matrix, x, 0, b);
}


void LargeMatrixDense::axpy_impl(Real alpha, std::shared_ptr<LargeMatrix> x)
{
  auto x_dense = std::dynamic_pointer_cast<LargeMatrixDense>(x);
  for (int j=0; j < getNLocal(); ++j)
    for (int i=0; i < getMLocal(); ++i)
      m_matrix(i, j) += x_dense->m_matrix(i, j)*alpha;    
}

AssemblerBasePtr LargeMatrixDense::getAssembler_impl(disc::StructuredDiscPtr disc) const
{
  return getAssembler(disc);
}


std::ostream& operator<<(std::ostream& os, const LargeMatrixDense& mat)
{
  for (GlobalDof i=0; i < mat.getMLocal(); ++i)
  {
    for (GlobalDof j=0; j < mat.getNLocal(); ++j)
    {
      os << mat.m_matrix(i, j);
      if (j != mat.getNLocal() - 1)
        os << ", ";
    }
    if (i != mat.getMLocal() - 1)
      os << std::endl;
  }

  return os;
}


}  // namespace
}