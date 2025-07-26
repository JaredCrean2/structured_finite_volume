#ifndef STRUCTURED_FINITE_VOLUME_TEST_UTILS_JACOBIAN_H
#define STRUCTURED_FINITE_VOLUME_TEST_UTILS_JACOBIAN_H

#include "physics/euler/typedefs.h"
#include "utils/array3.h"
#include "utils/matrix.h"
#include "utils/project_defs.h"
#include "utils/vec.h"
#include <limits>

namespace test_utils {

using namespace structured_fv;

using ScalarVectorPair = structured_fv::euler::ScalarVectorPair<Real>;
using VectorMatrixPair = structured_fv::euler::VectorMatrixPair<Real>;
using MatrixArray3Pair = std::pair<Matrix<Real, 4>, Array3<Real, 4>>;

inline Real get_tol(Real v1, Real v2, int num_ulp)
{
  Real min_val = std::min(std::abs(v1), std::abs(v2));
  Real eps = std::nextafter(min_val, std::numeric_limits<Real>::max()) - min_val;
  Real tol = std::max(num_ulp * eps, 1e-13);

  return tol;
}

#define EXPECT_DOUBLE_EQ_CUSTOM(v1, v2, num_ulp)  \
{                                                 \
  EXPECT_NEAR(v1, v2, test_utils::get_tol(v1, v2, num_ulp));  \
}                                                 \

inline Vec4<Complex> make_complex(const Vec4<Real>& q)
{
  Vec4<Complex> qc;
  for (UInt i=0; i < q.size(); ++i)
    qc[i] = q[i];

  return qc;
}

inline void zero_complex(Vec4<Complex>& qc)
{
  for (UInt i=0; i < qc.size(); ++i)
    qc[i].imag(0);
}

template <typename Func>
constexpr bool IsScalarFunc = std::is_invocable_r_v<Complex, Func, Vec4<Complex>>;

template <typename Func>
constexpr bool IsVectorFunc = std::is_invocable_r_v<Vec4<Complex>, Func, Vec4<Complex>>;

template <typename Func>
constexpr bool IsMatrixFunc = std::is_invocable_r_v<Matrix<Complex, 4>, Func, Vec4<Complex>>;

template <typename Func>
constexpr bool IsScalarJac = std::is_invocable_r_v<ScalarVectorPair, Func, Vec4<Real>>;

template <typename Func>
constexpr bool IsVectorJac = std::is_invocable_r_v<VectorMatrixPair, Func, Vec4<Real>>;

template <typename Func>
constexpr bool IsMatrixJac = std::is_invocable_r_v<MatrixArray3Pair, Func, Vec4<Real>>;

// Scalar func must be y = f(q), where y is a scalar.  Must be callable with complex
// numbers (only).
// func_jac must be a function y, dy/dq = func_jac(q, q_dot) that returns both the value and the jacobian
// of the function wrt q, with each entry scaled by the corresponding entry fo q_dot.  Must be callable with real numbers only
template <typename ScalarFunc, typename ScalarFuncJac>
void checkJacobianScalar(const Vec4<Real>& q, ScalarFunc func, ScalarFuncJac func_jac)
{
  static_assert(IsScalarFunc<ScalarFunc>);
  static_assert(IsScalarJac<ScalarFuncJac>);
  Vec4<Complex> qc = make_complex(q);
  Real h = 1e-80;
  Complex pert(0, h);

  Complex y = func(qc);

  Vec4<Real> dydq_cs;
  for (UInt i=0; i < q.size(); ++i)
  {
    qc[i] += pert;
    Complex yc = func(qc);
    dydq_cs[i] = yc.imag()/h;
    qc[i] -= pert;
  }

  const auto [y2, dydq_jac] = func_jac(q);
  EXPECT_NEAR(y.real(), y2, 1e-13);
  for (UInt i=0; i < dydq_cs.size(); ++i)
    EXPECT_DOUBLE_EQ(dydq_jac[i], dydq_cs[i]);
}

template <typename VectorFunc, typename VectorFuncJac>
void checkJacobianVector(const Vec4<Real>& q, VectorFunc func, VectorFuncJac func_jac)
{
  static_assert(IsVectorFunc<VectorFunc>);
  static_assert(IsVectorJac<VectorFuncJac>);

  Vec4<Complex> qc = make_complex(q);
  Real h = 1e-80;
  Complex pert(0, h);

  Vec4<Complex> y = func(qc);

  Matrix<Real, 4> dydq_cs;
  for (UInt k=0; k < q.size(); ++k)
  {
    qc[k] += pert;
    Vec4<Complex> yc = func(qc);
    for (UInt i=0; i < q.size(); ++i)
    {
      dydq_cs(i, k) = yc[i].imag()/h;
    }

    qc[k] -= pert;
  }

  const auto [y2, dydq_jac] = func_jac(q);
  for (UInt i=0; i < q.size(); ++i)
  {
    EXPECT_DOUBLE_EQ(y[i].real(), y2[i]);
    for (UInt j=0; j < q.size(); ++j)
    {
      //Real min_val = std::min(std::abs(dydq_jac(i, j)), std::abs(dydq_cs(i, j)));
      //Real eps = std::nextafter(min_val, std::numeric_limits<Real>::max()) - min_val;
      //Real tol = std::max(300 * eps, 1e-13);
      //EXPECT_NEAR(dydq_jac(i, j), dydq_cs(i, j), tol);
      EXPECT_DOUBLE_EQ_CUSTOM(dydq_jac(i, j), dydq_cs(i, j), 300);

    }
  }
}

template <typename MatrixFunc, typename MatrixFuncJac>
void checkJacobianMatrix(const Vec4<Real>& q, MatrixFunc func, MatrixFuncJac func_jac)
{
  static_assert(IsMatrixFunc<MatrixFunc>);
  static_assert(IsMatrixJac<MatrixFuncJac>);

  Vec4<Complex> qc = make_complex(q);
  Real h = 1e-80;
  Complex pert(0, h);

  Matrix<Complex, 4> y = func(qc);

  Array3<Real, 4> dydq_cs;
  for (UInt k=0; k < q.size(); ++k)
  {
    qc[k] += pert;
    Matrix<Complex, 4> yc = func(qc);
    for (UInt i=0; i < q.size(); ++i)
      for (UInt j=0; j < q.size(); ++j)
        dydq_cs(i, j, k) = yc(i, j).imag()/h;

    qc[k] -= pert;
  }

  const auto [y2, dydq_jac] = func_jac(q);
  for (UInt i=0; i < q.size(); ++i)
    for (UInt j=0; j < q.size(); ++j)
    {
      EXPECT_DOUBLE_EQ(y(i, j).real(), y2(i, j));
      for (UInt k=0; k < q.size(); ++k)
      {
        //Real min_val = std::min(std::abs(dydq_jac(i, j, k)), std::abs(dydq_cs(i, j, k)));
        //Real eps = std::nextafter(min_val, std::numeric_limits<Real>::max()) - min_val;
        //Real tol = std::max(8 * eps, 1e-13);
        //EXPECT_NEAR(dydq_jac(i, j, k), dydq_cs(i, j, k), tol);
        EXPECT_DOUBLE_EQ_CUSTOM(dydq_jac(i, j, k), dydq_cs(i, j, k), 300);

      }
    }
}

}

#endif