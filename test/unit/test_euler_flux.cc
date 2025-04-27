#include "gtest/gtest.h"
#include "physics/euler/euler_flux.h"

using namespace structured_fv;
using namespace structured_fv::euler;

namespace {
  using Complex = std::complex<Real>;
}

TEST(EulerFlux, Pressure)
{
  // standard atmosphere at 0 altitude
  Vec4<Real> prim_vars = {1.225, 0, 0, 288.16};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  EXPECT_NEAR(compute_pressure(q), 1.01325E5, 100.0);
}

TEST(EulerFlux, PressureJac)
{
  Vec4<Real> prim_vars = {1.225, 10, 20, 288.16};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  auto [p, p_jac] = compute_pressure_jac(q);
  Real h=1e-20;
  auto q_c = convert<Complex>(q);
  for (int i=0; i < 4; ++i)
  {
    q_c[i] += Complex(0, h);
    Complex p_dot = compute_pressure(q_c);
    EXPECT_DOUBLE_EQ(p_dot.imag()/h, p_jac[i]);
    q_c[i] -= Complex(0, h);
  }
}

TEST(EulerFlux, NormalVelocity)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  EXPECT_NEAR(compute_un(q, {2, 3}), 5*2 + 10*3, 1e-13);
}

TEST(EulerFlux, Xdirection)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  Real p = compute_pressure(q);
  auto flux = compute_euler_flux(q, {1, 0});

  Real E = prim_vars[0]*(Cv*prim_vars[3] + 0.5*(prim_vars[1]*prim_vars[1] + prim_vars[2]*prim_vars[2]));

  EXPECT_NEAR(flux[0], prim_vars[0]*prim_vars[1], 1e-13);
  EXPECT_NEAR(flux[1], prim_vars[0]*prim_vars[1]*prim_vars[1] + p, 1e-13);
  EXPECT_NEAR(flux[2], prim_vars[0]*prim_vars[1]*prim_vars[2], 1e-13);
  EXPECT_NEAR(flux[3], (E + p)*prim_vars[1], 1e-13);
}

TEST(EulerFlux, Ydirection)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  Real p = compute_pressure(q);
  auto flux = compute_euler_flux(q, {0, 1});

  Real E = prim_vars[0]*(Cv*prim_vars[3] + 0.5*(prim_vars[1]*prim_vars[1] + prim_vars[2]*prim_vars[2]));

  EXPECT_NEAR(flux[0], prim_vars[0]*prim_vars[2], 1e-13);
  EXPECT_NEAR(flux[1], prim_vars[0]*prim_vars[1]*prim_vars[2], 1e-13);
  EXPECT_NEAR(flux[2], prim_vars[0]*prim_vars[2]*prim_vars[2] + p, 1e-13);
  EXPECT_NEAR(flux[3], (E + p)*prim_vars[2], 1e-13);
}


TEST(EulerFlux, Jacobian)
{
  Vec4<Real> prim_vars = {1.225, 10, 20, 288.16};
  std::array<Real, 2> normal = {2, 3};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  auto [flux, flux_jac] = compute_euler_flux_jac(q, normal);

  Real h = 1e-20;
  auto q_c = convert<Complex>(q);
  for (int i=0; i < 4; ++i)
  {
    q_c[i] += Complex(0, h);
    Vec4<Complex> flux_dot = compute_euler_flux(q_c, normal);
    for (int j=0; j < 4; ++j)
    {
      EXPECT_DOUBLE_EQ(flux_dot[j].imag()/h, flux_jac(j, i));
    }
    q_c[i] -= Complex(0, h);
  }
}

TEST(EulerFlux, EigenDecompRightEigenVectors)
{
  Vec4<Real> prim_vars = {2, 10, 20, 300};
  std::array<Real, 2> normal = {2, 3};  //TODO: was 2, 3
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  
  Matrix<Real, 4, 4> R, Rinv;
  Vec4<Real> lambda;
  auto [flux, flux_jac] = compute_euler_flux_jac(q, normal);
  computeEigenDecomp(q, normal, R, lambda, Rinv);

  // The columns of R should be the right eigenvectors of the flux Jacobian,
  // which are defined as A*v = lambda*v
  for (int j=0; j < 4; ++j)
  {
    Vec4<Real> Av{0, 0, 0, 0};
    for (int i=0; i < 4; ++i)
      for (int k=0; k < 4; ++k)
        Av[i] += flux_jac(i, k) * R(k, j);

    for (int k=0; k < 4; ++k)
    {
      EXPECT_NEAR(Av[k], lambda[j]*R(k, j), 1e-6);
    }
  }
}


TEST(EulerFlux, EigenDecompLeftEigenVectors)
{
  Vec4<Real> prim_vars = {2, 0, 0, 300};
  std::array<Real, 2> normal = {1, 0};  //TODO: was 2, 3
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  
  Matrix<Real, 4, 4> R, Rinv;
  Vec4<Real> lambda;
  auto [flux, flux_jac] = compute_euler_flux_jac(q, normal);
  computeEigenDecomp(q, normal, R, lambda, Rinv);

  // The rows of Rinv should be the left eigenvectors of the flux Jacobian,
  // which are defined as v^T * A = lambda*v^T
  for (int j=0; j < 4; ++j)
  {
    Vec4<Real> vA{0, 0, 0, 0};
    for (int i=0; i < 4; ++i)
      for (int k=0; k < 4; ++k)
        vA[i] += Rinv(j, k) * flux_jac(k, i);

    for (int k=0; k < 4; ++k)
    {
      EXPECT_NEAR(vA[k], lambda[j]*Rinv(j, k), 1e-6);
    }
  }
}

TEST(EulerFlux, EigenDecompLeftEigenVectorsInv)
{
  Vec4<Real> prim_vars = {2, 20, 30, 300};
  std::array<Real, 2> normal = {2, 3};  //TODO: was 2, 3
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  
  Matrix<Real, 4, 4> R, Rinv;
  Vec4<Real> lambda;
  auto [flux, flux_jac] = compute_euler_flux_jac(q, normal);
  computeEigenDecomp(q, normal, R, lambda, Rinv);

  Matrix<Real, 4, 4> identity;

  for (int i=0; i < 4; ++i)
    for (int j=0; j < 4; ++j)
      for (int k=0; k < 4; ++k)
        identity(i, j) += Rinv(i, k) * R(k, j);

  for (int i=0; i < 4; ++i)
    for (int j=0; j < 4; ++j)
      if (i == j)
      {
        EXPECT_NEAR(identity(i, j), 1.0, 1e-13);
      } else
      {
        EXPECT_NEAR(identity(i, j), 0.0, 1e-13);
      }
}
