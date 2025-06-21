#include "gtest/gtest.h"
#include "physics/euler/euler_flux.h"
#include "jacobian.h"

using namespace structured_fv;
using namespace structured_fv::euler;

namespace {
  using Complex = std::complex<Real>;
}

TEST(EulerFlux, Pressure)
{
  // standard atmosphere at 0 altitude
  Vec4<Real> prim_vars = {1.225, 5, 10, 288.16};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  EXPECT_NEAR(compute_pressure(q), 1.01325E5, 100.0);
}

TEST(EulerFlux, PressureJac)
{
  //TODO: make velocity nonzero
  Vec4<Real> prim_vars = {1.225, 5, 10, 288.16};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  test_utils::checkJacobianScalar(q, &compute_pressure<Complex>, &compute_pressure_jac<Real>);
}

TEST(EulerFlux, Sos2Jac)
{
  Vec4<Real> prim_vars = {1.225, 5, 10, 288.16};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  test_utils::checkJacobianScalar(q, &compute_sos2<Complex>, &compute_sos2_jac<Real>);
}

TEST(EulerFlux, SosJac)
{
  Vec4<Real> prim_vars = {1.225, 5, 10, 288.16};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  test_utils::checkJacobianScalar(q, &compute_sos<Complex>, &compute_sos_jac<Real>);
}

TEST(EulerFlux, NormalMomenutmJac)
{
  Vec4<Real> prim_vars = {1.225, 5, 10, 288.16};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  auto func = [&](auto q) { return compute_momentum_n(q, {2, 3}); };
  auto func_jac = [&](auto q) { return compute_momentum_n_jac(q, {2, 3}); };
  test_utils::checkJacobianScalar(q, func, func_jac);
}


TEST(EulerFlux, ComputeTemperatureJac)
{
  Vec4<Real> prim_vars = {1.225, 5, 10, 288.16};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  test_utils::checkJacobianScalar(q, &compute_temperature<Complex>, &compute_temperature_jac<Real>);
}

TEST(EulerFlux, PrimitiveVariables)
{
  Vec4<Real> prim_vars = {1.225, 100, 200, 288.16};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());  
  Vec4<Real> prim_vars2 = compute_primitive_variables(q);
  for (UInt i=0; i < 4; ++i)
    EXPECT_NEAR(prim_vars[i], prim_vars2[i], 1e-13);
}


TEST(EulerFlux, ComputeConservativeVarsJac)
{
  Vec4<Real> prim_vars = {1.225, 5, 10, 288.16};
  auto func = [&](auto primitive) { return compute_conservative_variables(primitive, PrimitiveVarTag{}); };
  auto func_jac = [&](auto primitive) { return compute_conservative_variables_jac(primitive, PrimitiveVarTag{}); };
  test_utils::checkJacobianVector(prim_vars, func, func_jac);
}

TEST(EulerFlux, ComputePrimitiveVarsJac)
{
  Vec4<Real> prim_vars = {1.225, 5, 10, 288.16};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  auto func = [&](auto q) { return compute_primitive_variables(q, ConservativeVarTag{}); };
  auto func_jac = [&](auto q) { return compute_primitive_variables_jac(q, ConservativeVarTag{}); };
  test_utils::checkJacobianVector(q, func, func_jac);
}

/*
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
*/

TEST(EulerFlux, NormalVelocity)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  EXPECT_NEAR(compute_un(q, {2, 3}), 5*2 + 10*3, 1e-13);
}

TEST(EulerFlux, NormalVelocityJac)
{
  Vec4<Real> prim_vars = {1.225, 5, 10, 288.16};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  auto func = [&](auto q) { return compute_un(q, {2, 3}); };
  auto func_jac = [&](auto q) { return compute_un_jac(q, {2, 3}); };
  test_utils::checkJacobianScalar(q, func, func_jac);
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

TEST(EulerFlux, EulerFluxJac)
{
  Vec4<Real> prim_vars = {1.225, 5, 10, 288.16};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  auto func = [&](auto q) { return compute_euler_flux(q, {2, 3}); };
  auto func_jac = [&](auto q) { return compute_euler_flux_jac(q, {2, 3}); };
  test_utils::checkJacobianVector(q, func, func_jac);
}

TEST(EulerFlux, EigenDecompRightEigenVectors)
{
  Vec4<Real> prim_vars = {2, 10, 20, 300};
  FixedVec<Real, 2> normal = {2, 3};  //TODO: was 2, 3
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  
  Matrix<Real, 4, 4> R, Rinv;
  Vec4<Real> lambda;
  auto [flux, flux_jac] = compute_euler_flux_jac(q, normal);
  compute_eigen_decomp(q, normal, R, lambda, Rinv);

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

TEST(EulerFlux, EigenDecompRightEigenVectorsJac)
{
  Vec4<Real> prim_vars = {2, 10, 20, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());

  auto func = [](auto q)
  {
    Matrix<Complex, 4, 4> R, Rinv;
    Vec4<Complex> lambda;
    euler::compute_eigen_decomp(q, {3, 4}, R, lambda, Rinv);
    return R;
  };

  auto func_jac = [](auto q)
  {
    Matrix<Real, 4, 4> R, Rinv;
    Vec4<Real> lambda;

    Matrix<Real, 4, 4> lambda_jac;
    Array3<Real, 4> R_jac, Rinv_jac;
    euler::compute_eigen_decomp_jac(q, {3, 4}, R, R_jac, lambda, lambda_jac, Rinv, Rinv_jac);
    return std::make_pair(R, R_jac);
  };
  test_utils::checkJacobianMatrix(q, func, func_jac);
}


TEST(EulerFlux, EigenDecompLeftEigenVectors)
{
  Vec4<Real> prim_vars = {2, 20, 30, 300};
  FixedVec<Real, 2> normal = {1, 0};  //TODO: was 2, 3
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  
  Matrix<Real, 4, 4> R, Rinv;
  Vec4<Real> lambda;
  auto [flux, flux_jac] = compute_euler_flux_jac(q, normal);
  compute_eigen_decomp(q, normal, R, lambda, Rinv);

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

TEST(EulerFlux, EigenDecompLeftEigenVectorsJac)
{
  Vec4<Real> prim_vars = {2, 10, 20, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());

  auto func = [](auto q)
  {
    Matrix<Complex, 4, 4> R, Rinv;
    Vec4<Complex> lambda;
    euler::compute_eigen_decomp(q, {3, 4}, R, lambda, Rinv);
    return Rinv;
  };

  auto func_jac = [](auto q)
  {
    Matrix<Real, 4, 4> R, Rinv;
    Vec4<Real> lambda;

    Matrix<Real, 4, 4> lambda_jac;
    Array3<Real, 4> R_jac, Rinv_jac;
    euler::compute_eigen_decomp_jac(q, {3, 4}, R, R_jac, lambda, lambda_jac, Rinv, Rinv_jac);
    return std::make_pair(Rinv, Rinv_jac);
  };
  test_utils::checkJacobianMatrix(q, func, func_jac);
}

TEST(EulerFlux, EigenDecompLeftEigenVectorsInv)
{
  Vec4<Real> prim_vars = {2, 20, 30, 300};
  FixedVec<Real, 2> normal = {2, 3};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  
  Matrix<Real, 4, 4> R, Rinv;
  Vec4<Real> lambda;
  compute_eigen_decomp(q, normal, R, lambda, Rinv);

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

TEST(EulerFlux, EigenDecompEigenvaluesSorted)
{
  Vec4<Real> prim_vars = {2, 0, 0, 300};
  FixedVec<Real, 2> normal = {1, 0};  //TODO: was 2, 3
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  
  Matrix<Real, 4, 4> R, Rinv;
  Vec4<Real> lambda;
  compute_eigen_decomp(q, normal, R, lambda, Rinv);

  for (int i=1; i < 4; ++i)
    EXPECT_GE(lambda[i], lambda[i-1]); 
}

TEST(EulerFlux, EigenDecompLeftEigenvaluesJac)
{
  Vec4<Real> prim_vars = {2, 10, 20, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());

  auto func = [](auto q)
  {
    Matrix<Complex, 4, 4> R, Rinv;
    Vec4<Complex> lambda;
    euler::compute_eigen_decomp(q, {3, 4}, R, lambda, Rinv);
    return lambda;
  };

  auto func_jac = [](auto q)
  {
    Matrix<Real, 4, 4> R, Rinv;
    Vec4<Real> lambda;

    Matrix<Real, 4, 4> lambda_jac;
    Array3<Real, 4> R_jac, Rinv_jac;
    euler::compute_eigen_decomp_jac(q, {3, 4}, R, R_jac, lambda, lambda_jac, Rinv, Rinv_jac);
    return std::make_pair(lambda, lambda_jac);
  };
  test_utils::checkJacobianVector(q, func, func_jac);
}

namespace {

void solveLinearizedRiemannProblem(const Vec4<Real>& qL, const Vec4<Real>& qR,
                                   const Vec2<Real>& normal)
{
  Matrix<Real, 4, 4> R, Rinv;
  Vec4<Real> lambda;
  compute_eigen_decomp(qL, normal, R, lambda, Rinv);

  auto delta_q = qR - qL;
  auto alpha = Rinv * delta_q;
  auto fL = compute_euler_flux(qL, normal);

  std::cout << "lambdas = " << lambda << std::endl;
  std::cout << "alphas = " << alpha << std::endl;
  auto q_i = qL;
  auto f_i = fL;
  for (UInt i=0; i < 5; ++i)
  {
    std::cout << "\nstate " << i << ":" << std::endl;

    if (i > 0)
    {
      q_i += alpha[i-1] * R(Column{i-1});
      f_i += lambda[i-1] * alpha[i-1] * R(Column{i-1});
    }

    auto f_qi = compute_euler_flux(q_i, normal);



    std::cout << "q = " << q_i << std::endl;
    std::cout << "flux = " << f_i << std::endl;
    std::cout << "f(q) = " << f_qi << std::endl;
    std::cout << "diff = " << f_i - f_qi << std::endl;

    if (i == 4)
    {
      std::cout << "qR = " << qR << std::endl;
      std::cout << "diff = " << qR - q_i << std::endl;

      auto fR = compute_euler_flux(qR, normal);
      std::cout << "fR = " << fR << std::endl;
      std::cout << "diff = " << f_i - fR << std::endl;
    }
  }

}
}

TEST(EulerFlux, LinearizedRiemannProblem)
{
  FixedVec<Real, 2> normal = {1, 0};

  Vec4<Real> prim_varsL = {2, 400, 0, 300};
  Vec4<Real> prim_varsR = {2, 500, 0, 300};

  auto qL = compute_conservative_variables(prim_varsL, PrimitiveVarTag());
  auto qR = compute_conservative_variables(prim_varsR, PrimitiveVarTag());

  solveLinearizedRiemannProblem(qL, qR, normal);
}


namespace acoustics {

// acoustics coupled with advection
// the state vector is [p, u, phi], where p is the pressure,
// u is velocity, and phi is a passive tracer being advected
// by the fluid

struct AcousticsParams
{
  Real rho;
  Real K;
  Real u0;
};

FixedVec<Real, 3> compute_acoustics_flux(const AcousticsParams& params, const FixedVec<Real, 3>& u)
{
  return {params.u0 * u[0] + params.K*u[1], u[0]/params.rho + params.u0*u[1], params.u0*u[2]};
}

void compute_eigen_decomp(const AcousticsParams& params,
                        Matrix<Real, 3, 3>& R, FixedVec<Real, 3>& lambda)
{
  
  Real c0 = std::sqrt(params.K/params.rho);
  Real z0 = params.rho*c0;

  R(0, 0) = -z0;
  R(0, 1) = 0;
  R(0, 2) = z0;

  R(1, 0) = 1;
  R(1, 1) = 0;
  R(1, 2) = 1;

  R(2, 0) = 0;
  R(2, 1) = 1;
  R(2, 2) = 0;

  lambda[0] = params.u0 - c0;
  lambda[1] = params.u0;
  lambda[2] = params.u0 + c0;
}

FixedVec<Real, 3> solveForAlpha(const AcousticsParams& params, const FixedVec<Real, 3>& uL, const FixedVec<Real, 3>& uR)
{
  Real c0 = std::sqrt(params.K/params.rho);
  Real z0 = params.rho*c0;  

  Real alpha1 = (-(uR[0] - uL[0]) + z0*(uR[1] - uL[1]))/(2*z0);
  Real alpha2 = uR[2] - uL[2];
  Real alpha3 = ( (uR[0] - uL[0]) + z0*(uR[1] - uL[1]))/(2*z0);

  return {alpha1, alpha2, alpha3};
}

void solveLinearizedRiemannProblem(const AcousticsParams& params, const FixedVec<Real, 3>& uL, const FixedVec<Real, 3>& uR)
{
  Matrix<Real, 3, 3> R;
  FixedVec<Real, 3> lambda;
  compute_eigen_decomp(params, R, lambda);
  FixedVec<Real, 3> alpha = solveForAlpha(params, uL, uR);

  std::cout << "lambda = " << lambda << std::endl;
  std::cout << "alpha = " << alpha << std::endl;

  auto fL = compute_acoustics_flux(params, uL);

  auto u_i = uL;
  auto f_i = fL;

  for (UInt i=0; i < 4; ++i)
  {
    std::cout << "\nstate " << i << ":" << std::endl;
    if (i > 0)
    {
      u_i += alpha[i-1] * R(Column{i-1});
      f_i += lambda[i-1] * alpha[i-1] * R(Column{i-1});
    }

    auto f_qi = compute_acoustics_flux(params, u_i);

    std::cout << "q = " << u_i << std::endl;
    std::cout << "flux = " << f_i << std::endl;
    std::cout << "f(q) = " << f_qi << std::endl;
    std::cout << "diff = " << f_i - f_qi << std::endl;

    if (i == 3)
    {
      std::cout << "qR = " << uR << std::endl;
      std::cout << "diff = " << uR - u_i << std::endl;

      auto fR = compute_acoustics_flux(params, uR);
      std::cout << "fR = " << fR << std::endl;
      std::cout << "diff = " << f_i - fR << std::endl;
    }
  }
}

TEST(AcousticsFlux, LinearizedRiemannProblem)
{
  AcousticsParams params{2, 100, 5};
  FixedVec<Real, 3> uL = {10, 20, 2};
  FixedVec<Real, 3> uR = {20, 5, 1};
  solveLinearizedRiemannProblem(params, uL, uR);
}

}
