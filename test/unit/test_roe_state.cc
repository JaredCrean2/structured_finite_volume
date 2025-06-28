#include "gtest/gtest.h"
#include "physics/euler/euler_flux.h"
#include "physics/euler/roe_state.h"
#include "physics/euler/typedefs.h"
#include "jacobian.h"

using namespace structured_fv;
using namespace structured_fv::euler;

TEST(RoeAvgState, SameState)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  RoeAvgState<Real> state = compute_roe_avg(q, q);
  Real p = compute_pressure(q);

  EXPECT_NEAR(state.rho, prim_vars[0], 1e-13);
  EXPECT_NEAR(state.u, prim_vars[1], 1e-13);
  EXPECT_NEAR(state.v, prim_vars[2], 1e-13);
  EXPECT_NEAR(state.H, (q[3] + p)/q[0], 1e-10);
  EXPECT_NEAR(compute_pressure(state), p, 1e-10);
  EXPECT_NEAR(compute_sos(state), compute_sos(q), 1e-10);

  auto q2 = compute_conservative_variables(state, RoeStateTag());
  for (UInt i=0; i < DofsPerCell; ++i)
    EXPECT_DOUBLE_EQ(q[i], q2[i]);
}

TEST(RoeAvgState, SpeedOfSound2)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  RoeAvgState<Real> state = compute_roe_avg(q, q);
  
  EXPECT_DOUBLE_EQ(compute_sos2(state), compute_sos2(q));
}

TEST(RoeAvgState, SpeedOfSound2Jac)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  RoeAvgState<Real> state = compute_roe_avg(q, q);

  auto func = [](const Vec4<Complex>& roe_avg)
  {
    return compute_sos2(RoeAvgState<Complex>{roe_avg[0], roe_avg[1], roe_avg[2], roe_avg[3]});
  };

  auto jac = [](const Vec4<Real>& roe_avg)
  {
    return compute_sos2_jac(RoeAvgState<Real>{roe_avg[0], roe_avg[1], roe_avg[2], roe_avg[3]});
  }; 

  test_utils::checkJacobianScalar({state.rho, state.u, state.v, state.H}, func, jac);
}

TEST(RoeAvgState, SpeedOfSound)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  RoeAvgState<Real> state = compute_roe_avg(q, q);
  
  EXPECT_DOUBLE_EQ(compute_sos(state), compute_sos(q));
}

TEST(RoeAvgState, SpeedOfSoundJac)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  RoeAvgState<Real> state = compute_roe_avg(q, q);

  auto func = [](const Vec4<Complex>& roe_avg)
  {
    return compute_sos(RoeAvgState<Complex>{roe_avg[0], roe_avg[1], roe_avg[2], roe_avg[3]});
  };

  auto jac = [](const Vec4<Real>& roe_avg)
  {
    return compute_sos_jac(RoeAvgState<Real>{roe_avg[0], roe_avg[1], roe_avg[2], roe_avg[3]});
  }; 

  test_utils::checkJacobianScalar({state.rho, state.u, state.v, state.H}, func, jac);
}

TEST(RoeAvgState, Pressure)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  RoeAvgState<Real> state = compute_roe_avg(q, q);
  
  EXPECT_DOUBLE_EQ(compute_pressure(state), compute_pressure(q));
}

TEST(RoeAvgState, PressureJac)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  RoeAvgState<Real> state = compute_roe_avg(q, q);

  auto func = [](const Vec4<Complex>& roe_avg)
  {
    return compute_pressure(RoeAvgState<Complex>{roe_avg[0], roe_avg[1], roe_avg[2], roe_avg[3]});
  };

  auto jac = [](const Vec4<Real>& roe_avg)
  {
    return compute_pressure_jac(RoeAvgState<Real>{roe_avg[0], roe_avg[1], roe_avg[2], roe_avg[3]});
  };

  test_utils::checkJacobianScalar({state.rho, state.u, state.v, state.H}, func, jac);
}

TEST(RoeAvgState, NormalVelocity)
{
  RoeAvgState<Real> state{1, 5, 10, 4};
  EXPECT_NEAR(compute_un(state, {2, 3}), 5*2 + 10*3, 1e-13);
}

TEST(RoeAvgState, ComputeUn)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  RoeAvgState<Real> state = compute_roe_avg(q, q);
  Vec2<Real> normal = {2, 3};

  auto func = [&](const Vec4<Complex>& roe_avg)
  {
    return compute_un(RoeAvgState<Complex>{roe_avg[0], roe_avg[1], roe_avg[2], roe_avg[3]}, normal);
  };

  auto jac = [&](const Vec4<Real>& roe_avg)
  {
    return compute_un_jac(RoeAvgState<Real>{roe_avg[0], roe_avg[1], roe_avg[2], roe_avg[3]}, normal);
  };

  test_utils::checkJacobianScalar({state.rho, state.u, state.v, state.H}, func, jac);
}

TEST(RoeAvgState, ComputeConservativeVars)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  RoeAvgState<Real> state = compute_roe_avg(q, q);
  Vec4<Real> q2 = compute_conservative_variables(state, RoeStateTag{});

  for (UInt i=0; i < 4; ++i)
  {
    EXPECT_DOUBLE_EQ(q[i], q2[i]);
  }
}

TEST(RoeAvgState, ComputeConservativeVarsJac)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  RoeAvgState<Real> state = compute_roe_avg(q, q);

  auto func = [&](const Vec4<Complex>& roe_avg)
  {
    return compute_conservative_variables(RoeAvgState<Complex>{roe_avg[0], roe_avg[1], roe_avg[2], roe_avg[3]}, RoeStateTag{});
  };

  auto jac = [&](const Vec4<Real>& roe_avg)
  {
    return compute_conservative_variables_jac(RoeAvgState<Real>{roe_avg[0], roe_avg[1], roe_avg[2], roe_avg[3]}, RoeStateTag{});
  };

  test_utils::checkJacobianVector({state.rho, state.u, state.v, state.H}, func, jac);
}

TEST(RoeAvgState, DifferentState)
{
  Vec4<Real> prim_varsL = {4, 6, 8, 1.0/Cv};
  Vec4<Real> prim_varsR = {16, 8, 10, 1.0/Cv};
  auto qL = compute_conservative_variables(prim_varsL, PrimitiveVarTag());
  auto qR = compute_conservative_variables(prim_varsR, PrimitiveVarTag());
  RoeAvgState state = compute_roe_avg(qL, qR);
  Real pL = compute_pressure(qL);
  Real pR = compute_pressure(qR);

  EXPECT_NEAR(state.rho, 8, 1e-13);
  EXPECT_NEAR(state.u, 44.0/6, 1e-13);
  EXPECT_NEAR(state.v, 56.0/6, 1e-13);
  EXPECT_NEAR(state.H, (204 + pL)/12 + ((1328 + pR)/24), 1e-13);
}

TEST(RoeAvgState, DifferentStateJac)
{
  Vec4<Real> prim_varsL = {4, 6, 8, 1.0/Cv};
  Vec4<Real> prim_varsR = {16, 8, 10, 1.0/Cv};
  auto qL = compute_conservative_variables(prim_varsL, PrimitiveVarTag());
  auto qR = compute_conservative_variables(prim_varsR, PrimitiveVarTag());
  Vec4<Complex> qLc = test_utils::make_complex(qL);
  Vec4<Complex> qRc = test_utils::make_complex(qR);

  auto funcL = [&](auto qL)
  {
    RoeAvgState<Complex> roe_avg = compute_roe_avg(qL, qRc);
    return Vec4<Complex>{roe_avg.rho, roe_avg.u, roe_avg.v, roe_avg.H};
  };

  auto jacL = [&](auto qL)
  {
    Matrix<Real, 4> roe_avg_dotL, roe_avg_dotR;
    RoeAvgState<Real> roe_avg = compute_roe_avg_jac(qL, qR, roe_avg_dotL, roe_avg_dotR);
    Vec4<Real> avg_state_vec {roe_avg.rho, roe_avg.u, roe_avg.v, roe_avg.H};
    return std::make_pair(avg_state_vec, roe_avg_dotL);
  };

  auto funcR = [&](auto qR)
  {
    RoeAvgState<Complex> roe_avg = compute_roe_avg(qLc, qR);
    return Vec4<Complex>{roe_avg.rho, roe_avg.u, roe_avg.v, roe_avg.H};
  };

  auto jacR = [&](auto qR)
  {
    Matrix<Real, 4> roe_avg_dotL, roe_avg_dotR;
    RoeAvgState<Real> roe_avg = compute_roe_avg_jac(qL, qR, roe_avg_dotL, roe_avg_dotR);
    Vec4<Real> avg_state_vec {roe_avg.rho, roe_avg.u, roe_avg.v, roe_avg.H};
    return std::make_pair(avg_state_vec, roe_avg_dotR);
  };  

  test_utils::checkJacobianVector(qL, funcL, jacL);
  test_utils::checkJacobianVector(qL, funcR, jacR);
}