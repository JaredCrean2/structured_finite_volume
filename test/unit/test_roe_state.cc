#include "gtest/gtest.h"
#include "physics/euler/euler_flux.h"
#include "physics/euler/roe_state.h"
#include "physics/euler/typedefs.h"

using namespace structured_fv;
using namespace structured_fv::euler;

TEST(RoeAvgState, SameState)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  RoeAvgState state = compute_roe_avg(q, q);
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

TEST(RoeAvgState, NormalVelocity)
{
  RoeAvgState state{1, 5, 10, 4};
  EXPECT_NEAR(compute_un(state, {2, 3}), 5*2 + 10*3, 1e-13);
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