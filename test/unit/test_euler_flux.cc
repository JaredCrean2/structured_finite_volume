#include "gtest/gtest.h"
#include "physics/euler/euler_flux.h"

using namespace structured_fv;
using namespace structured_fv::euler;

TEST(EulerFlux, Pressure)
{
  // standard atmosphere at 0 altitude
  Vec4<Real> prim_vars = {1.225, 0, 0, 288.16};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  EXPECT_NEAR(compute_pressure(q), 1.01325E5, 100.0);
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