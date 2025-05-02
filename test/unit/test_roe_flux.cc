#include "gtest/gtest.h"
#include "physics/euler/roe_flux.h"
#include "physics/euler/typedefs.h"

using namespace structured_fv;
using namespace structured_fv::euler;

TEST(RoeFlux, SameState)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  Vec2<Real> normal = {2, 3};

  RoeFlux flux(0.0);
  Vec4<Real> flux_val = flux(q, q, normal);
  Vec4<Real> flux_expected = compute_euler_flux(q, normal);

  for (UInt i=0; i < DofsPerCell; ++i)
    EXPECT_NEAR(flux_val[i], flux_expected[i], 1e-8);
}

TEST(RoeFlux, Supersonic)
{
  Vec4<Real> prim_varsL = {2, 400, 500, 300};
  Vec4<Real> prim_varsR = {2, 500, 700, 300};
  auto qL = compute_conservative_variables(prim_varsL, PrimitiveVarTag());
  auto qR = compute_conservative_variables(prim_varsR, PrimitiveVarTag());
  Vec2<Real> normalR = {1.0, 0.0};
  Vec2<Real> normalL = {-1.0, 0.0};

  RoeFlux flux(0.0);

  {
    Vec4<Real> flux_val = flux(qL, qR, normalR);
    Vec4<Real> flux_expected = compute_euler_flux(qL, normalR);
    for (UInt i=0; i < DofsPerCell; ++i)
      EXPECT_NEAR(flux_val[i], flux_expected[i], 1e-8);
  }


  {
    Vec4<Real> flux_val = flux(qL, qR, normalL);
    Vec4<Real> flux_expected = compute_euler_flux(qR, normalL);
    for (UInt i=0; i < DofsPerCell; ++i)
      EXPECT_NEAR(flux_val[i], flux_expected[i], 1e-8);
  }
}