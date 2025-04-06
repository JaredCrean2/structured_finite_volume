#include "gtest/gtest.h"
#include "physics/euler/hlle_flux.h"
#include "physics/euler/typedefs.h"

using namespace structured_fv;
using namespace structured_fv::euler;

TEST(HLLEFlux, SameState)
{
  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  Vec2<Real> normal = {2, 3};

  HLLEFlux flux;
  Vec4<Real> flux_val = flux(q, q, normal);
  Vec4<Real> flux_expected = compute_euler_flux(q, normal);

  for (UInt i=0; i < DofsPerCell; ++i)
    EXPECT_NEAR(flux_val[i], flux_expected[i], 1e-8);
}