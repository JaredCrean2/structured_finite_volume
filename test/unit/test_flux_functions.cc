#include "gtest/gtest.h"
#include "physics/euler/flux_function_enums.h"
#include "physics/euler/roe_flux.h"
#include "physics/euler/roe_hh_flux.h"
#include "physics/euler/hlle_flux.h"
#include "physics/euler/lax_friedrich_flux.h"
#include "physics/euler/hllc_flux.h"


namespace {

using namespace structured_fv;
using namespace structured_fv::euler;


template <typename T>
class FluxFunctionTester : public testing::Test
{};

using FluxFunctions = ::testing::Types<euler::RoeFlux,
                                       euler::RoeHHFlux,
                                       euler::HLLEFlux,
                                       euler::LaxFriedrichFlux,
                                       euler::HLLCFlux>;

class NameGenerator
{
  public:
    template <typename T>
    static std::string GetName(int)
    {
      if constexpr (std::is_same_v<T, euler::RoeFlux>)
        return euler::get_name(euler::FluxFunction::Roe);
      if constexpr (std::is_same_v<T, euler::RoeHHFlux>)
        return euler::get_name(euler::FluxFunction::RoeHH);
      if constexpr (std::is_same_v<T, euler::HLLEFlux>)
        return euler::get_name(euler::FluxFunction::HLLE);
      if constexpr (std::is_same_v<T, euler::LaxFriedrichFlux>)
        return euler::get_name(euler::FluxFunction::LLF);
      if constexpr (std::is_same_v<T, euler::HLLCFlux>)
        return euler::get_name(euler::FluxFunction::HLLC);
      else
        throw std::runtime_error("unhandled enum");       
    }
};

TYPED_TEST_SUITE(FluxFunctionTester, FluxFunctions, NameGenerator);
}


TYPED_TEST(FluxFunctionTester, SameState)
{
  using FluxFunc = TypeParam;

  Vec4<Real> prim_vars = {2, 5, 10, 300};
  auto q = compute_conservative_variables(prim_vars, PrimitiveVarTag());
  Vec2<Real> normal = {2, 3};

  FluxFunc flux;
  Vec4<Real> flux_val = flux(q, q, normal);
  Vec4<Real> flux_expected = compute_euler_flux(q, normal);

  for (UInt i=0; i < DofsPerCell; ++i)
  {
    std::cout << "\ni = " << i << std::endl;
    EXPECT_NEAR(flux_val[i], flux_expected[i], 1e-8);
  }
}

TYPED_TEST(FluxFunctionTester, Supersonic)
{
  using FluxFunc = TypeParam;

  if (std::is_same_v<FluxFunc, euler::LaxFriedrichFlux>)
    GTEST_SKIP();

  Vec4<Real> prim_varsL = {2, 400, 500, 300};
  Vec4<Real> prim_varsR = {2, 500, 700, 300};
  auto qL = compute_conservative_variables(prim_varsL, PrimitiveVarTag());
  auto qR = compute_conservative_variables(prim_varsR, PrimitiveVarTag());
  Vec2<Real> normalR = {1.0, 0.0};
  Vec2<Real> normalL = {-1.0, 0.0};

  FluxFunc flux;

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

TYPED_TEST(FluxFunctionTester, Symmetry)
{
  using FluxFunc = TypeParam;

  Vec4<Real> qL = {1.819813647055912, 761.5114870354085, -12.50056587286352, 542762.079360307};
  Vec4<Real> qR = {1.847883122092081, 771.4842523148619, 2.673569354817001, 552845.7433148436};
  Vec2<Real> normal = {0, 0.3333333333333334};
  FluxFunc flux;

  auto flux1 = flux(qL, qR, normal);
  auto flux2 = flux(qR, qL, -normal);
  for (int i=0; i < 4; ++i)
    EXPECT_NEAR(flux1[i], -flux2[i], 1e-8);
}

