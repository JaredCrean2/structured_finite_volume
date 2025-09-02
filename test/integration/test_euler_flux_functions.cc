#include "gtest/gtest.h"
#include "physics/euler/flux_function_enums.h"
#include "physics/euler/roe_flux.h"
#include "physics/euler/roe_hh_flux.h"
#include "physics/euler/hlle_flux.h"
#include "physics/euler/lax_friedrich_flux.h"
#include "physics/euler/hllc_flux.h"
#include "utils/cpu_timer.h"


namespace {

using namespace structured_fv;
using namespace structured_fv::euler;



template <typename T>
class FluxFunctionJacTester : public testing::Test
{};

using FluxFunctionsJac = ::testing::Types<euler::LaxFriedrichFlux,
                                          euler::HLLCFlux,
                                          euler::HLLEFlux,
                                          euler::RoeFlux,
                                          euler::RoeHHFlux>;

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


TYPED_TEST_SUITE(FluxFunctionJacTester, FluxFunctionsJac, NameGenerator);

}

TYPED_TEST(FluxFunctionJacTester, Performance)
{
  using FluxFunc = TypeParam;
  using Dual8 = Dual<Real, 8>;

  constexpr UInt NumIters = 10000;

  FluxFunc flux;
  std::vector<std::pair<Vec4<Real>, Vec4<Real>>> states =
  {
     {{2,  400.0/2,  500.0/3, 300}, {2,  500.0/2,  700.0/3, 300}},
     {{2, -400.0/2, -500.0/3, 300}, {2, -500.0/2, -700.0/3, 300}},
     {{2,   40.0/2,   50.0/3, 300}, {2,   50.0/2,   70.0/3, 300}},  // HLLC: uses SL and SR
     {{2,   40.0/2,   50.0/3, 300}, {2,   50.0/2,   70.0/3, 600}},  // HLLC: uses SL_avg and SR, also s_star < 0
     {{2,   40.0/2,   50.0/3, 600}, {2,   50.0/2,   70.0/3, 300}},  // HLLC: uses sR_avg
     {{2,  500.0/2,   50.0/3, 700}, {2,  500.0/2,   50.0/3, 300}},  // RoeHH: trip entropy fix
     {{2,   40.0/2,   50.0/3, 300}, {2,   40.0/2,   50.0/3, 300}},
     {{2,   40.0/2,   50.0/3, 300}, {2,   40,       10.0/3, 300}},  // LLF: max wave speed is same on both sides of interface
  };
  Vec2<Real> normal = {2.0, 3.0};


  std::vector<std::pair<Vec4<Real>, Vec4<Real>>> states_conservative;
  for (auto [prim_varsL, prim_varsR] : states)
  {
    auto qL = compute_conservative_variables(prim_varsL, PrimitiveVarTag());
    auto qR = compute_conservative_variables(prim_varsR, PrimitiveVarTag());
    states_conservative.push_back(std::make_pair(qL, qR));
  }

  std::string desc_hand_coded = NameGenerator::GetName<FluxFunc>(0) + " hand coded jacobian";
  {
    ScopedCPUTimer timer(desc_hand_coded);
    for (UInt i=0; i < NumIters; ++i)
    {
      for (const auto& [qL, qR] : states)
      {
        Matrix<Real, 4> flux_dotL, flux_dotR;
        flux(qL, qR, normal, flux_dotL, flux_dotR);      
      }
    }
  }

  std::string desc_dual = NameGenerator::GetName<FluxFunc>(0) + " Dual jacobian";
  {
    ScopedCPUTimer timer(desc_dual);
    for (UInt i=0; i < NumIters; ++i)
    {
      for (const auto& [qL, qR] : states)
      {
        Vec4<Dual8> qL_dual, qR_dual;
        for (UInt j=0; j < 4; ++j)
        {
          qL_dual[j] = Dual8(qL[j]);
          qL_dual[j](j) = 1;

          qR_dual[j] = Dual8(qR[j]);
          qR_dual[j](j + 4) = 1;
        }

        Vec4<Dual8> flux_dual = flux(qL_dual, qR_dual, normal);

        Matrix<Real, 4> flux_dotL, flux_dotR;
        for (UInt j=0; j < 4; ++j)
          for (UInt k=0; k < 4; ++k)
          {
            flux_dotL(j, k) = flux_dual[j](k);
            flux_dotR(j, k) = flux_dual[j](k+4);
          }
      }
    }
  } 
}