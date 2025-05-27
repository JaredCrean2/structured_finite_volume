#include "gtest/gtest.h"
#include "physics/common/flux_limiters.h"
#include "physics/common//slope_limiters.h"

namespace  {

using namespace structured_fv::common;

template <typename T>
class LimiterTestFixture : public testing::Test
{};

using MyTypes = ::testing::Types<FluxLimiterMinMod,
                                 FluxLimiterSuperBee,
                                 FluxLimiterVanAlba,
                                 FluxLimiterVanLeer,
                                 SlopeLimiterMinMod,
                                 SlopeLimiterSuperBee,
                                 SlopeLimiterVanAlba,
                                 SlopeLimiterVanLeer>;

class NameGenerator
{
  public:
  template <typename T>
  static std::string GetName(int)
  {
    auto enum_val = get_enum(T());
    return get_name(enum_val);
  }
};

TYPED_TEST_SUITE(LimiterTestFixture, MyTypes, NameGenerator);


struct Foo
{};

}


TYPED_TEST(LimiterTestFixture, SpecialValues)
{
  // all TVD slope and flux limiters should have the properties
  // f(-1) = 0
  // f(0) = 0
  // f(1) = 1

  using Limiter = TypeParam;
  Limiter limiter;
  EXPECT_DOUBLE_EQ(limiter(-1.0), 0.0);
  EXPECT_DOUBLE_EQ(limiter(0.0),  0.0);
  EXPECT_DOUBLE_EQ(limiter(1.0),  1.0);
}

TEST(Limiters, IsSlopeLimiter)
{
  static_assert(IsSlopeLimiter<SlopeLimiterMinMod>);
  EXPECT_FALSE(IsSlopeLimiter<FluxLimiterMinMod>);
  EXPECT_FALSE(IsSlopeLimiter<Foo>);
}