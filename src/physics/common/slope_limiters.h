#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_COMMON_SLOPE_LIMITERS_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_COMMON_SLOPE_LIMITERS_H

#include "flux_limiters.h"
#include <algorithm>
#include <string>
#include <iostream>

namespace structured_fv {
namespace common {

// Following Murmans paper "Analysis of Slope Limiters on Irregular Grids",
// slope limiters are used as:
// u_i+1/2^L = u_i + 0.5*phi(R)(u_i+1 - u_i-1)/2
// u_i-1/2^R = u_i - 0.5*phi(R)(u_i+1 - u_i-1)/2
// Note that Murman defines R = 1/r, where r is Sweby's r, but
// phi(R) = phi(1/R), so it doesnt matter


enum class SlopeLimiter
{
  FirstOrder,
  MinMod,
  SuperBee,
  VanAlba,
  VanLeer
};

std::string get_name(SlopeLimiter limiter);

inline std::ostream& operator<<(std::ostream& os, SlopeLimiter limiter)
{
  os << get_name(limiter);
  return os;
}

class SlopeLimiterBase
{};

template <typename SlopeLimiterType>
constexpr bool IsSlopeLimiter = std::is_base_of_v<SlopeLimiterBase, SlopeLimiterType>;

template <typename FluxLimiterType>
class SlopeLimiterAdaptor : public SlopeLimiterBase
{
  public:

  template <typename T>
  constexpr T operator()(const T& r) const
  {
    // all TVD flux limiters have phi(r) = 0 for r < 0
    // The conversion factor is 2/(1+r), so we need to guard against r=-1
    T factor = T(2)/(1 + std::max(r, T(0)));
    return factor*m_flux_limiter(r);
  }

  private:
    FluxLimiterType m_flux_limiter;
};

using SlopeLimiterFirstOrder = SlopeLimiterAdaptor<FluxLimiterFirstOrder>;
using SlopeLimiterMinMod     = SlopeLimiterAdaptor<FluxLimiterMinMod>;
using SlopeLimiterSuperBee   = SlopeLimiterAdaptor<FluxLimiterSuperBee>;
using SlopeLimiterVanAlba    = SlopeLimiterAdaptor<FluxLimiterVanAlba>;
using SlopeLimiterVanLeer    = SlopeLimiterAdaptor<FluxLimiterVanLeer>;

constexpr SlopeLimiter get_enum(SlopeLimiterFirstOrder)
{
  return SlopeLimiter::FirstOrder;
}

constexpr SlopeLimiter get_enum(SlopeLimiterMinMod)
{
  return SlopeLimiter::MinMod;
}

constexpr SlopeLimiter get_enum(SlopeLimiterSuperBee)
{
  return SlopeLimiter::SuperBee;
}

constexpr SlopeLimiter get_enum(SlopeLimiterVanAlba)
{
  return SlopeLimiter::VanAlba;
}

constexpr SlopeLimiter get_enum(SlopeLimiterVanLeer)
{
  return SlopeLimiter::VanLeer;
}


template <typename FluxLimiter>
std::ostream& operator<<(std::ostream& os, const SlopeLimiterAdaptor<FluxLimiter>& limiter)
{
  os << get_name(get_enum(limiter));
  return os;
}
}
}

#endif