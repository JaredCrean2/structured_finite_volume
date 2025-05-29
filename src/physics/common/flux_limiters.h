#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_COMMON_FLUX_LIMITERS_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_COMMON_FLUX_LIMITERS_H

#include <algorithm>
#include <string>

namespace structured_fv {
namespace common {

// for all these flux limiters, r_i = (u_i - u_i-1)/(u_i+1 - u_i)
// and are intended to be used (in the context of linear advection) as :
// u_i+1/2^L = u_i + 0.5*psi(1/r_i)(u_i - u_i-1)  (right endpoint of cell)
// u_i-1/2^R = u_i - 0.5*psi(r_i)(u_i+1 - u_i)    (left endpoint of cell)
// (note: in Murmans paper "Analysis of Slope Limiters on Irregular Grids" he 
//  defines R = 1/r)

enum class FluxLimiter
{
  FirstOrder,
  MinMod,
  SuperBee,
  VanAlba,
  VanLeer
};

std::string get_name(FluxLimiter limiter);

// sets psi = 0 everywhere, giving elementwise-constant reconstruction
class FluxLimiterFirstOrder
{
  public:

  template <typename T>
  constexpr T operator()(const T& r) const
  {
    return T(0);
  }  
};


class FluxLimiterMinMod
{
  public:

  template <typename T>
  constexpr T operator()(const T& r) const
  {
    return std::max(T(0), std::min(T(1), r));
  }
};

class FluxLimiterSuperBee
{
  public:

  template <typename T>
  constexpr T operator()(const T& r) const
  {
    T term1 = std::min(2*r, T(1));
    T term2 = std::min(r, T(2));
    return std::max(T(0), std::max(term1, term2));
  }  
};

class FluxLimiterVanAlba
{
  public:

  template <typename T>
  constexpr T operator()(const T& r) const
  {
    // the textbook form of the van Alba limiter is
    //  psi = (r*r + r)/(r*r + 1) for r >= 0 and 
    //  psi = 0 for r < 0.
    // To get psi = 0 for r < 0, modify it as
    T num = r + 1;
    return 0.5*(r*num + std::abs(r)*num)/(r*r + 1);
  }  
};

class FluxLimiterVanLeer
{
  public:

  template <typename T>
  constexpr T operator()(const T& r) const
  {
    return (r + std::abs(r))/(1 + std::abs(r));
  }  
};

constexpr FluxLimiter get_enum(FluxLimiterFirstOrder)
{
  return FluxLimiter::FirstOrder;
}

constexpr FluxLimiter get_enum(FluxLimiterMinMod)
{
  return FluxLimiter::MinMod;
}

constexpr FluxLimiter get_enum(FluxLimiterSuperBee)
{
  return FluxLimiter::SuperBee;
}

constexpr FluxLimiter get_enum(FluxLimiterVanAlba)
{
  return FluxLimiter::VanAlba;
}

constexpr FluxLimiter get_enum(FluxLimiterVanLeer)
{
  return FluxLimiter::VanLeer;
}

}
}

#endif