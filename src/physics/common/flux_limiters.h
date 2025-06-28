#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_COMMON_FLUX_LIMITERS_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_COMMON_FLUX_LIMITERS_H

#include <algorithm>
#include <string>
#include "utils/math.h"

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

  template <typename T>
  constexpr std::pair<T, T> operator()(const T& r, const T& r_dot) const
  {
    return {T(0), 0};
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

  template <typename T>
  constexpr std::pair<T, T> operator()(const T& r, const T& r_dot) const
  {
    T t1 = std::min(T(1), r);
    T t1_dot = T(1) < r ? 0.0 : r_dot;

    T psi = std::max(T(0), t1);
    T psi_dot = t1 < T(0) ? 0.0 : t1_dot;

    return {psi, psi_dot};
  }  
};

class FluxLimiterSuperBee
{
  public:

  template <typename T>
  constexpr T operator()(const T& r) const
  {
    T term1 = std::min(2.0*r, T(1));
    T term2 = std::min(r, T(2));
    return std::max(T(0), std::max(term1, term2));
  }

  template <typename T>
  constexpr std::pair<T, T> operator()(const T& r, const T& r_dot) const
  {
    T term1 = std::min(2*r, T(1));
    T term1_dot = 2*r < T(1) ? 2*r_dot : 0.0;

    T term2 = std::min(r, T(2));
    T term2_dot = r < T(2) ? r_dot : 0.0;

    T maxt = std::max(term1, term2);
    T maxt_dot = term1 > term2 ? term1_dot : term2_dot;

    T psi = std::max(T(0), maxt);
    T psi_dot = T(0) > maxt ? 0.0 : maxt_dot;

    return {psi, psi_dot};
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
    T num = r + T(1);
    return 0.5*(r*num + smoothAbs(r)*num)/(r*r + T(1));
  }

  template <typename T>
  constexpr std::pair<T, T> operator()(const T& r, const T& r_dot) const
  {
    T num = r + T(1);
    T num_dot = r_dot;

    T t1 = r*num;
    T t1_dot = r_dot*num + num_dot*r;

    auto [r_abs, r_abs_dot] = smoothAbs_dot(r, r_dot);
    T t2 = r_abs*num;
    T t2_dot = r_abs*num_dot + num*r_abs_dot;
  
    T t3 = t1 + t2;
    T t3_dot = t1_dot + t2_dot;

    T den = r*r + T(1);
    T den_dot = 2*r*r_dot;

    T psi = 0.5*t3/den;
    T psi_dot = 0.5*(t3_dot/(r*r+1) - t3*den_dot/(den*den));

    return {psi, psi_dot};
  }
};

class FluxLimiterVanLeer
{
  public:

  template <typename T>
  constexpr T operator()(const T& r) const
  {
    T abs_r = smoothAbs(r);
    return (r + abs_r)/(T(1) + abs_r);
  }

  template <typename T>
  constexpr std::pair<T, T> operator()(const T& r, const T& r_dot) const
  {
    auto [abs_r, abs_r_dot] = smoothAbs_dot(r, r_dot);
    
    T den = T(1) + abs_r;
    T den_dot = abs_r_dot;

    T psi = (r + abs_r)/den;
    T psi_dot = (r_dot + abs_r_dot)/den - (r + abs_r)*den_dot/(den*den);

    return {psi, psi_dot};
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