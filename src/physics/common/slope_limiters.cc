#include "slope_limiters.h"
#include <stdexcept>

namespace structured_fv {
namespace common {

std::string get_name(SlopeLimiter limiter)
{
  switch (limiter)
  {
    case SlopeLimiter::FirstOrder: { return "Slope_Limiter_FirstOrder"; }
    case SlopeLimiter::MinMod:     { return "Slope_Limiter_MinMod"; }
    case SlopeLimiter::SuperBee:   { return "Slope_Limiter_SuperMod"; }
    case SlopeLimiter::VanAlba:    { return "Slope_Limiter_VanAlba"; }
    case SlopeLimiter::VanLeer:    { return "Slope_Limiter_VanLeer"; }
    default:
    {
      throw std::runtime_error("unhandled FluxLimiter enum");
    }
  }
}

}
}