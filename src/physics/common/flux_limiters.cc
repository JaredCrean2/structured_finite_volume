#include "flux_limiters.h"
#include <stdexcept>

namespace structured_fv {
namespace common {

std::string get_name(FluxLimiter limiter)
{
  switch (limiter)
  {
    case FluxLimiter::FirstOrder: { return "Flux_Limiter_FirstOrder"; }
    case FluxLimiter::MinMod:     { return "Flux_Limiter_MinMod"; }
    case FluxLimiter::SuperBee:   { return "Flux_Limiter_SuperBee"; }
    case FluxLimiter::VanAlba:    { return "Flux_Limiter_VanAlba"; }
    case FluxLimiter::VanLeer:    { return "Flux_Limiter_VanLeer"; }
    default:
    {
      throw std::runtime_error("unhandled FluxLimiter enum");
    }
  }
}

}
}