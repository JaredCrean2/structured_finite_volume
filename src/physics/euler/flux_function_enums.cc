#include "flux_function_enums.h"
#include <string>
#include <stdexcept>

namespace structured_fv {
namespace euler {


std::string get_name(FluxFunction flux)
{
  switch(flux)
  {
    case FluxFunction::Roe: { return "Roe"; }
    case FluxFunction::RoeHH: { return "RoeHH"; }
    case FluxFunction::HLLE: { return "HLLE"; }
    case FluxFunction::LLF: { return "Local_Lax_Friedrich"; }
    default:  throw std::runtime_error("Unhandled FluxFunction enum");
  }
}

std::ostream& operator<<(std::ostream& os, FluxFunction flux)
{
  os << get_name(flux);
  return os;
}

}
}