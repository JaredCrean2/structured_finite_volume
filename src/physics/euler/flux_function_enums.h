#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_FLUX_ENUMS_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_FLUX_ENUMS_H

#include <iosfwd>

namespace structured_fv {
namespace euler {

enum class FluxFunction
{
  Roe,
  RoeHH,
  LLF,
  HLLE,
  HLLC,
};

std::string get_name(FluxFunction flux);

std::ostream& operator<<(std::ostream& os, FluxFunction flux);

}
}

#endif