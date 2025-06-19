#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_NUMERIECAL_FLUX_BASE_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_NUMERIECAL_FLUX_BASE_H

#include <type_traits>

namespace structured_fv {
namespace euler {

class NumericalFlux
{};

template <typename Flux>
constexpr bool IsNumericalFlux = std::is_base_of_v<NumericalFlux, Flux>;

}
}

#endif