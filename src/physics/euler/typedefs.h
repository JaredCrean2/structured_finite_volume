#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_TYPEDEFS_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_TYPEDEFS_H

#include "disc/disc_vector.h"
#include "disc/discretization.h"
#include "utils/math.h"

namespace structured_fv {
namespace euler {

using disc::ElementFieldPtr;
using disc::StructuredDiscPtr;
using disc::DiscVectorPtr;
using disc::StructuredBlock;
using disc::StructuredBlockInterface;

constexpr UInt DofsPerCell = 4;

template <typename T>
using Vec4 = std::array<T, 4>;

// f(x, y, t)
using Fxyt = std::function<Vec4<Real>(Real, Real, Real)>;

constexpr Real Gamma = 1.4;
constexpr Real Gamma_m1 = Gamma - 1;
constexpr Real Cv = 718.0;  // J/kg*K
constexpr Real Cp = Gamma * Cv;
constexpr Real R = Cp - Cv;

}
}


#endif