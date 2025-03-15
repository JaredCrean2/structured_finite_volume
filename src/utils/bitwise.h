#ifndef STRUCTURED_FINITE_VOLUME_UTILS_BITWISE_H
#define STRUCTURED_FINITE_VOLUME_UTILS_BITWISE_H

#include "project_defs.h"

namespace structured_fv {

constexpr UInt UIntMax = std::numeric_limits<UInt>::max();

// returns a UInt with all bits set to 1 if a == b
constexpr UInt mask_if_equal(UInt a, UInt b)
{
  return (a != b) + UIntMax;
}

constexpr UInt mask_if_notequal(UInt a, UInt b)
{
  return ~mask_if_equal(a, b);
}

}

#endif