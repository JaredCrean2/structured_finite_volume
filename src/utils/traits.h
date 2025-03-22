#ifndef STRUCTURED_FINITE_VOLUME_UTILS_TRAITS_H
#define STRUCTURED_FINITE_VOLUME_UTILS_TRAITS_H

#include "project_defs.h"

namespace structured_fv {

// used for SFINAE when Func is a function func(Real, Real)
template <typename Func>
using IsFuncXY_t = std::enable_if_t<std::is_invocable<Func, Real, Real>::value, bool>;

}

#endif