#ifndef STRUCTURED_FINITE_VOLUME_NLSOLVERS_EXPLICIT_TIMESTEPPER_OPTS_H
#define STRUCTURED_FINITE_VOLUME_NLSOLVERS_EXPLICIT_TIMESTEPPER_OPTS_H

#include "utils/project_defs.h"

namespace structured_fv {
namespace nlsolvers {

struct ExplicitTimestepperOpts
{
  Real delta_t = 0;
  Real t_start = 0;
  Real t_end = 1;
  Real residual_tol = -1;  // stop if residual norm is less than this
                           // set to negative value to disable
  UInt itermax = std::numeric_limits<UInt>::max();

  ExplicitTimestepperOpts& operator=(const ExplicitTimestepperOpts&) = default;
};

void checkOpts(const ExplicitTimestepperOpts& opts);

}
}

#endif