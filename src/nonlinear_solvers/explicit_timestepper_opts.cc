#include "explicit_timestepper_opts.h"

namespace structured_fv {
namespace nlsolvers {

void checkOpts(const ExplicitTimestepperOpts& opts)
{
  if (opts.delta_t == 0)
    throw std::runtime_error("delta_t cannot be zero");

  if (opts.t_end <= opts.t_start)
    throw std::runtime_error("t_end must be greater than t_start");
}

}
}