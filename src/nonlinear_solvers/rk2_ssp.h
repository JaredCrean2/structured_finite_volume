#ifndef STRUCTURED_FINITE_VOLUME_NLSOLVERS_RK2_SSP_H
#define STRUCTURED_FINITE_VOLUME_NLSOLVERS_RK2_SSP_H

#include "physics/physics_model.h"
#include "utils/project_defs.h"
#include "explicit_timestepper_opts.h"

namespace structured_fv {
namespace nlsolvers {

// a 2 stage, 2nd order strong stability preserving (SSP) Runge Kutta method
struct RK2SSPOpts : public ExplicitTimestepperOpts
{
  using ExplicitTimestepperOpts::operator=;
};

void rk2ssp(const RK2SSPOpts& opts, PhysicsModelPtr model, disc::DiscVectorPtr<Real> sol);

}
}

#endif