#ifndef STRUCTURED_FINITE_VOLUME_NLSOLVERS_EXPLICIT_EULER_H
#define STRUCTURED_FINITE_VOLUME_NLSOLVERS_EXPLICIT_EULER_H

#include "physics/physics_model.h"
#include "utils/project_defs.h"
#include "explicit_timestepper_opts.h"

namespace structured_fv {
namespace nlsolvers {

struct ExplicitEulerOpts : public ExplicitTimestepperOpts
{
  using ExplicitTimestepperOpts::operator=;
};

void explicitEuler(const ExplicitEulerOpts& opts, PhysicsModelPtr model, disc::DiscVectorPtr<Real> sol);

}
}

#endif