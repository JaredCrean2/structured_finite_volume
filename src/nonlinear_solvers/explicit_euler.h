#ifndef STRUCTURED_FINITE_VOLUME_NLSOLVERS_EXPLICIT_EULER_H
#define STRUCTURED_FINITE_VOLUME_NLSOLVERS_EXPLICIT_EULER_H

#include "physics/physics_model.h"
#include "utils/project_defs.h"
#include <stdexcept>
#include <iostream>

namespace structured_fv {
namespace nlsolvers {

struct ExplicitEulerOpts
{
  Real delta_t = 0;
  Real t_start = 0;
  Real t_end = 1;
  Real residual_tol = -1;  // stop if residual norm is less than this
                           // set to negative value to disable
  UInt itermax = std::numeric_limits<UInt>::max();
};

void checkOpts(const ExplicitEulerOpts& opts);

void explicitEuler(const ExplicitEulerOpts& opts, PhysicsModelPtr model, disc::DiscVectorPtr<Real> sol);


}
}

#endif