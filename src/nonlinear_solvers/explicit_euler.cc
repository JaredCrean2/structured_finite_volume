#include "explicit_euler.h"

namespace structured_fv {
namespace nlsolvers {

void checkOpts(const ExplicitEulerOpts& opts)
{
  if (opts.delta_t == 0)
    throw std::runtime_error("delta_t cannot be zero");

  if (opts.t_end <= opts.t_start)
    throw std::runtime_error("t_end must be greater than t_start");
}


void explicitEuler(const ExplicitEulerOpts& opts, PhysicsModelPtr model, disc::DiscVectorPtr<Real> sol)
{
  checkOpts(opts);

  auto dudt = std::make_shared<disc::DiscVector<Real>>(model->getDisc(), "residual");

  Real t = opts.t_start;
  bool converged = false;
  UInt iter = 0;
  while (t < opts.t_end && !converged && iter != opts.itermax)
  {
    std::cout << "\niter " << iter << std::endl;
    model->evaluateRhs(sol, t, dudt);

    
    if (opts.residual_tol >= 0)
    {
      Real norm = model->computeRhsNorm(dudt);
      converged = norm < opts.residual_tol;
      std::cout << "norm = " << norm << std::endl;
    }

    if (!converged)
    {
      for (GlobalDof i=0; i < sol->size(); ++i)
        (*sol)(i) += opts.delta_t * (*dudt)(i);

      t = std::min(t + opts.delta_t, opts.t_end);
    }

    iter++;
  }

  if (t < opts.t_end && (opts.residual_tol >= 0 && !converged))
    throw std::runtime_error("ExplicitEuler failed");
}

}
}