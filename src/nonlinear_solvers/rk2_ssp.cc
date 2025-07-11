#include "rk2_ssp.h"

namespace structured_fv {
namespace nlsolvers {

void rk2ssp(const RK2SSPOpts& opts, PhysicsModelPtr model, disc::DiscVectorPtr<Real> sol_ptr)
{
  checkOpts(opts);

  auto q_star_ptr = std::make_shared<disc::DiscVector<Real>>(model->getDisc(), "q_star");
  auto residual_ptr = std::make_shared<disc::DiscVector<Real>>(model->getDisc(), "residual");

  auto& sol      = *sol_ptr;
  auto& q_star   = *q_star_ptr;
  auto& residual = *residual_ptr;

  Real t = opts.t_start;
  bool converged = false;
  UInt iter = 0;
  while (t < opts.t_end && iter != opts.itermax)
  {
    std::cout << "\niter " << iter << std::endl;
    Real delta_t = std::min(opts.delta_t, opts.t_end - t);
    model->evaluateRhs(sol_ptr, t, residual_ptr);

    if (opts.residual_tol >= 0)
    {
      Real norm = model->computeRhsNorm(residual_ptr);
      converged = norm < opts.residual_tol;
      std::cout << "norm = " << norm << std::endl;
      if (converged)
        break;
    }

    Real max_q_star = 0.0;
    for (UInt i=0; i < q_star.size(); ++i)
    {
      q_star(i) = sol(i) + delta_t * residual(i);
      max_q_star = std::max(max_q_star, std::abs(q_star(i)));
    }
    std::cout << "max_q_star = " << max_q_star << std::endl;

    model->evaluateRhs(q_star_ptr, t + delta_t, residual_ptr);

    Real max_delta_u = 0.0;
    for (UInt i=0; i < sol.size(); ++i)
    {
      Real q_double_star = q_star(i) + delta_t * residual(i);
      sol(i) = 0.5*(sol(i) + q_double_star);
      max_delta_u = std::max(2*std::abs(q_double_star), max_delta_u);
    }
    std::cout << "max_delta_u = " << max_delta_u << std::endl;

    t += delta_t;
    iter++;
  }

  if (t < opts.t_end && (opts.residual_tol >= 0 && !converged))
    throw std::runtime_error("RK2 SSP failed");
}

}
}