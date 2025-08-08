#include "newton.h"

namespace structured_fv {
namespace nlsolvers {

NewtonsMethod::NewtonsMethod(const NewtonOpts& opts, PhysicsModelPtr model, disc::DiscVectorPtr<Real> sol,
              std::ostream& log) :
  m_opts(opts),
  m_model(model),
  m_sol(sol),
  m_rhs(std::make_shared<disc::DiscVector<Real>>(model->getDisc(), "residual")),
  m_delta_u(std::make_shared<disc::DiscVector<Real>>(model->getDisc(), "delta_u")),
  m_mat(linear_system::largeMatrixFactory(opts.mat_type, opts.mat_opts, model->getSparsityPattern())),
  m_assembler(m_mat->getAssembler(model->getDisc())),
  m_log(log)
{}

void NewtonsMethod::solve()
{
  std::string convergence_reason;

  m_model->evaluateRhs(m_sol, 0.0, m_rhs);
  Real rhs_norm = m_model->computeRhsNorm(m_rhs);
  m_rhs_norm0 = rhs_norm;

  if (isConverged(rhs_norm, convergence_reason))
  {
    m_log << "initial solution satisfied convergence criteria:\n"
          << convergence_reason << std::endl;
    return;
  }

  for (UInt iter=0; iter < m_opts.itermax; ++iter)
  {
    m_mat->zeroMatrix();
    m_model->evaluateJacobian(m_sol, 0.0, m_rhs, m_assembler);

    m_mat->finishMatrixAssembly();
    m_mat->factor();
    m_mat->solve(m_rhs->getData(), m_delta_u->getData());
    updateSolution(m_delta_u, m_sol);

    m_model->evaluateRhs(m_sol, 0.0, m_rhs);
    rhs_norm = m_model->computeRhsNorm(m_rhs);

    if (isConverged(rhs_norm, convergence_reason))
    {
      m_log << "iteration " << iter << " converged:\n"
            << convergence_reason << std::endl;

      return;
    }
  }

  throw std::runtime_error(getNonconvergenceMessage(rhs_norm));
}


void NewtonsMethod::updateSolution(DiscVectorPtr delta_u, DiscVectorPtr u) const
{
  for (UInt i=0; i < u->size(); ++i)
    (*u)(i) -= (*delta_u)(i);
}

bool NewtonsMethod::isConverged(Real rhs_norm, std::string& reason) const
{
  std::stringstream ss;
  ss << std::scientific;

  bool is_converged = false;
  Real rel_residual = rhs_norm/m_rhs_norm0 ;
  if (rhs_norm <= m_opts.nl_abs_tol)
  {
    is_converged = true;
    ss << "abs tolerance satisfied: " << rhs_norm << " <= " << m_opts.nl_abs_tol << "\n";
  }

  if (rel_residual < m_opts.nl_rel_tol)
  {
    is_converged = true;
    ss << "rel tol satisfied: " << rel_residual << " <= " << m_opts.nl_rel_tol << "\n";
  }

  reason = ss.str();
  return is_converged;
}

std::string NewtonsMethod::getNonconvergenceMessage(Real rhs_norm) const
{
  std::stringstream ss;
  ss << std::scientific;
  ss << "Newtons method failed to converge in " << m_opts.itermax << " iterations:\n"
      << "final abs residual = " << rhs_norm << " > " << m_opts.nl_abs_tol << "\n"
      << "final rel residual = " << rhs_norm/m_rhs_norm0 << " > " << m_opts.nl_rel_tol << "\n";

  return ss.str();
}

}
}