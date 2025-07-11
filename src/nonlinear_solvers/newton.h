#ifndef STRUCTURED_FINITE_VOLUME_NLSOLVERS_NEWTON_H
#define STRUCTURED_FINITE_VOLUME_NLSOLVERS_NEWTON_H

#include "physics/physics_model.h"
#include "utils/project_defs.h"
#include "linear_system/large_matrix_factory.h"

namespace structured_fv {
namespace nlsolvers {

using DiscVectorPtr = disc::DiscVectorPtr<Real>;

// a 2 stage, 2nd order strong stability preserving (SSP) Runge Kutta method
struct NewtonOpts 
{
  Real nl_rel_tol = 1e-12;
  Real nl_abs_tol = 1e-12;
  Real linear_rel_tol = 1e-8;
  Real linear_abs_tol = 1e-13;
  UInt itermax = std::numeric_limits<UInt>::max();
  linear_system::LargeMatrixType mat_type = linear_system::LargeMatrixType::Petsc;
  std::shared_ptr<linear_system::LargeMatrixOpts> mat_opts;
};

class NewtonsMethod
{
  public:
    NewtonsMethod(const NewtonOpts& opts, PhysicsModelPtr model, disc::DiscVectorPtr<Real> sol, std::ostream& log = std::cout);

    void solve();

  private:

    void updateSolution(DiscVectorPtr delta_u, DiscVectorPtr u) const;

    bool isConverged(Real rhs_norm, std::string& reason) const;

    std::string getNonconvergenceMessage(Real rhs_norm) const;

    NewtonOpts m_opts;
    PhysicsModelPtr m_model;
    DiscVectorPtr m_sol;
    DiscVectorPtr m_rhs;
    DiscVectorPtr m_delta_u;
    linear_system::LargeMatrixPtr m_mat;
    linear_system::AssemblerBasePtr m_assembler;
    Real m_rhs_norm0 = 0;
    std::ostream& m_log;
};


}
}

#endif