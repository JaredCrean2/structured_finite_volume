#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_ADVECTION_MODEL_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_ADVECTION_MODEL_H

#include "disc/disc_interface.h"
#include "physics/physics_model.h"
#include "disc/discretization.h"
#include "utils/math.h"

namespace structured_fv {
namespace advection {

using disc::ElementFieldPtr;
using disc::StructuredDiscPtr;
using disc::DiscVectorPtr;
using disc::StructuredBlock;
using disc::StructuredBlockInterface;

// f(x, y, t)
using Fxyt = std::function<Real(Real, Real, Real)>;


struct AdvectionOpts
{
  Vec2<Real> adv_velocity{0, 0};
};

class AdvectionModel : public PhysicsModel
{
  public:
    AdvectionModel(const AdvectionOpts& opts, StructuredDiscPtr disc,
                   const std::vector<Fxyt>& bc_functions,
                   Fxyt source_func);

    // evaluate the right hand side of the equation:
    // dq/dt = R(q, t)
    void evaluateRhs(DiscVectorPtr<Real> q, Real t, DiscVectorPtr<Real> residual) override;

    Real computeRhsNorm(disc::DiscVectorPtr<Real> residual) override;


    const disc::StructuredDiscPtr& getDisc() const override { return m_disc; }

  private:

    Fxyt& getBCFunction(UInt block_id);

    void setBCValues(ElementFieldPtr<Real> solution, Real t);

    template <typename Flux, typename Tag>
    void evaluateInterfaceTerms(const ElementFieldPtr<Real>& solution, Real t,
                                Flux& flux, Tag tag, ElementFieldPtr<Real> resdual);

    void evaluateSourceTerm(Real t, ElementFieldPtr<Real> residual);

    AdvectionOpts m_opts;
    StructuredDiscPtr m_disc;
    ElementFieldPtr<Real> m_solution;
    ElementFieldPtr<Real> m_residual;
    std::vector<Fxyt> m_bc_functions;
    Fxyt m_source_func;
};


}
}

#endif