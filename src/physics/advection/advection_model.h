#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_ADVECTION_MODEL_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_ADVECTION_MODEL_H

#include "disc/disc_interface.h"
#include "physics/common/slope_limiters.h"
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
  common::SlopeLimiter limiter = common::SlopeLimiter::FirstOrder;
};

template <typename T>
struct Fields
{
  Fields(const disc::StructuredDisc& disc) :
    solution(std::make_shared<disc::ElementField<T>>(disc, 1)),
    residual(std::make_shared<disc::ElementField<T>>(disc, 1))
  {}
  
  ElementFieldPtr<T> solution;
  ElementFieldPtr<T> residual;
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

    void evaluateJacobian(disc::DiscVectorPtr<Real> q, Real t, disc::DiscVectorPtr<Real> residual, linear_system::AssemblerBasePtr assembler) override;

    // computes h = dR/dq * v without forming dR/dq
    void computeJacVecProduct(disc::DiscVectorPtr<Real> q, Real t, disc::DiscVectorPtr<Real> v, disc::DiscVectorPtr<Real> h) override;

    std::shared_ptr<linear_system::SparsityPattern> getSparsityPattern() const override;

    const disc::StructuredDiscPtr& getDisc() const override { return m_disc; }

  private:

    template <typename T>
    void evaluateRhsT(Fields<T> fields, Real t);

    Fxyt& getBCFunction(UInt block_id);

    template <typename T>
    void setBCValues(ElementFieldPtr<T> solution, Real t);


    AdvectionOpts m_opts;
    StructuredDiscPtr m_disc;
    Fields<Real> m_fields_real;
    Fields<Complex> m_fields_complex;
    std::vector<Fxyt> m_bc_functions;
    Fxyt m_source_func;
};

}
}

#endif