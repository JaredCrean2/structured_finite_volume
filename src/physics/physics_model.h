#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_PHYSICS_MODEL_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_PHYSICS_MODEL_H

#include "disc/disc_vector.h"
#include "disc/discretization.h"
#include "linear_system/assembler_base.h"
#include "linear_system/sparsity_pattern.h"


namespace structured_fv {

class PhysicsModel
{
  public:

    virtual ~PhysicsModel() = default;
    
    // evaluate the right hand side of the equation:
    // dq/dt = R(q, t)
    // source terms are defined as being on the right hand side:
    // dq/dt = R(q, t) + S(x, t)
    // where S(x, t) is the source term
    virtual void evaluateRhs(disc::DiscVectorPtr<Real> q, Real t, disc::DiscVectorPtr<Real> residual) = 0;

    virtual void evaluateJacobian(disc::DiscVectorPtr<Real> q, Real t, disc::DiscVectorPtr<Real> residual, linear_system::AssemblerBasePtr assembler) = 0;

    // computes h = dR/dq * v without forming dR/dq
    virtual void computeJacVecProduct(disc::DiscVectorPtr<Real> q, Real t, disc::DiscVectorPtr<Real> v, disc::DiscVectorPtr<Real> h) = 0;

    virtual Real computeRhsNorm(disc::DiscVectorPtr<Real> residual);

    virtual std::shared_ptr<linear_system::SparsityPattern> getSparsityPattern() const = 0;

    virtual const disc::StructuredDiscPtr& getDisc() const = 0;

};

using PhysicsModelPtr = std::shared_ptr<PhysicsModel>;


}

#endif