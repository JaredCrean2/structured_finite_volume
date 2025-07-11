#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_PHYSICS_MODEL_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_PHYSICS_MODEL_H

#include "disc/disc_block.h"
#include "disc/disc_vector.h"
#include "disc/discretization.h"
#include "linear_system/assembler_base.h"
#include "linear_system/sparsity_pattern.h"
#include "disc/elem_field.h"


#include <iostream>

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


template <typename T>
void vecToField(disc::StructuredDiscPtr disc, const disc::DiscVectorPtr<T> vec, disc::ElementFieldPtr<T> field, bool update_ghosts=true)
{
  const UInt vals_per_cell = field->getNumValsPerElement();
  for (UInt block_id=0; block_id < disc->getNumRegularBlocks(); ++block_id)
  {
    const disc::StructuredBlock& block = disc->getBlock(block_id);
    auto& field_data = field->getData(block_id);
    auto& dof_nums = disc->getDofNumbering()->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        for (UInt k=0; k < vals_per_cell; ++k)
          field_data(i, j, k) = (*vec)(dof_nums(i, j, k));
  }

  if (update_ghosts)
    field->updateGhostValues();
}

// updates field for computing jacobian-vector products, where
// vec is the state the jacobian will be evaluated at
// and vec_dot is the vector the jacobian will be multiplied by
void vecToFieldDot(disc::StructuredDiscPtr disc, 
                   const disc::DiscVectorPtr<Real> vec, 
                   const disc::DiscVectorPtr<Real> vec_dot,
                   disc::ElementFieldPtr<Complex> field,
                   Real h=1e-40);

// copies owned values from the field into corresponding entries in the vector.
// Does not sum non-owned values
template <typename T>
void fieldToVec(disc::StructuredDiscPtr disc, const disc::ElementFieldPtr<T> field, disc::DiscVectorPtr<T> vec)
{
  const UInt vals_per_cell = field->getNumValsPerElement();
  for (UInt block_id=0; block_id < disc->getNumRegularBlocks(); ++block_id)
  {
    const disc::StructuredBlock& block = disc->getBlock(block_id);
    auto& field_data = field->getData(block_id);
    auto& dof_nums = disc->getDofNumbering()->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        for (UInt k=0; k < vals_per_cell; ++k)
          (*vec)(dof_nums(i, j, k)) = field_data(i, j, k); 
  }  
}

// given a field with the derivative stored in the complex part of each value,
// computes the vector form of the jacobian vector product, ie.
// vec_dot = dR/dq * v
void fieldToVecDot(disc::StructuredDiscPtr disc,
                   disc::ElementFieldPtr<Complex> field,
                   disc::DiscVectorPtr<Real> vec_dot,
                   Real h = 1e-40,
                   bool sum_ghosts=false);

}

#endif