#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_PHYSICS_MODEL_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_PHYSICS_MODEL_H

#include "disc/disc_block.h"
#include "disc/disc_vector.h"
#include "disc/discretization.h"
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

    virtual Real computeRhsNorm(disc::DiscVectorPtr<Real> residual);

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

}

#endif