#ifndef STRUCTURED_FV_PHYSICS_COMMON_VECFIELD_H
#define STRUCTURED_FV_PHYSICS_COMMON_VECFIELD_H

#include "disc/disc_vector.h"
#include "disc/elem_field.h"
#include <utils/dual_number.h>

namespace structured_fv {
namespace common {

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
                   disc::ElementFieldPtr<Dual<Real, 1>> field);

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
                   disc::ElementFieldPtr<Dual<Real, 1>> field,
                   disc::DiscVectorPtr<Real> vec_dot,
                   bool sum_ghosts=false);

}
}

#endif