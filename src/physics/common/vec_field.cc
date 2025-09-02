#include "vec_field.h"

namespace structured_fv {
namespace common {

void vecToFieldDot(disc::StructuredDiscPtr disc, 
                   const disc::DiscVectorPtr<Real> vec, 
                   const disc::DiscVectorPtr<Real> vec_dot,
                   disc::ElementFieldPtr<Dual<Real, 1>> field)
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
        {
          GlobalDof dof = dof_nums(i, j, k);
          field_data(i, j, k) = Dual<Real, 1>((*vec)(dof), {(*vec_dot)(dof)});
        }
  }

  field->updateGhostValues();  
}

void fieldToVecDot(disc::StructuredDiscPtr disc,
                   disc::ElementFieldPtr<Dual<Real, 1>> field,
                   disc::DiscVectorPtr<Real> vec_dot,
                   bool sum_ghosts)
{
  if (sum_ghosts)
    field->reduceGhostValuesToOwner(std::plus());
  
  const UInt vals_per_cell = field->getNumValsPerElement();
  for (UInt block_id=0; block_id < disc->getNumRegularBlocks(); ++block_id)
  {
    const disc::StructuredBlock& block = disc->getBlock(block_id);
    auto& field_data = field->getData(block_id);
    auto& dof_nums = disc->getDofNumbering()->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        for (UInt k=0; k < vals_per_cell; ++k)
          (*vec_dot)(dof_nums(i, j, k)) = field_data(i, j, k).get_deriv(0);
  }  
}

}
}