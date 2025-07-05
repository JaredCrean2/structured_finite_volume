#include "physics_model.h"

namespace structured_fv {

Real PhysicsModel::computeRhsNorm(disc::DiscVectorPtr<Real> residual) 
{
  // compute integral of the residual over the domain
  Real residual_norm = 0.0;
  const disc::StructuredDiscPtr& disc = getDisc();
  UInt ncomp = disc->getNumDofsPerNode();

  for (UInt block_id : disc->getRegularBlocksIds())
  {
    const disc::StructuredBlock& block = disc->getBlock(block_id);
    const auto& inv_cell_volumes = disc->getInvCellVolumeField()->getData(block_id);
    const auto& dof_nums = disc->getDofNumbering()->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
      {
        for (UInt k=0; k < ncomp; ++k)
        {
          Real volume = 1.0/inv_cell_volumes(i, j, 0);
          Real integral_residual = volume * (*residual)(dof_nums(i, j, k));
          residual_norm += integral_residual * integral_residual;
        }
      }
  }

  return std::sqrt(residual_norm);
}

void vecToFieldDot(disc::StructuredDiscPtr disc, 
                   const disc::DiscVectorPtr<Real> vec, 
                   const disc::DiscVectorPtr<Real> vec_dot,
                   disc::ElementFieldPtr<Complex> field,
                   Real h)
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
          field_data(i, j, k) = Complex((*vec)(dof), h*(*vec_dot)(dof));
        }
  }

  field->updateGhostValues();  
}

void fieldToVecDot(disc::StructuredDiscPtr disc,
                   disc::ElementFieldPtr<Complex> field,
                   disc::DiscVectorPtr<Real> vec_dot,
                   Real h,
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
          (*vec_dot)(dof_nums(i, j, k)) = field_data(i, j, k).imag()/h; 
  }  
}

}