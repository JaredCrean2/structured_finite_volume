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

  return sqrt(residual_norm);
}

}