#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_ADVECTION_SOURCE_TERM_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_ADVECTION_SOURCE_TERM_H

#include "utils/project_defs.h"
#include "disc/discretization.h"
#include "disc/elem_field.h"
#include "disc/vert_field.h"


namespace structured_fv {
namespace advection {

using disc::ElementFieldPtr;
using disc::StructuredDiscPtr;

template <typename T, typename SourceFunc>
void evaluateSourceTerm(Real t, StructuredDiscPtr disc, ElementFieldPtr<T> residual, SourceFunc func)
{
  for (UInt block_id : disc->getRegularBlocksIds())
  {
    const disc::StructuredBlock& block = disc->getBlock(block_id);
    const auto& vert_coords = disc->getCoordField()->getData(block_id);
    auto& res = residual->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
      {
        Vec2<Real> x = disc::computeCellCentroid(vert_coords, i, j);
        res(i, j, 0) += func(x[0], x[1], t);
      }
  }
}

}
}

#endif