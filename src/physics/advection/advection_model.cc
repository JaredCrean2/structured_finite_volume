#include "advection_model.h"
#include "disc/discretization.h"
#include "utils/face_iter_per_direction.h"

#include "utils/math.h"
#include "disc/elem_field.h"
#include "utils/face_iter_per_direction.h"
#include "disc/face_field.h"
#include "utils/neighbor_direction.h"
#include <cassert>

#include <iostream>

namespace structured_fv {
namespace advection {

class FluxFunctionUpwind
{
  public:
    FluxFunctionUpwind(const Vec2<Real>& advection_velocity) :
      m_adv_velocity(advection_velocity)
    {}

    Real operator()(Real qL, Real qR, const Vec2<Real>& normal)
    {
      //TODO: do this more efficiently
      Real a_normal = dot(m_adv_velocity, normal);
      return a_normal > 0 ? a_normal * qL : a_normal * qR;
    }

  private:
    Vec2<Real> m_adv_velocity;
};

AdvectionModel::AdvectionModel(const AdvectionOpts& opts, StructuredDiscPtr disc,
                               const std::vector<Fxyt>& bc_functions,
                               Fxyt source_func) :
  m_opts(opts),
  m_disc(disc),
  m_solution(std::make_shared<disc::ElementField<Real>>(*disc, 1)),
  m_residual(std::make_shared<disc::ElementField<Real>>(*disc, 1)),
  m_bc_functions(bc_functions),
  m_source_func(source_func)
{
  if (disc->getNumDofsPerNode() != 1)
    throw std::runtime_error("Advection can only have 1 dof per node");

  if (m_bc_functions.size() != disc->getNumGhostBCBlocks())
    throw std::runtime_error("number of BC functions must be equal to number of BC blocks");
}

// evaluate the right hand side of the equation:
// dq/dt = R(q, t)
void AdvectionModel::evaluateRhs(DiscVectorPtr<Real> q, Real t, 
                                 DiscVectorPtr<Real> residual)
{
  vecToField(m_disc, q, m_solution);
  setBCValues(m_solution, t);
  m_solution->updateGhostValues();
  m_residual->set(0);

  FluxFunctionUpwind flux(m_opts.adv_velocity);
  evaluateInterfaceTerms(m_solution, t, flux, XDirTag(), m_residual);
  evaluateInterfaceTerms(m_solution, t, flux, YDirTag(), m_residual);

  evaluateSourceTerm(t, m_residual);

  fieldToVec(m_disc, m_residual, residual);
}

Fxyt& AdvectionModel::getBCFunction(UInt block_id)
{
  assert(in(m_disc->getGhostBCBlocksIds(), block_id));
  return m_bc_functions[block_id - m_disc->getNumRegularBlocks()];
}

Real AdvectionModel::computeRhsNorm(disc::DiscVectorPtr<Real> residual) 
{
  // compute integral of the residual over the domain
  Real residual_norm = 0.0;

  for (UInt block_id : m_disc->getRegularBlocksIds())
  {
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    const auto& inv_cell_volumes = m_disc->getInvCellVolumeField()->getData(block_id);
    const auto& dof_nums = m_disc->getDofNumbering()->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
      {
        Real volume = inv_cell_volumes(i, j, 0);
        Real integral_residual = volume * (*residual)(dof_nums(i, j, 0));
        residual_norm += integral_residual * integral_residual;
      }
  }

  return std::sqrt(residual_norm);
}


void AdvectionModel::setBCValues(ElementFieldPtr<Real> solution, Real t)
{
  // if we have an analyical solution that extends beyond the domain,
  // putting exact values in the ghost BC cells works.
  // If we only have the value on the boundary, something more sophisticated
  // may be required to compute the ghost cell value s.t. the value at the
  // face comes out correct.
  for (UInt block_id : m_disc->getGhostBCBlocksIds())
  {
    const StructuredBlock& block = m_disc->getBlock(block_id);
    const auto& vert_coords = m_disc->getCoordField()->getData(block_id);
    auto& sol = solution->getData(block_id);
    Fxyt& func = getBCFunction(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
      {
        std::array<Real, 2> x = disc::computeCellCentroid(vert_coords, i, j);
        sol(i, j, 0) = func(x[0], x[1], t);
      }
  }
}


template <typename Flux, typename Tag>
void AdvectionModel::evaluateInterfaceTerms(const ElementFieldPtr<Real>& solution, Real t,
                                            Flux& flux_func, Tag dir_tag, ElementFieldPtr<Real> residual)
{
  NeighborDirection dir = toNeighborDirection(dir_tag);

  for (UInt block_id : m_disc->getRegularBlocksIds())
  {
    const StructuredBlock& block = m_disc->getBlock(block_id);
    const auto& sol              = solution->getData(block_id);
    const auto& cell_inv_volume  = m_disc->getInvCellVolumeField()->getData(block_id);
    const auto& normals          = m_disc->getNormalField()->getData(block_id, dir);
    auto& res                    = residual->getData(block_id);

    FaceRangePerDirection faces = block.getOwnedFaces();
    for (UInt i : faces.getXRange(dir_tag))
      for (UInt j : faces.getYRange(dir_tag))
      {

        FaceId face_id = faces.getFaceId(dir_tag, i, j);
        Real qL        = sol(face_id.cell_i_left, face_id.cell_j_left, 0);
        Real qR        = sol(face_id.cell_i_right, face_id.cell_j_right, 0);
        Vec2<Real> normal{normals(i, j, 0), normals(i, j, 1)};
        Real flux      = flux_func(qL, qR, normal);

        res(face_id.cell_i_left, face_id.cell_j_left, 0)   -= cell_inv_volume(face_id.cell_i_left, face_id.cell_j_left, 0) * flux;
        res(face_id.cell_i_right, face_id.cell_j_right, 0) += cell_inv_volume(face_id.cell_i_right, face_id.cell_j_right, 0) * flux;
      }     
  }
}

void AdvectionModel::evaluateSourceTerm(Real t, ElementFieldPtr<Real> residual)
{
  for (UInt block_id : m_disc->getRegularBlocksIds())
  {
    const StructuredBlock& block = m_disc->getBlock(block_id);
    const auto& vert_coords = m_disc->getCoordField()->getData(block_id);
    auto& res = residual->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
      {
        Vec2<Real> x = disc::computeCellCentroid(vert_coords, i, j);
        res(i, j, 0) += m_source_func(x[0], x[1], t);
      }
  }
}

}
}