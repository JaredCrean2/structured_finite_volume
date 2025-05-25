#include "advection_model.h"
#include "disc/discretization.h"
#include "utils/face_iter_per_direction.h"
#include "physics/common/slope_limiters.h"

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

    constexpr Real operator()(Real qL, Real qR, const Vec2<Real>& normal) const
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

  if (m_disc->getNumGhostCells() < 2)
    throw std::runtime_error("must have at least 2 ghost cells");
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
  if (m_opts.limiter == common::SlopeLimiter::FirstOrder)
  {
    common::SlopeLimiterFirstOrder limiter;
    evaluateInterfaceTerms(m_solution, t, flux, limiter, XDirTag(), m_residual);
    evaluateInterfaceTerms(m_solution, t, flux, limiter, YDirTag(), m_residual);
  } else if (m_opts.limiter == common::SlopeLimiter::MinMod)
  {
    common::SlopeLimiterMinMod limiter;
    evaluateInterfaceTerms(m_solution, t, flux, limiter, XDirTag(), m_residual);
    evaluateInterfaceTerms(m_solution, t, flux, limiter, YDirTag(), m_residual);
  } else if (m_opts.limiter == common::SlopeLimiter::SuperBee)
  {
    common::SlopeLimiterSuperBee limiter;
    evaluateInterfaceTerms(m_solution, t, flux, limiter, XDirTag(), m_residual);
    evaluateInterfaceTerms(m_solution, t, flux, limiter, YDirTag(), m_residual);
  } else if (m_opts.limiter == common::SlopeLimiter::VanAlba)
  {
    common::SlopeLimiterVanAlba limiter;
    evaluateInterfaceTerms(m_solution, t, flux, limiter, XDirTag(), m_residual);
    evaluateInterfaceTerms(m_solution, t, flux, limiter, YDirTag(), m_residual);
  } else if (m_opts.limiter == common::SlopeLimiter::VanLeer)
  {
    common::SlopeLimiterVanLeer limiter;
    evaluateInterfaceTerms(m_solution, t, flux, limiter, XDirTag(), m_residual);
    evaluateInterfaceTerms(m_solution, t, flux, limiter, YDirTag(), m_residual);
  } else
  {
    throw std::runtime_error("unsupported SlopeLimiter: " + common::get_name(m_opts.limiter));
  }

  evaluateSourceTerm(t, m_residual);

  fieldToVec(m_disc, m_residual, residual);
}

Fxyt& AdvectionModel::getBCFunction(UInt block_id)
{
  assert(in(m_disc->getGhostBCBlocksIds(), block_id));
  return m_bc_functions[block_id - m_disc->getNumRegularBlocks()];
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


template <typename Flux, typename SlopeLimiter, typename Tag>
void AdvectionModel::evaluateInterfaceTerms(const ElementFieldPtr<Real>& solution, Real t,
                                            const Flux& flux_func, const SlopeLimiter& limiter, 
                                            Tag dir_tag, ElementFieldPtr<Real> residual)
{
  constexpr double epsilon = 1e-15;
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
        const auto [cell_im1_left, cell_jm1_left] = increment(dir_tag, face_id.cell_i_left, face_id.cell_j_left, -1);
        const auto [cell_ip1_right, cell_jp1_right] = increment(dir_tag, face_id.cell_i_right, face_id.cell_j_right, 1);
        Real qLm1      = sol(cell_im1_left, cell_jm1_left, 0);
        Real qL        = sol(face_id.cell_i_left, face_id.cell_j_left, 0);
        Real qR        = sol(face_id.cell_i_right, face_id.cell_j_right, 0);
        Real qRp1      = sol(cell_ip1_right, cell_jp1_right, 0);

        Real rL = (qL - qLm1)/(qR - qL + epsilon);
        Real rR = (qR - qL)/(qRp1 - qR + epsilon);
        Real slopeL = (qR - qLm1)/2;
        Real slopeR = (qRp1 - qL)/2;

        Real qLhalf = qL + 0.5*limiter(rL)*slopeL;
        Real qRhalf = qR - 0.5*limiter(rR)*slopeR;

        Vec2<Real> normal{normals(i, j, 0), normals(i, j, 1)};
        Real flux      = flux_func(qLhalf, qRhalf, normal);

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