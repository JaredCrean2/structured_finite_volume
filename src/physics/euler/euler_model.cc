#include "euler_model.h"

namespace structured_fv {
namespace euler {

void EulerModel::evaluateRhs(DiscVectorPtr<Real> q, Real t, 
                             DiscVectorPtr<Real> residual)
{
  vecToField(m_disc, q, m_solution);
  setBCValues(m_solution, t);
  m_solution->updateGhostValues();
  checkPositivity(m_solution);
  m_residual->set(0);

  //HLLEFlux flux;
  LaxFriedrichFlux flux;
  evaluateInterfaceTerms(m_solution, t, flux, XDirTag(), m_residual);
  evaluateInterfaceTerms(m_solution, t, flux, YDirTag(), m_residual);

  evaluateSourceTerm(t, m_residual);

  fieldToVec(m_disc, m_residual, residual);
}


Fxyt& EulerModel::getBCFunction(UInt block_id)
{
  assert(in(m_disc->getGhostBCBlocksIds(), block_id));
  return m_bc_functions[block_id - m_disc->getNumRegularBlocks()];
}

void EulerModel::setBCValues(ElementFieldPtr<Real> solution, Real t)
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
        auto vals = func(x[0], x[1], t);
        for (UInt k=0; k < DofsPerCell; ++k)
          sol(i, j, k) = vals[k];
      }
  }
}

template <typename Flux, typename Tag>
void EulerModel::evaluateInterfaceTerms(const ElementFieldPtr<Real>& solution, Real t,
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
        Vec4<Real> qL        = getValues(sol, face_id.cell_i_left, face_id.cell_j_left);
        Vec4<Real> qR        = getValues(sol, face_id.cell_i_right, face_id.cell_j_right);
        Vec2<Real> normal{normals(i, j, 0), normals(i, j, 1)};
        Vec4<Real> flux      = flux_func(qL, qR, normal);
        Real inv_volL = cell_inv_volume(face_id.cell_i_left, face_id.cell_j_left, 0);
        Real inv_volR = cell_inv_volume(face_id.cell_i_right, face_id.cell_j_right, 0);

        for (UInt k=0; k < DofsPerCell; ++k)
        {
          res(face_id.cell_i_left, face_id.cell_j_left, k)   -=  inv_volL * flux[k];
          res(face_id.cell_i_right, face_id.cell_j_right, k) +=  inv_volR * flux[k];
        }
      }
  }
}

void EulerModel::evaluateSourceTerm(Real t, ElementFieldPtr<Real> residual)
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
        auto vals = m_source_func(x[0], x[1], t);
        for (UInt k=0; k < DofsPerCell; ++k)
          res(i, j, k) += vals[k];
      }
  }
}

void EulerModel::checkPositivity(const ElementFieldPtr<Real>& solution)
{
  std::stringstream ss;
  for (UInt block_id : m_disc->getRegularBlocksIds())
  {
    const StructuredBlock& block = m_disc->getBlock(block_id);
    const auto& sol = solution->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
      {
        auto q = getValues(sol, i, j);
        if (q[0] < 0)
        {
          ss << "negative density found in block " << block_id << ", cell "
             << i << ", " << j << " with q = " << q;
          throw std::runtime_error(ss.str());
        }

        if (compute_pressure(q) < 0)
        {
          ss << "negative pressure found in block " << block_id << ", cell "
             << i << ", " << j << " with q = " << q
             << ", p = " << compute_pressure(q);
          throw std::runtime_error(ss.str());          
        }

        if (compute_temperature(q) < 0)
        {
          ss << "negative temperature found in block " << block_id << ", cell "
             << i << ", " << j << " with q = " << q
             << ", T = " << compute_temperature(q) << std::endl;
          throw std::runtime_error(ss.str());          
        }

      }

  }

}

}
}