#include "euler_model.h"
#include <iostream>

#include "physics/euler/euler_flux.h"
#include "physics/euler/typedefs.h"
#include "interface_term.h"

namespace structured_fv {
namespace euler {


void EulerModel::evaluateRhs(DiscVectorPtr<Real> q, Real t, 
                             DiscVectorPtr<Real> residual)
{
  vecToField(m_disc, q, m_solution);
  setBCValues(m_solution, t);
  m_solution->updateGhostValues();
  checkPositivity(m_solution);
  std::cout << "max wave speed = " << computeMaxWaveSpeed(m_solution) << std::endl;
  m_residual->set(0);

  evaluateInterfaceTermsEntryPoint(m_opts, m_solution, t, m_disc, m_residual);

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


Real EulerModel::computeMaxWaveSpeed(ElementFieldPtr<Real> solution)
{
  Real max_wave_speed = 0.0;
  for (UInt block_id : m_disc->getRegularBlocksIds())
  {
    const StructuredBlock& block = m_disc->getBlock(block_id);
    const auto& sol = solution->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
      {
        Vec4<Real> q = getValues(sol, i, j);
        Real Umag = std::sqrt(q[1]*q[1] + q[2]*q[2])/q[0];
        Real a = compute_sos(q);
        /*
        if (Umag + a > max_wave_speed)
        {
          std::cout << "found new max wave speed" << std::endl;
          std::cout << "q = " << q << std::endl;
          std::cout << "U = " << Umag << ", a = " << a << std::endl;
          std::cout << "p = " << compute_pressure(q) << std::endl;
        }
        */
        max_wave_speed = std::max(max_wave_speed, Umag + a);
      }
  }

  return max_wave_speed;  
}


}
}