#include "euler_model.h"
#include <iostream>

#include "physics/common/slope_limiters.h"
#include "physics/common/vec_field.h"
#include "physics/euler/euler_flux.h"
#include "physics/euler/typedefs.h"
#include "linear_system/sparsity_pattern_mesh.h"
#include "interface_term.h"
#include "interface_term_jac.h"

namespace structured_fv {
namespace euler {

void EulerModel::evaluateRhs(DiscVectorPtr<Real> q, Real t, 
                             DiscVectorPtr<Real> residual)
{
  Fields<Real>& fields = m_fields_real;
  common::vecToField(m_disc, q, fields.solution);
  evaluateRhsT(fields, t);
  common::fieldToVec(m_disc, fields.residual, residual);
}

void EulerModel::evaluateJacobian(disc::DiscVectorPtr<Real> q, Real t, disc::DiscVectorPtr<Real> residual,
                                  linear_system::AssemblerBasePtr assembler)
{
  Fields<Real> fields = m_fields_real;
  common::vecToField(m_disc, q, fields.solution);
  setBCValues(fields.solution, t);
  fields.solution->updateGhostValues();
  fields.residual->set(0);
  evaluateInterfaceTermsJac(m_opts, fields, t, m_disc, assembler);
}

void EulerModel::computeJacVecProduct(disc::DiscVectorPtr<Real> q, Real t, disc::DiscVectorPtr<Real> v,
                                      disc::DiscVectorPtr<Real> h)
{
  Fields<Dual1> fields = m_fields_dual;
  common::vecToFieldDot(m_disc, q, v, fields.solution);
  evaluateRhsT(fields, t);
  common::fieldToVecDot(m_disc, fields.residual, h);
}

template <typename T>
void EulerModel::evaluateRhsT(Fields<T>& fields, Real t)
{
  setBCValues(fields.solution, t);
  fields.solution->updateGhostValues();
  checkPositivity(fields.solution);
  //std::cout << "max wave speed = " << computeMaxWaveSpeed(fields.solution) << std::endl;
  fields.residual->set(0);

  evaluateInterfaceTerms(m_opts, fields, t, m_disc);

  evaluateSourceTerm(t, fields);
}

std::shared_ptr<linear_system::SparsityPattern> EulerModel::getSparsityPattern() const
{
  int stencilsize = m_opts.limiter == common::SlopeLimiter::FirstOrder ? 1 : 2;
  return std::make_shared<linear_system::SparsityPatternDisc>(m_disc, stencilsize);
}


Fxyt& EulerModel::getBCFunction(UInt block_id)
{
  assert(in(m_disc->getGhostBCBlocksIds(), block_id));
  return m_bc_functions[block_id - m_disc->getNumRegularBlocks()];
}

template <typename T>
void EulerModel::setBCValues(ElementFieldPtr<T> solution, Real t)
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
        FixedVec<Real, 2> x = disc::computeCellCentroid(vert_coords, i, j);
        auto vals = func(x[0], x[1], t);
        for (UInt k=0; k < DofsPerCell; ++k)
          sol(i, j, k) = vals[k];
      }
  }
}


template <typename T>
void EulerModel::evaluateSourceTerm(Real t, Fields<T>& fields)
{
  for (UInt block_id : m_disc->getRegularBlocksIds())
  {
    const StructuredBlock& block = m_disc->getBlock(block_id);
    const auto& vert_coords = m_disc->getCoordField()->getData(block_id);
    auto& res = fields.residual->getData(block_id);

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

template <typename T>
void EulerModel::checkPositivity(const ElementFieldPtr<T>& solution)
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

template <typename T>
T EulerModel::computeMaxWaveSpeed(ElementFieldPtr<T> solution)
{
  T max_wave_speed = 0.0;
  for (UInt block_id : m_disc->getRegularBlocksIds())
  {
    const StructuredBlock& block = m_disc->getBlock(block_id);
    const auto& sol = solution->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
      {
        Vec4<T> q = getValues(sol, i, j);
        T Umag = std::sqrt(q[1]*q[1] + q[2]*q[2])/q[0];
        T a = compute_sos(q);
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