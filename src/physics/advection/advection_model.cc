#include "advection_model.h"
#include "disc/discretization.h"
#include "utils/face_iter_per_direction.h"
#include "physics/common/slope_limiters.h"
#include "utils/math.h"
#include "disc/elem_field.h"
#include "utils/face_iter_per_direction.h"
#include "disc/face_field.h"
#include "utils/neighbor_direction.h"
#include "linear_system/sparsity_pattern_mesh.h"

#include "interface_term.h"
#include "interface_term_jac.h"
#include "source_term.h"

#include <iostream>

namespace structured_fv {
namespace advection {


AdvectionModel::AdvectionModel(const AdvectionOpts& opts, StructuredDiscPtr disc,
                               const std::vector<Fxyt>& bc_functions,
                               Fxyt source_func) :
  m_opts(opts),
  m_disc(disc),
  m_fields_real(*disc),
  m_fields_complex(*disc),
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

void AdvectionModel::evaluateRhs(DiscVectorPtr<Real> q, Real t, 
                                 DiscVectorPtr<Real> residual)
{
  Fields<Real> fields = m_fields_real;
  vecToField(m_disc, q, fields.solution);
  evaluateRhsT(fields, t);
  fieldToVec(m_disc, fields.residual, residual);
}

// evaluate the right hand side of the equation:
// dq/dt = R(q, t)
template <typename T>
void AdvectionModel::evaluateRhsT(Fields<T> fields, Real t)
{
  setBCValues(fields.solution, t);
  fields.solution->updateGhostValues();
  fields.residual->set(0);

  evaluateInterfaceTerm(m_opts, fields, t, m_disc);
  evaluateSourceTerm(t, m_disc, fields.residual, m_source_func);
}

void AdvectionModel::evaluateJacobian(disc::DiscVectorPtr<Real> q, Real t, disc::DiscVectorPtr<Real> residual, linear_system::AssemblerBasePtr assembler)
{
  Fields<Real> fields = m_fields_real;
  vecToField(m_disc, q, fields.solution);
  setBCValues(fields.solution, t);
  fields.solution->updateGhostValues();
  fields.residual->set(0);
  evaluateInterfaceTermJac(m_opts, fields, t, m_disc, assembler);
}

void AdvectionModel::computeJacVecProduct(disc::DiscVectorPtr<Real> q, Real t, disc::DiscVectorPtr<Real> v, disc::DiscVectorPtr<Real> h)
{
  Fields<Complex> fields = m_fields_complex;
  vecToFieldDot(m_disc, q, v, fields.solution);
  evaluateRhsT(fields, t);
  fieldToVecDot(m_disc, fields.residual, h);
}

std::shared_ptr<linear_system::SparsityPattern> AdvectionModel::getSparsityPattern() const
{
  // currently we write to the matrix entries for a stencil size of 2, even when
  // using the first order scheme.
  return std::make_shared<linear_system::SparsityPatternDisc>(m_disc, 2);
}

Fxyt& AdvectionModel::getBCFunction(UInt block_id)
{
  assert(in(m_disc->getGhostBCBlocksIds(), block_id));
  return m_bc_functions[block_id - m_disc->getNumRegularBlocks()];
}


template <typename T>
void AdvectionModel::setBCValues(ElementFieldPtr<T> solution, Real t)
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
        sol(i, j, 0) = func(x[0], x[1], t);
      }
  }
}

/*
template <typename Flux, typename SlopeLimiter, typename Assembler, typename Tag>
void AdvectionModel::evaluateInterfaceTermsJac(const ElementFieldPtr<Real>& solution, Real t,
                                            const Flux& flux_func, const SlopeLimiter& limiter, 
                                            Tag dir_tag, ElementFieldPtr<Real> residual, Assembler& assembler)
{
  using Indices = linear_system::Indices;

  constexpr double epsilon = 1e-15;
  NeighborDirection dir = toNeighborDirection(dir_tag);
  FixedVec<Indices, 1> row_indices;
  FixedVec<Indices, 4> col_indices;
  Matrix<Real, 1, 4> jac;

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

        Real deltaRL   = qR - qL + epsilon;
        Real rL        = (qL - qLm1)/deltaRL;
        Real rL_dotLm1 = -1/deltaRL;
        Real rL_dotL   = 1/deltaRL + (qL - qLm1)/(deltaRL*deltaRL);
        Real rL_dotR   = -(qL - qLm1)/deltaRL;

        Real deltaRp1R = qRp1 - qR + epsilon;
        Real rR        = (qR - qL)/deltaRp1R;
        Real rR_dotL   = -1/deltaRp1R;
        Real rR_dotR   = 1/deltaRp1R + (qR - qL)/(deltaRp1R*deltaRp1R);
        Real rR_dotRp1 = -(qR - qL)/(deltaRp1R*deltaRp1R);


        Real slopeL = (qR - qLm1)/2;
        Real slopeL_dotR = 1.0/2;
        Real slopeL_dotLm1 = -1.0/2;

        Real slopeR = (qRp1 - qL)/2;
        Real slopeR_dotRp1 = 1.0/2;
        Real slopeR_dotL = -1.0/2;

        auto [phiL, phiL_dot] = limiter(rL, Real(1));
        Real phiL_dotLm1 = phiL_dot * rL_dotLm1;
        Real phiL_dotL   = phiL_dot * rL_dotL;
        Real phiL_dotR   = phiL_dot * rL_dotR;

        auto [phiR, phiR_dot] = limiter(rR, Real(1));
        Real phiR_dotL   = phiR_dot * rR_dotL;
        Real phiR_dotR   = phiR_dot * rR_dotR;
        Real phiR_dotRp1 = phiR_dot * rR_dotRp1;

        Real qLhalf = qL + 0.5*phiL*slopeL;
        Real qLhalf_dotLm1 = 0.5*(phiL_dotL*slopeL + phiL*slopeL_dotLm1);
        Real qLhalf_dotL   = 1 + 0.5*(phiL_dotL*slopeL);
        Real qLhalf_dotR   = 0.5*(phiL_dotL*slopeL + phiL*slopeL_dotR);

        Real qRhalf      = qR - 0.5*phiR*slopeR;
        Real qRhalf_dotL = -0.5*(phiR_dotL*slopeR + phiR*slopeR_dotL);
        Real qRhalf_dotR = 1 - 0.5*(phiR_dotR*slopeR);
        Real qRhalf_dotRp1 = -0.5*(phiR_dotRp1*slopeR + phiR*slopeR_dotRp1);

        Vec2<Real> normal{normals(i, j, 0), normals(i, j, 1)};
        Real flux_dotqLhalf = 0, flux_dotqRhalf = 0;
        Real flux      = flux_func(qLhalf, qRhalf, normal, flux_dotqLhalf, flux_dotqRhalf);

        Real flux_dotLm1 = flux_dotqLhalf*qLhalf_dotLm1;
        Real flux_dotL   = flux_dotqLhalf*qLhalf_dotL + flux_dotqRhalf*qRhalf_dotL;
        Real flux_dotR   = flux_dotqLhalf*qLhalf_dotR + flux_dotqRhalf*qRhalf_dotR;
        Real flux_dotRp1 = flux_dotqRhalf*qRhalf_dotRp1;


        Real inv_volL = cell_inv_volume(face_id.cell_i_left, face_id.cell_j_left, 0);
        Real inv_volR = cell_inv_volume(face_id.cell_i_right, face_id.cell_j_right, 0);
        res(face_id.cell_i_left, face_id.cell_j_left, 0)   -= inv_volL * flux;
        res(face_id.cell_i_right, face_id.cell_j_right, 0) += inv_volR * flux;

        jac(0, 0) = inv_volL * flux_dotLm1;
        jac(0, 1) = inv_volL * flux_dotL;
        jac(0, 2) = inv_volL * flux_dotR;
        jac(0, 3) = inv_volL * flux_dotRp1;

        row_indices[0] = {face_id.cell_i_left, face_id.cell_j_left};
        for (UInt k=0; k < 4; ++k)
        {
          col_indices[0] = increment(dir_tag, cell_im1_left, cell_jm1_left, k);
        }
        assembler.assembleValues(row_indices, col_indices, jac);


        jac(0, 0) = inv_volR * flux_dotLm1;
        jac(0, 1) = inv_volR * flux_dotL;
        jac(0, 2) = inv_volR * flux_dotR;
        jac(0, 3) = inv_volR * flux_dotRp1;
        row_indices[0] = {face_id.cell_i_right, face_id.cell_j_right};
        assembler.assembleValues(row_indices, col_indices, jac);

        // Note: this works for dirichlet BCs, but for any kind of characteristic BC
        //       where the ghost cell value is a function of the interior value, that
        //       contribution needs to be included
      }     
  }
}

*/

}
}