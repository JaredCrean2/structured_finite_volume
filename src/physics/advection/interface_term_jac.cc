

#include "interface_term_jac.h"
#include "interface_term.h"
#include "linear_system/assembler.h"
#include "linear_system/large_matrix_dense.h"
#include "linear_system/large_matrix_petsc.h"
#include "utils/matrix.h"

namespace structured_fv {
namespace advection {

template <typename SlopeLimiter, typename Assembler, typename Tag>
void evaluateInterfaceTermsJacImpl(Fields<Real> fields, Real t, StructuredDiscPtr disc,
                                   const FluxFunctionUpwind& flux_func, const SlopeLimiter& limiter, 
                                   Tag dir_tag, Assembler& assembler)
{
  using Indices = linear_system::Indices;

  constexpr double epsilon = 1e-15;
  NeighborDirection dir = toNeighborDirection(dir_tag);
  FixedVec<Indices, 1> row_indices;
  FixedVec<Indices, 4> col_indices;
  Matrix<Real, 1, 4> jac;

  for (UInt block_id : disc->getRegularBlocksIds())
  {
    const StructuredBlock& block = disc->getBlock(block_id);
    const auto& sol              = fields.solution->getData(block_id);
    const auto& cell_inv_volume  = disc->getInvCellVolumeField()->getData(block_id);
    const auto& normals          = disc->getNormalField()->getData(block_id, dir);
    auto& res                    = fields.residual->getData(block_id);
    
    assembler.setBlock(block_id);
    FaceRangePerDirection faces = block.getOwnedFaces();
    Range2D owned_cells = block.getOwnedCells();
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
        Real rL_dotR   = -(qL - qLm1)/(deltaRL*deltaRL);

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
        Real qLhalf_dotLm1 = 0.5*(phiL_dotLm1*slopeL + phiL*slopeL_dotLm1);
        Real qLhalf_dotL   = 1 + 0.5*(phiL_dotL*slopeL);
        Real qLhalf_dotR   = 0.5*(phiL_dotR*slopeL + phiL*slopeL_dotR);

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

        row_indices[0] = {face_id.cell_i_left, face_id.cell_j_left};
        for (UInt k=0; k < 4; ++k)
        {
          const auto [cell_i, cell_j] = increment(dir_tag, cell_im1_left, cell_jm1_left, k);
          col_indices[k] = {cell_i, cell_j};
        }

        if (in(owned_cells, row_indices[0].i, row_indices[0].j))
        {
          jac(0, 0) = -inv_volL * flux_dotLm1;
          jac(0, 1) = -inv_volL * flux_dotL;
          jac(0, 2) = -inv_volL * flux_dotR;
          jac(0, 3) = -inv_volL * flux_dotRp1;

          assembler.assembleValues(row_indices, col_indices, jac);
        }

        row_indices[0] = {face_id.cell_i_right, face_id.cell_j_right};
        if (in(owned_cells, row_indices[0].i, row_indices[0].j))
        {
          jac(0, 0) = inv_volR * flux_dotLm1;
          jac(0, 1) = inv_volR * flux_dotL;
          jac(0, 2) = inv_volR * flux_dotR;
          jac(0, 3) = inv_volR * flux_dotRp1;
          assembler.assembleValues(row_indices, col_indices, jac);
        }

        // Note: this works for dirichlet BCs, but for any kind of characteristic BC
        //       where the ghost cell value is a function of the interior value, that
        //       contribution needs to be included
      }     
  }
}

template <typename SlopeLimiter, typename Assembler>
void evaluateInterfaceTermsJacAddDirection(Fields<Real> fields, Real t, StructuredDiscPtr disc,
                                            FluxFunctionUpwind flux, SlopeLimiter limiter,
                                           Assembler& assembler)
{
  evaluateInterfaceTermsJacImpl(fields, t, disc, flux, limiter, XDirTag(), assembler);
  evaluateInterfaceTermsJacImpl(fields, t, disc, flux, limiter, YDirTag(), assembler);
}

template <typename SlopeLimiter>
void evaluateInterfaceTermJacAddAssembler(Fields<Real> fields, Real t, StructuredDiscPtr disc,
                                          FluxFunctionUpwind flux, SlopeLimiter limiter,
                                          linear_system::AssemblerBasePtr assembler)
{
  using AssemblerDense = linear_system::Assembler<linear_system::LargeMatrixDense>;
  using AssemblerPetsc = linear_system::Assembler<linear_system::LargeMatrixPetsc>;

  auto assembler_dense = std::dynamic_pointer_cast<AssemblerDense>(assembler);
  auto assembler_petsc = std::dynamic_pointer_cast<AssemblerPetsc>(assembler);

  if (assembler_dense)
  {
    evaluateInterfaceTermsJacAddDirection(fields, t, disc, flux, limiter, *assembler_dense);
  } else if (assembler_petsc)
  {
    evaluateInterfaceTermsJacAddDirection(fields, t, disc, flux, limiter, *assembler_petsc);
  } else
  {
    throw std::runtime_error("unhandled Assembler type");
  }
}


void evaluateInterfaceTermJac(const AdvectionOpts& opts, Fields<Real> fields, Real t, StructuredDiscPtr disc,
                              linear_system::AssemblerBasePtr assembler)
{
  FluxFunctionUpwind flux(opts.adv_velocity);  
  if (opts.limiter == common::SlopeLimiter::FirstOrder)
  {
    common::SlopeLimiterFirstOrder limiter;
    evaluateInterfaceTermJacAddAssembler(fields, t, disc, flux, limiter, assembler);
  } else if (opts.limiter == common::SlopeLimiter::MinMod)
  {
    common::SlopeLimiterMinMod limiter;
    evaluateInterfaceTermJacAddAssembler(fields, t, disc, flux, limiter, assembler);
  } else if (opts.limiter == common::SlopeLimiter::SuperBee)
  {
    common::SlopeLimiterSuperBee limiter;
    evaluateInterfaceTermJacAddAssembler(fields, t, disc, flux, limiter, assembler);
  } else if (opts.limiter == common::SlopeLimiter::VanAlba)
  {
    common::SlopeLimiterVanAlba limiter;
    evaluateInterfaceTermJacAddAssembler(fields, t, disc, flux, limiter, assembler);
  } else if (opts.limiter == common::SlopeLimiter::VanLeer)
  {
    common::SlopeLimiterVanLeer limiter;
    evaluateInterfaceTermJacAddAssembler(fields, t, disc, flux, limiter, assembler);
  } else
  {
    throw std::runtime_error("unsupported SlopeLimiter: " + common::get_name(opts.limiter));
  }
}

}
}