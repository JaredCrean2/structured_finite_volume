#include "euler_model.h"
#include "linear_system/assembler_base.h"
#include "physics/euler/numerical_flux_base.h"
#include "physics/euler/reconstruction_enum.h"
#include "physics/euler/typedefs.h"
#include "physics/common/slope_limiters.h"
#include "physics/euler/reconstruction.h"
#include "roe_flux.h"
#include "roe_hh_flux.h"
#include "hlle_flux.h"
#include "lax_friedrich_flux.h"
#include "hllc_flux.h"
#include "disc/face_field.h"
#include "linear_system/assembler.h"
#include "linear_system/large_matrix_dense.h"
#include "linear_system/large_matrix_petsc.h"

namespace structured_fv {
namespace euler {

using linear_system::AssemblerBasePtr;

template <typename T, typename Reconstruction, typename Flux, typename Tag, typename Assembler>
void evaluateInterfaceTermsJacImpl(Fields<T>& fields, Real t,
                                  const StructuredDiscPtr disc, 
                                  const Reconstruction& recon, const Flux& flux_func,
                                  Tag dir_tag, Assembler& assembler)
{
  static_assert(IsReconstruction<Reconstruction>);
  static_assert(IsNumericalFlux<Flux>);
  static_assert(linear_system::IsAssembler<Assembler>);

  using Indices = linear_system::Indices;
  NeighborDirection dir = toNeighborDirection(dir_tag);
  constexpr bool is_first_order = std::is_same_v<typename Reconstruction::Limiter, common::SlopeLimiterFirstOrder>;

  for (UInt block_id : disc->getRegularBlocksIds())
  {
    const StructuredBlock& block = disc->getBlock(block_id);
    const auto& sol              = fields.solution->getData(block_id);
    const auto& cell_inv_volume  = disc->getInvCellVolumeField()->getData(block_id);
    const auto& normals          = disc->getNormalField()->getData(block_id, dir);
    auto& res                    = fields.residual->getData(block_id);
    assembler.setBlock(block_id);

    Range2D owned_cells = block.getOwnedCells();
    FaceRangePerDirection faces = block.getOwnedFaces();
    for (UInt i : faces.getXRange(dir_tag))
      for (UInt j : faces.getYRange(dir_tag))
      {
        FaceId face_id = faces.getFaceId(dir_tag, i, j);
        Vec4<T> flux;
        if constexpr (is_first_order)
        {
          FaceId face_id = faces.getFaceId(dir_tag, i, j);
          Vec4<T> qL        = getValues(sol, face_id.cell_i_left, face_id.cell_j_left);
          Vec4<T> qR        = getValues(sol, face_id.cell_i_right, face_id.cell_j_right);
          Vec2<Real> normal{normals(i, j, 0), normals(i, j, 1)};
          Matrix<T, 4> flux_dotL, flux_dotR;
          flux = flux_func(qL, qR, normal, flux_dotL, flux_dotR);

          Real inv_volL = cell_inv_volume(face_id.cell_i_left, face_id.cell_j_left, 0);
          Real inv_volR = cell_inv_volume(face_id.cell_i_right, face_id.cell_j_right, 0);
          Matrix<T, 2*DofsPerCell, 2*DofsPerCell> jac;

          Vec2<Indices> row_idx = {{{face_id.cell_i_left, face_id.cell_j_left},
                                    {face_id.cell_i_right, face_id.cell_j_right}}};
          Vec2<Indices>& col_idx = row_idx;

          //TODO: this is not a great solution in  parallel because it requires the matrix
          //      to do parallel communication during assembly for values that are all zeros
          UInt owned_maskL = in(owned_cells, face_id.cell_i_left, face_id.cell_j_left) ? 1 : 0;
          UInt owned_maskR = in(owned_cells, face_id.cell_i_right, face_id.cell_j_right) ? 1 : 0;
          for (UInt k=0; k < DofsPerCell; ++k)
          {
            res(face_id.cell_i_left, face_id.cell_j_left, k)   -=  inv_volL * flux[k];
            res(face_id.cell_i_right, face_id.cell_j_right, k) +=  inv_volR * flux[k];

            for (UInt j=0; j < DofsPerCell; ++j)
            {
              jac(k, j)                             = owned_maskL * -inv_volL * flux_dotL(k, j);
              jac(k, j + DofsPerCell)               = owned_maskL * -inv_volL * flux_dotR(k, j);
              jac(k + DofsPerCell, j)               = owned_maskR *  inv_volR * flux_dotL(k, j);
              jac(k + DofsPerCell, j + DofsPerCell) = owned_maskR *  inv_volR * flux_dotR(k, j);
            }
          }

          // Note: this works because the dofs array used by the assembler has a negative value
          //       for all non-owned entries, meaning Petsc will ignore them.  If that
          //       wasn't true, we would have to split this into two jacobians, one for the
          //       left and one for the right element.
          assembler.assembleValues(row_idx, col_idx, jac);

        } else        
        {
          const auto [cell_im1, cell_jm1] = increment(dir_tag, face_id.cell_i_left, face_id.cell_j_left, -1);
          const auto [cell_ip1, cell_jp1] = increment(dir_tag, face_id.cell_i_right, face_id.cell_j_right, 1);
          //TODO: this uses a lot of registers.  Make a non-owning vector class that refers to global memory?
          Vec4<T> qLm1      = getValues(sol, cell_im1, cell_jm1);
          Vec4<T> qL        = getValues(sol, face_id.cell_i_left, face_id.cell_j_left);
          Vec4<T> qR        = getValues(sol, face_id.cell_i_right, face_id.cell_j_right);
          Vec4<T> qRp1      = getValues(sol, cell_ip1, cell_jp1);
          Vec2<Real> normal{normals(i, j, 0), normals(i, j, 1)};

          Matrix<T, DofsPerCell> qLhalf_dotLm1, qLhalf_dotL, qLhalf_dotR;
          Matrix<T, DofsPerCell> qRhalf_dotL, qRhalf_dotR, qRhalf_dotRp1;

          Vec4<T> qLhalf = recon(qLm1, qLhalf_dotLm1, qL, qLhalf_dotL, qR,   qLhalf_dotR, 1);
          Vec4<T> qRhalf = recon(qL,   qRhalf_dotL,   qR, qRhalf_dotR, qRp1, qRhalf_dotRp1, -1);

          Matrix<T, DofsPerCell> flux_dotqLhalf, flux_dotqRhalf;
          flux      = flux_func(qLhalf, qRhalf, normal, flux_dotqLhalf, flux_dotqRhalf);

          Real inv_volL = cell_inv_volume(face_id.cell_i_left, face_id.cell_j_left, 0);
          Real inv_volR = cell_inv_volume(face_id.cell_i_right, face_id.cell_j_right, 0);

          Matrix<T, 2*DofsPerCell, 4*DofsPerCell> jac;
          Vec2<Indices> row_idx = {{{face_id.cell_i_left, face_id.cell_j_left},
                                    {face_id.cell_i_right, face_id.cell_j_right}}};

          Vec4<Indices> col_idx = {{{cell_im1, cell_jm1}, 
                                    {face_id.cell_i_left, face_id.cell_j_left},
                                    {face_id.cell_i_right, face_id.cell_j_right},
                                    {cell_ip1, cell_jp1}}};

          UInt owned_maskL = in(owned_cells, face_id.cell_i_left, face_id.cell_j_left) ? 1 : 0;
          UInt owned_maskR = in(owned_cells, face_id.cell_i_right, face_id.cell_j_right) ? 1 : 0;    

          for (UInt k=0; k < DofsPerCell; ++k)
          {
            res(face_id.cell_i_left, face_id.cell_j_left, k)   -=  inv_volL * flux[k];
            res(face_id.cell_i_right, face_id.cell_j_right, k) +=  inv_volR * flux[k];

            for (UInt p=0; p < DofsPerCell; ++p)
            {
              T dR_dqLm1 = 0; //, dRR_dqLm1 = 0;
              T dR_dqL   = 0; //, dRR_dqL   = 0;
              T dR_dqR   = 0; //, dRR_dqR   = 0;
              T dR_dqRp1 = 0; //, dRR_dqRp1 = 0;
              for (UInt q=0; q < DofsPerCell; ++q)
              {
                dR_dqLm1 += flux_dotqLhalf(k, q) * qLhalf_dotLm1(q, p);
                dR_dqL   += flux_dotqLhalf(k, q) * qLhalf_dotL(q, p) + flux_dotqRhalf(k, q) * qRhalf_dotL(q, p);
                dR_dqR   += flux_dotqLhalf(k, q) * qLhalf_dotR(q, p) + flux_dotqRhalf(k, q) * qRhalf_dotR(q, p);
                dR_dqRp1 +=                                            flux_dotqRhalf(k, q) * qRhalf_dotRp1(q, p);
              }

              //TODO: because of the striding in p, this can't be a simd store
              jac(k, p)                 = owned_maskL * -inv_volL * dR_dqLm1;
              jac(k, p + DofsPerCell)   = owned_maskL * -inv_volL * dR_dqL;
              jac(k, p + 2*DofsPerCell) = owned_maskL * -inv_volL * dR_dqR;
              jac(k, p + 3*DofsPerCell) = owned_maskL * -inv_volL * dR_dqRp1;

              jac(k + DofsPerCell, p)                 = owned_maskR * inv_volR * dR_dqLm1;
              jac(k + DofsPerCell, p + DofsPerCell)   = owned_maskR * inv_volR * dR_dqL;
              jac(k + DofsPerCell, p + 2*DofsPerCell) = owned_maskR * inv_volR * dR_dqR;
              jac(k + DofsPerCell, p + 3*DofsPerCell) = owned_maskR * inv_volR * dR_dqRp1;
            }
          }

          assembler.assembleValues(row_idx, col_idx, jac);
        }
      }
  }
}

template <typename T, typename Reconstruction, typename Flux, typename Assembler>
void evaluateInterfaceTermsJacAddDirection(Fields<T>& fields, Real t, const StructuredDiscPtr disc, 
                                           const Reconstruction& recon, const Flux& flux_func,
                                           Assembler& assembler)
{
  static_assert(IsReconstruction<Reconstruction>);
  static_assert(IsNumericalFlux<Flux>);
  static_assert(linear_system::IsAssembler<Assembler>);

  evaluateInterfaceTermsJacImpl(fields, t, disc, recon, flux_func, XDirTag(), assembler);
  evaluateInterfaceTermsJacImpl(fields, t, disc, recon, flux_func, YDirTag(), assembler);
}

template <typename T, typename Reconstruction, typename Flux>
void evaluateInterfaceTermsJacAddAssembler(Fields<T>& fields, Real t, const StructuredDiscPtr disc, 
                                           const Reconstruction& recon, const Flux& flux_func,
                                           AssemblerBasePtr assembler)
{
  static_assert(IsReconstruction<Reconstruction>);
  static_assert(IsNumericalFlux<Flux>);

  using AssemblerDense = linear_system::Assembler<linear_system::LargeMatrixDense>;
  using AssemblerPetsc = linear_system::Assembler<linear_system::LargeMatrixPetsc>;

  auto assembler_dense = std::dynamic_pointer_cast<AssemblerDense>(assembler);
  auto assembler_petsc = std::dynamic_pointer_cast<AssemblerPetsc>(assembler);

  if (assembler_dense)
  {
    evaluateInterfaceTermsJacAddDirection(fields, t, disc, recon, flux_func, *assembler_dense);
  } else if (assembler_petsc)
  {
    evaluateInterfaceTermsJacAddDirection(fields, t, disc, recon, flux_func, *assembler_petsc);
  } else
  {
    throw std::runtime_error("unhandled Assembler type");
  }
}

template <typename T, typename Reconstruction>
void evaluateInterfaceTermsJacFluxFunction(
       const EulerOpts& opts, Fields<T>& fields, Real t,
       const StructuredDiscPtr disc, const Reconstruction& recon,
       AssemblerBasePtr assembler)
{
  static_assert(IsReconstruction<Reconstruction>);

  if (opts.flux == FluxFunction::Roe)
  {
    RoeFlux flux(opts.roe_efix_delta);
    evaluateInterfaceTermsJacAddAssembler(fields, t, disc, recon, flux, assembler);
  } else if (opts.flux == FluxFunction::RoeHH)
  {
    RoeHHFlux flux;
    evaluateInterfaceTermsJacAddAssembler(fields, t, disc, recon, flux, assembler);
  } else if (opts.flux == FluxFunction::HLLE)
  {
    HLLEFlux flux;
    evaluateInterfaceTermsJacAddAssembler(fields, t, disc, recon, flux, assembler);
  } else if (opts.flux == FluxFunction::LLF)
  {
    LaxFriedrichFlux flux;
    evaluateInterfaceTermsJacAddAssembler(fields, t, disc, recon, flux, assembler);
  } else if (opts.flux == FluxFunction::HLLC)
  {
    HLLCFlux flux;
    evaluateInterfaceTermsJacAddAssembler(fields, t, disc, recon, flux, assembler);
  } else
  {
    throw std::runtime_error("unhandled FluxFunction enum");
  }  
}


template <typename T, typename SlopeLimiter>
void evaluateInterfaceTermsJacAddReconstruction(
       const EulerOpts& opts, Fields<T>& fields, Real t,
       const StructuredDiscPtr disc, const SlopeLimiter& limiter,
       AssemblerBasePtr assembler)
{
  static_assert(common::IsSlopeLimiter<SlopeLimiter>);

  if (opts.recon == Reconstruction::Conservative) {
    ReconstructionConservative<SlopeLimiter> recon(limiter);  
    evaluateInterfaceTermsJacFluxFunction(opts, fields, t, disc, recon, assembler);
  } else if (opts.recon == Reconstruction::Primitive)
  {
    ReconstructionPrimitive<SlopeLimiter> recon(limiter);
    evaluateInterfaceTermsJacFluxFunction(opts, fields, t, disc, recon, assembler); 
  } else
  {
    throw std::runtime_error("unhandled reconstruction enum " + get_name(opts.recon));
  }
}


template <typename T>
void evaluateInterfaceTermsJacAddLimiter(const EulerOpts& opts, Fields<T>& fields, Real t, const StructuredDiscPtr disc,
                                         AssemblerBasePtr assembler)
{
  if (opts.limiter == common::SlopeLimiter::FirstOrder)
  {
    common::SlopeLimiterFirstOrder limiter;   
    evaluateInterfaceTermsJacAddReconstruction(opts, fields, t, disc, limiter, assembler);
  } else if (opts.limiter == common::SlopeLimiter::MinMod)
  {
    common::SlopeLimiterMinMod limiter;
    evaluateInterfaceTermsJacAddReconstruction(opts, fields, t, disc, limiter, assembler);
  } else if (opts.limiter == common::SlopeLimiter::SuperBee)
  {
    common::SlopeLimiterSuperBee limiter;
    evaluateInterfaceTermsJacAddReconstruction(opts, fields, t, disc, limiter, assembler);
  } else if (opts.limiter == common::SlopeLimiter::VanAlba)
  {
    common::SlopeLimiterVanAlba limiter;
    evaluateInterfaceTermsJacAddReconstruction(opts, fields, t, disc, limiter, assembler);
  } else if (opts.limiter == common::SlopeLimiter::VanLeer)
  {
    common::SlopeLimiterVanLeer limiter;
    evaluateInterfaceTermsJacAddReconstruction(opts, fields, t, disc, limiter, assembler); 
  } else
  {
    throw std::runtime_error("unhandled SlopeLimiter enum " + get_name(opts.limiter));
  }
}


void evaluateInterfaceTermsJac(const EulerOpts& opts, Fields<Real>& fields, Real t, const StructuredDiscPtr disc,
                               AssemblerBasePtr assembler)
{
  // Unlike a factory, which takes concrete types and returns a pointer to an abstract
  // type, this code builds up the template arguments one by one until everything
  // has a concrete type.
  // Each function in the chain uses an enum from EulerOpts and figures out the
  // concrete type (similar to a factory), but then calls the next function in the chain.
  evaluateInterfaceTermsJacAddLimiter(opts, fields, t, disc, assembler);
}


}
}