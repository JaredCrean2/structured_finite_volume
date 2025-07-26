#include "euler_model.h"
#include "physics/euler/numerical_flux_base.h"
#include "physics/euler/typedefs.h"
#include "physics/common/slope_limiters.h"
#include "physics/euler/reconstruction.h"
#include "roe_flux.h"
#include "roe_hh_flux.h"
#include "hlle_flux.h"
#include "lax_friedrich_flux.h"
#include "hllc_flux.h"
#include "disc/face_field.h"
#include "interface_term.h"


namespace structured_fv {
namespace euler {


template <typename T, typename Reconstruction, typename Flux, typename Tag>
void evaluateInterfaceTermImpl(Fields<T>& fields, Real t,
                               const StructuredDiscPtr disc, 
                               const Reconstruction& recon, const Flux& flux_func,
                               Tag dir_tag)
{
  NeighborDirection dir = toNeighborDirection(dir_tag);
  constexpr bool is_first_order = std::is_same_v<typename Reconstruction::Limiter, common::SlopeLimiterFirstOrder>;

  for (UInt block_id : disc->getRegularBlocksIds())
  {
    const StructuredBlock& block = disc->getBlock(block_id);
    const auto& sol              = fields.solution->getData(block_id);
    const auto& cell_inv_volume  = disc->getInvCellVolumeField()->getData(block_id);
    const auto& normals          = disc->getNormalField()->getData(block_id, dir);
    auto& res                    = fields.residual->getData(block_id);

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
          flux      = flux_func(qL, qR, normal);
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

          Vec4<T> qLhalf = recon(qLm1, qL, qR, 1);
          Vec4<T> qRhalf = recon(qL, qR, qRp1, -1);
          flux      = flux_func(qLhalf, qRhalf, normal);
        }

        T inv_volL = cell_inv_volume(face_id.cell_i_left, face_id.cell_j_left, 0);
        T inv_volR = cell_inv_volume(face_id.cell_i_right, face_id.cell_j_right, 0);

        for (UInt k=0; k < DofsPerCell; ++k)
        {
          res(face_id.cell_i_left, face_id.cell_j_left, k)   -=  inv_volL * flux[k];
          res(face_id.cell_i_right, face_id.cell_j_right, k) +=  inv_volR * flux[k];
        }
      }
  }
}


template <typename T, typename Reconstruction, typename Flux>
void evaluateInterfaceTermsAddDirection(Fields<T>& fields, Real t, const StructuredDiscPtr disc, 
                                        const Reconstruction& recon, const Flux& flux_func)
{
  static_assert(IsReconstruction<Reconstruction>);
  static_assert(IsNumericalFlux<Flux>);

  evaluateInterfaceTermImpl(fields, t, disc, recon, flux_func, XDirTag());
  evaluateInterfaceTermImpl(fields, t, disc, recon, flux_func, YDirTag());  
}


template <typename T, typename Reconstruction>
void evaluateInterfaceTermsFluxFunction(
       const EulerOpts& opts, Fields<T>& fields, Real t,
       const StructuredDiscPtr disc, const Reconstruction& recon)
{
  static_assert(IsReconstruction<Reconstruction>);

  if (opts.flux == FluxFunction::Roe)
  {
    RoeFlux flux(opts.roe_efix_delta);
    evaluateInterfaceTermsAddDirection(fields, t, disc, recon, flux);
  } else if (opts.flux == FluxFunction::RoeHH)
  {
    RoeHHFlux flux;
    evaluateInterfaceTermsAddDirection(fields, t, disc, recon, flux);
  } else if (opts.flux == FluxFunction::HLLE)
  {
    HLLEFlux flux;
    evaluateInterfaceTermsAddDirection(fields, t, disc, recon, flux);
  } else if (opts.flux == FluxFunction::LLF)
  {
    LaxFriedrichFlux flux;
    evaluateInterfaceTermsAddDirection(fields, t, disc, recon, flux);
  } else if (opts.flux == FluxFunction::HLLC)
  {
    HLLCFlux flux;
    evaluateInterfaceTermsAddDirection(fields, t, disc, recon, flux);
  } else
  {
    throw std::runtime_error("unhandled FluxFunction enum");
  }  
}

template <typename T, typename SlopeLimiter>
void evaluateInterfaceTermsAddReconstruction(
       const EulerOpts& opts, Fields<T>& fields, Real t,
       const StructuredDiscPtr disc, const SlopeLimiter& limiter)
{
  static_assert(common::IsSlopeLimiter<SlopeLimiter>);

  if (opts.recon == Reconstruction::Conservative) {
    ReconstructionConservative<SlopeLimiter> recon(limiter);  
    evaluateInterfaceTermsFluxFunction(opts, fields, t, disc, recon);
  } else if (opts.recon == Reconstruction::Primitive)
  {
    ReconstructionPrimitive<SlopeLimiter> recon(limiter);
    evaluateInterfaceTermsFluxFunction(opts, fields, t, disc, recon); 
  } else
  {
    throw std::runtime_error("unhandled reconstruction enum " + get_name(opts.recon));
  }
}


template <typename T>
void evaluateInterfaceTermsAddLimiter(const EulerOpts& opts, Fields<T>& fields, Real t, 
                                      const StructuredDiscPtr disc)
{
  if (opts.limiter == common::SlopeLimiter::FirstOrder)
  {
    common::SlopeLimiterFirstOrder limiter;   
    evaluateInterfaceTermsAddReconstruction(opts, fields, t, disc, limiter);
  } else if (opts.limiter == common::SlopeLimiter::MinMod)
  {
    common::SlopeLimiterMinMod limiter;
    evaluateInterfaceTermsAddReconstruction(opts, fields, t, disc, limiter);
  } else if (opts.limiter == common::SlopeLimiter::SuperBee)
  {
    common::SlopeLimiterSuperBee limiter;
    evaluateInterfaceTermsAddReconstruction(opts, fields, t, disc, limiter);
  } else if (opts.limiter == common::SlopeLimiter::VanAlba)
  {
    common::SlopeLimiterVanAlba limiter;
    evaluateInterfaceTermsAddReconstruction(opts, fields, t, disc, limiter);
  } else if (opts.limiter == common::SlopeLimiter::VanLeer)
  {
    common::SlopeLimiterVanLeer limiter;
    evaluateInterfaceTermsAddReconstruction(opts, fields, t, disc, limiter); 
  } else
  {
    throw std::runtime_error("unhandled SlopeLimiter enum " + get_name(opts.limiter));
  }
}

template <typename T>
void evaluateInterfaceTerms(const EulerOpts& opts, Fields<T>& fields, Real t, const StructuredDiscPtr disc)
{
  // Unlike a factory, which takes concrete types and returns a pointer to an abstract
  // type, this code builds up the template arguments one by one until everything
  // has a concrete type.
  // Each function in the chain uses an enum from EulerOpts and figures out the
  // concrete type (similar to a factory), but then calls the next function in the chain.
  evaluateInterfaceTermsAddLimiter(opts, fields, t, disc);
}



// ETI
template void evaluateInterfaceTerms<Real>(const EulerOpts& opts, Fields<Real>& fields, Real t, StructuredDiscPtr disc);
template void evaluateInterfaceTerms<Complex>(const EulerOpts& opts, Fields<Complex>& fields, Real t, StructuredDiscPtr disc);


}
}