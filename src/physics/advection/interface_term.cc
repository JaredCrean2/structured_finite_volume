#include "interface_term.h"
#include "advection_model.h"
#include "disc/discretization.h"
#include "disc/face_field.h"
#include "disc/elem_field.h"
#include "utils/neighbor_direction.h"


namespace structured_fv {
namespace advection {

template <typename T, typename Flux, typename SlopeLimiter, typename Tag>
void evaluateInterfaceTermsImpl(Fields<T> fields, Real t, StructuredDiscPtr disc,
                                const Flux& flux_func, const SlopeLimiter& limiter, 
                                Tag dir_tag)
{
  constexpr double epsilon = 1e-15;
  NeighborDirection dir = toNeighborDirection(dir_tag);

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
        const auto [cell_im1_left, cell_jm1_left] = increment(dir_tag, face_id.cell_i_left, face_id.cell_j_left, -1);
        const auto [cell_ip1_right, cell_jp1_right] = increment(dir_tag, face_id.cell_i_right, face_id.cell_j_right, 1);
                         
        T qLm1      = sol(cell_im1_left, cell_jm1_left, 0);
        T qL        = sol(face_id.cell_i_left, face_id.cell_j_left, 0);
        T qR        = sol(face_id.cell_i_right, face_id.cell_j_right, 0);
        T qRp1      = sol(cell_ip1_right, cell_jp1_right, 0);

        T rL = (qL - qLm1)/(qR - qL + epsilon);
        T rR = (qR - qL)/(qRp1 - qR + epsilon);
        T slopeL = (qR - qLm1)/2.0;
        T slopeR = (qRp1 - qL)/2.0;

        T qLhalf = qL + 0.5*limiter(rL)*slopeL;
        T qRhalf = qR - 0.5*limiter(rR)*slopeR;

        Vec2<Real> normal{normals(i, j, 0), normals(i, j, 1)};
        T flux      = flux_func(qLhalf, qRhalf, normal);

        res(face_id.cell_i_left, face_id.cell_j_left, 0)   -= cell_inv_volume(face_id.cell_i_left, face_id.cell_j_left, 0) * flux;
        res(face_id.cell_i_right, face_id.cell_j_right, 0) += cell_inv_volume(face_id.cell_i_right, face_id.cell_j_right, 0) * flux;
      }     
  }
}


template <typename T, typename SlopeLimiter>
void evaluateInterfaceTermAddDirection(Fields<T> fields, Real t, StructuredDiscPtr disc,
                                      const FluxFunctionUpwind& flux,
                                      SlopeLimiter limiter)
{
  evaluateInterfaceTermsImpl(fields, t, disc, flux, limiter, XDirTag());
  evaluateInterfaceTermsImpl(fields, t, disc, flux, limiter, YDirTag());
}


template <typename T>
void evaluateInterfaceTerm(const AdvectionOpts& opts, Fields<T> fields, Real t, StructuredDiscPtr disc)
{
  FluxFunctionUpwind flux(opts.adv_velocity);  
  if (opts.limiter == common::SlopeLimiter::FirstOrder)
  {
    common::SlopeLimiterFirstOrder limiter;
    evaluateInterfaceTermAddDirection(fields, t, disc, flux, limiter);
  } else if (opts.limiter == common::SlopeLimiter::MinMod)
  {
    common::SlopeLimiterMinMod limiter;
    evaluateInterfaceTermAddDirection(fields, t, disc, flux, limiter);
  } else if (opts.limiter == common::SlopeLimiter::SuperBee)
  {
    common::SlopeLimiterSuperBee limiter;
    evaluateInterfaceTermAddDirection(fields, t, disc, flux, limiter);
  } else if (opts.limiter == common::SlopeLimiter::VanAlba)
  {
    common::SlopeLimiterVanAlba limiter;
    evaluateInterfaceTermAddDirection(fields, t, disc, flux, limiter);
  } else if (opts.limiter == common::SlopeLimiter::VanLeer)
  {
    common::SlopeLimiterVanLeer limiter;
    evaluateInterfaceTermAddDirection(fields, t, disc, flux, limiter);
  } else
  {
    throw std::runtime_error("unsupported SlopeLimiter: " + common::get_name(opts.limiter));
  }
}

// ETI
template void evaluateInterfaceTerm<Real>(const AdvectionOpts& opts, Fields<Real> fields, Real t, StructuredDiscPtr disc);
template void evaluateInterfaceTerm<Complex>(const AdvectionOpts& opts, Fields<Complex> fields, Real t, StructuredDiscPtr disc);

}
}
