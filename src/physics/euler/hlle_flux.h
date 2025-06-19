#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_HLLE_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_HLLE_H

#include "euler_flux.h"
#include "numerical_flux_base.h"
#include "roe_state.h"

namespace structured_fv {
namespace euler {

class HLLEFlux final : public NumericalFlux
{
  public:
    constexpr Vec4<Real> operator()(const Vec4<Real>& qL, const Vec4<Real>& qR, 
                                    const Vec2<Real>& normal) const
    {
      Real n_mag = std::sqrt(dot(normal, normal));
      Real aL = compute_sos(qL);
      Real aR = compute_sos(qR);

      RoeAvgState avg_state = compute_roe_avg(qL, qR);
      Real aAvg = compute_sos(avg_state);

      Real unL = compute_un(qL, normal);
      Real unR = compute_un(qR, normal);
      Real unAvg = compute_un(avg_state, normal);

      Real s1 = std::min(unL - aL*n_mag, unAvg - aAvg*n_mag);
      Real s2 = std::max(unR + aR*n_mag, unAvg + aAvg*n_mag);

      Vec4<Real> fL = compute_euler_flux(qL, normal);
      Vec4<Real> fR = compute_euler_flux(qR, normal);

      if (s1 < 0 && s2 > 0)
      {
        Vec4<Real> flux{};
        for (UInt i=0; i < 4; ++i)
          flux[i] = (s2*fL[i] - s1*fR[i] + s1*s2*(qR[i] - qL[i]))/(s2 - s1);
        return flux;
      } else if (s1 > 0)
      {
        return fL;
      } else  // s1 < s2 < 0
      {
        return fR;
      }
    }
    
};

}
}

#endif