#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_HLLE_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_HLLE_H

#include "euler_flux.h"
#include "roe_state.h"

#include <iostream>

namespace structured_fv {
namespace euler {

class HLLEFlux
{
  public:
    /*constexpr*/ Vec4<Real> operator()(const Vec4<Real>& qL, const Vec4<Real>& qR, 
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

      //TODO: optimize this: we only need the most negative and most positive
      //      wave speeds
      Real lambda1L = std::min(unL - aL*n_mag, unAvg - aAvg*n_mag);
      Real lambda2L = std::min(unL, unAvg);
      Real lambda3L = std::min(unL + aL*n_mag, unAvg + aAvg*n_mag);

      Real lambda1R = std::max(unR - aR*n_mag, unAvg - aAvg*n_mag);
      Real lambda2R = std::max(unR, unAvg);
      Real lambda3R = std::max(unR + aR*n_mag, unAvg + aAvg*n_mag);

      Real s1 = std::min(lambda1L, std::min(lambda2L, lambda3L));
      Real s2 = std::max(lambda1R, std::max(lambda2R, lambda3R));


      Vec4<Real> fL = compute_euler_flux(qL, normal);
      Vec4<Real> fR = compute_euler_flux(qR, normal);

      
      if (s1 < 0 && s2 > 0)
      {
        return (s2*fL - s1*fR + s1*s2*(qR - qL))/(s2 - s1);
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