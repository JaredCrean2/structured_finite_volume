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

      //std::cout << "s1 = " << s1 << ", s2 = " << s2 << std::endl;

      Vec4<Real> fL = compute_euler_flux(qL, normal);
      Vec4<Real> fR = compute_euler_flux(qR, normal);

      //TODO: optimize this
      Vec4<Real> q_star = (fR - fL - s2*qR + s1*qL)/(s1 - s2);
      //std::cout << "qL = " << qL << std::endl;
      //std::cout << "qR = " << qR << std::endl;
      //std::cout << "q_star = " << q_star << std::endl;

      if (s1 > 0)
      {
        //std::cout << "in region 1" << std::endl;
        return compute_euler_flux(qL, normal);
      } else if (s2 < 0)
      {
        //std::cout << "in region 3" << std::endl;
        return compute_euler_flux(qR, normal);
      } else
      {
        //std::cout << "in region 2" << std::endl;
        return compute_euler_flux(q_star, normal);
      }
    }
};

}
}

#endif