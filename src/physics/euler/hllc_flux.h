#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_HLLC_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_HLLC_H

#include "euler_flux.h"
#include "numerical_flux_base.h"
#include "roe_state.h"


namespace structured_fv {
namespace euler {

class HLLCFlux final : public NumericalFlux
{
  public:
    /*constexpr*/ Vec4<Real> operator()(const Vec4<Real>& qL, const Vec4<Real>& qR, 
                                    const Vec2<Real>& normal) const
    {
      Real n_mag2 = dot(normal, normal);
      Real n_mag = std::sqrt(n_mag2);
      Real aL = compute_sos(qL);
      Real aR = compute_sos(qR);

      RoeAvgState avg_state = compute_roe_avg(qL, qR);
      Real aAvg = compute_sos(avg_state);

      Real unL = compute_un(qL, normal);
      Real unR = compute_un(qR, normal);
      Real unAvg = compute_un(avg_state, normal);

      Real sL = std::min(unL - aL*n_mag, unAvg - aAvg*n_mag);
      Real sR = std::max(unR + aR*n_mag, unAvg + aAvg*n_mag);


      //std::cout << "sL = " << sL << std::endl;
      //std::cout << "sR = " << sR << std::endl;

      if (sL > 0)
      {
        return compute_euler_flux(qL, normal);
      } else if (sR < 0)
      {
        return compute_euler_flux(qR, normal);
      } else
      {
        Real pL = compute_pressure(qL);
        Real pR = compute_pressure(qR);
        Real deltaSL = sL - unL;
        Real deltaSR = sR - unR;

        Real s_star = ((pR - pL)*n_mag2 + qL[0]*unL*deltaSL - qR[0]*unR*deltaSR)/(qL[0]*deltaSL - qR[0]*deltaSR);
        //std::cout << "s_star = " << s_star << std::endl;

        Vec4<Real> fL = compute_euler_flux(qL, normal);
        Vec4<Real> fR = compute_euler_flux(qR, normal);

        const auto& qk = s_star > 0 ? qL : qR;
        const auto& fk = s_star > 0 ? fL : fR;
        const auto sk  = s_star > 0 ? sL : sR;
        const auto pk  = s_star > 0 ? pL : pR;
        const auto unk = s_star > 0 ? unL : unR;

        Vec4<Real> d{0, normal[0]/n_mag, normal[1]/n_mag, s_star/n_mag};
        Real pressure_term = sk*(pk*n_mag2 + qk[0]*(sk - unk)*(s_star - unk))/n_mag;

        Vec4<Real> flux{0, 0, 0, 0};
        for (UInt i=0; i < 4; ++i)
          flux[i] = (s_star*(sk*qk[i] - fk[i]) + pressure_term*d[i])/(sk - s_star);

        return flux;
      }
    }
    
};

}
}

#endif