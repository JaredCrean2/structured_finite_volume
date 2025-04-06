#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_LFFLUX_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_LFFLUX_H

#include "euler_flux.h"

namespace structured_fv {
namespace euler {

class LaxFriedrichFlux
{
  public:
    constexpr Vec4<Real> operator()(const Vec4<Real>& qL, const Vec4<Real>& qR, 
                                    const Vec2<Real>& normal) const
    {
      Real n_mag = std::sqrt(dot(normal, normal));
      Real aL = compute_sos(qL);
      Real aR = compute_sos(qR);
      Real unL = compute_un(qL, normal);
      Real unR = compute_un(qR, normal);

      Real lambda_max = std::max(std::abs(unL) + aL*n_mag, std::abs(unR) + aR*n_mag);
      Vec4<Real> fL = compute_euler_flux(qL, normal);
      Vec4<Real> fR = compute_euler_flux(qR, normal);

      Vec4<Real> flux{0, 0, 0, 0};
      for (UInt i=0; i < 4; ++i)
        flux[i] = 0.5*(fL[i] + fR[i]) - 0.5*lambda_max*(qR[i] - qL[i]);

      return flux;
    }
};

}
}

#endif