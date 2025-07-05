#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_ADVECTION_INTERFACE_TERM_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_ADVECTION_INTERFACE_TERM_H

#include "advection_model.h"
#include "disc/discretization.h"
#include "disc/face_field.h"
#include "disc/elem_field.h"
#include "utils/neighbor_direction.h"


namespace structured_fv {
namespace advection {

class FluxFunctionUpwind
{
  public:
    FluxFunctionUpwind(const Vec2<Real>& advection_velocity) :
      m_adv_velocity(advection_velocity)
    {}

    template <typename T>
    constexpr T operator()(T qL, T qR, const Vec2<Real>& normal) const
    {
      //TODO: do this more efficiently
      T a_normal = dot(m_adv_velocity, normal);
      return a_normal > 0 ? a_normal * qL : a_normal * qR;
    }

    template <typename T>
    constexpr T operator()(T qL, T qR, const Vec2<Real>& normal, T& flux_dotL, T& flux_dotR) const
    {
      T a_normal = dot(m_adv_velocity, normal);
      T flux = a_normal > 0 ? a_normal * qL : a_normal * qR;
      flux_dotL = a_normal > 0 ? a_normal : 0;
      flux_dotR = a_normal > 0 ? 0 : a_normal;

      return flux;
    }    

  private:
    Vec2<Real> m_adv_velocity;
};

template <typename T>
void evaluateInterfaceTerm(const AdvectionOpts& opts, Fields<T> fields, Real t, StructuredDiscPtr disc);


}
}
#endif