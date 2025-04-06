#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_FLUX_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_FLUX_H

#include "typedefs.h"

namespace structured_fv {
namespace euler {

// the solution vector is [rho, rho u, rho v, E]
// where E is the total energy per unit volume
// Some constituitive relations:
//   - e is the specific internal energy (energy per unit mass)
//   - rho e is the internal energy per unit volume
//   - 0.5(u^2 + v^2) is specific kinetic energy (kinetic energy per unit mass)
//   - 0.5 rho(u^2 + v^2) is kinetic energy per unit volume
//   - E = rho e + 0.5 rho (u^2 + v^2)
//   - e = cv T
//   - h = cp T =  e + p/rho (specific enthaly aka enthalpy per unit mass)
//   - gamma = cp/cv
//   - e = cv T = cv (p/(rho R)) = (cv/R)(p/rho) = p/((gamma-1) rho)
//   - this E = rho e + 0.5 rho (u^2 + v^2) = p/(gamma-1) + 0.5 rho(u^2 + v^2)
// Note that some authors define E to be total energy per unit mass and
// have the 4th component of the solution vector = rho E
constexpr Real compute_pressure(const Vec4<Real>& sol)
{
  return Gamma_m1 * (sol[3] - 0.5*(sol[1]*sol[1] + sol[2]*sol[2])/sol[0]);
}

constexpr Real compute_sos2(const Vec4<Real>& sol)
{
  return Gamma * compute_pressure(sol)/sol[0];

}

constexpr Real compute_sos(const Vec4<Real>& sol)
{
  return std::sqrt(compute_sos2(sol));
}

constexpr Real compute_un(const Vec4<Real>& sol, const Vec2<Real>& normal)
{
  return (sol[1]*normal[0] + sol[2]*normal[1])/sol[0];
}

constexpr Real compute_temperature(const Vec4<Real>& sol)
{
  Real rho_e = sol[3] - 0.5*(sol[1]*sol[1] + sol[2]*sol[2])/sol[0];
  return rho_e/(sol[0]*Cv);
}

constexpr Vec4<Real> compute_euler_flux(const Vec4<Real>& sol, const Vec2<Real>& normal)
{
  Real pressure = compute_pressure(sol);
  Real momentum_n = (sol[1]*normal[0] + sol[2]*normal[1]);
  Real velocity_n = momentum_n/sol[0];

  Vec4<Real> flux{0, 0, 0, 0};
  flux[0] = momentum_n;
  flux[1] = sol[1]*velocity_n + pressure*normal[0];
  flux[2] = sol[2]*velocity_n + pressure*normal[1];
  flux[3] = (sol[3] + pressure)*velocity_n;
  
  return flux;
}

struct PrimitiveVarTag {};

constexpr Vec4<Real> compute_conservative_variables(const Vec4<Real>& primitive_vars, PrimitiveVarTag)
{
  Real e = Cv*primitive_vars[3];
  Real E = primitive_vars[0]*(e + 0.5*(primitive_vars[1]*primitive_vars[1] + primitive_vars[2]*primitive_vars[2]));

  return {primitive_vars[0], primitive_vars[0]*primitive_vars[1], primitive_vars[0]*primitive_vars[2], E};
}

// For now, rotate solution into the x direction and rotate the flux
// back to the wall normal direction.
// Eventually work out the eigen decomposition of the flux jacobian in the
// wall normal direction 
/*
class Rotator
{
  public:
    constexpr Rotator(const Vec2<Real>& normal) :
      m_normal(normal),
      m_x_normal{std::sqrt(normal[0]*normal[0] + normal[1]*normal[1]), 0}
    {}

    constexpr Vec4<Real> rotateForward(const Vec4<Real>& sol) const
    {
      Real mag = m_x_normal[0];
      return {sol[0], (m_normal[0]*sol[1] + m_normal[1]*sol[1])/mag,
              (-m_normal[1]*sol[1] + m_normal[0]*sol[2])/mag, sol[3]};
    }

    constexpr Vec4<Real> rotateBack(const Vec4<Real>& flux) const
    {
      Real mag = m_x_normal[0];
      return {flux[0], (m_normal[0]*flux[1] - m_normal[1]*flux[2])/mag,
              (m_normal[1]*flux[1] + m_normal[0]*flux[0])/mag, flux[3]};
    }

    constexpr const Vec2<Real>& getRotatedNormal() const { return m_x_normal; }

  private:
    const Vec2<Real> m_normal;
    const Vec2<Real> m_x_normal;
};
*/




}
}

#endif