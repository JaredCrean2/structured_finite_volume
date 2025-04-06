#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_ROE_STATE_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_ROE_STATE_H

#include "typedefs.h"
#include "utils/project_defs.h"
#include "euler_flux.h"

namespace structured_fv {
namespace euler {

struct RoeAvgState
{
  Real rho = 0;
  Real u = 0;
  Real v = 0;
  Real H = 0;
};

constexpr RoeAvgState compute_roe_avg(const Vec4<Real>& qL, const Vec4<Real>& qR)
{
  Real sqrt_rho_l = std::sqrt(qL[0]); Real sqrt_rho_r = std::sqrt(qR[0]);
  Real pL = compute_pressure(qL);     Real pR = compute_pressure(qR);
  Real rho_sum = sqrt_rho_l + sqrt_rho_r;

  Real rho_avg = std::sqrt(qL[0]*qR[0]);
  Real u_avg = (qL[1]/sqrt_rho_l + qR[1]/sqrt_rho_r)/rho_sum;
  Real v_avg = (qL[2]/sqrt_rho_l + qR[2]/sqrt_rho_r)/rho_sum;
  Real H_avg = ((qL[3] + pL)/sqrt_rho_l + (qR[3] + pR)/sqrt_rho_r)/rho_sum;

  return RoeAvgState{rho_avg, u_avg, v_avg, H_avg};
}

constexpr Real compute_sos2(const RoeAvgState& state)
{
  Real ke = 0.5*(state.u*state.u + state.v*state.v);
  return Gamma_m1*(state.H - ke);
}

constexpr Real compute_sos(const RoeAvgState& state)
{
  return std::sqrt(compute_sos2(state));
}

constexpr Real compute_pressure(const RoeAvgState& avg_state)
{
  return avg_state.rho * compute_sos2(avg_state)/Gamma;
}

constexpr Real compute_un(const RoeAvgState& avg_state, const Vec2<Real>& normal)
{
  return (avg_state.u*normal[0] + avg_state.v*normal[1]);
}

struct RoeStateTag {};

constexpr Vec4<Real> compute_conservative_variables(const RoeAvgState& avg_state, RoeStateTag)
{
  Real pressure = compute_pressure(avg_state);
  Real E = avg_state.H * avg_state.rho - pressure;
  return {avg_state.rho, avg_state.rho*avg_state.u , avg_state.rho*avg_state.v, E};
}

}
}

#endif