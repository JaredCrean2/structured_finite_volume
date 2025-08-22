#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_ROE_STATE_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_ROE_STATE_H

#include "euler_flux.h"
#include "typedefs.h"
#include "utils/project_defs.h"
#include <iostream>

namespace structured_fv {
namespace euler {

template <typename T>
struct RoeAvgState {
  T rho = 0;
  T u = 0;
  T v = 0;
  T H = 0;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const RoeAvgState<T>& state)
{
  os << "rho = " << state.rho << ", u = " << state.u << ", v = " << state.v << ", H = " << state.H;
  return os;
}


template <typename T>
constexpr RoeAvgState<T> compute_roe_avg(const Vec4<T> &qL, const Vec4<T> &qR)
{
  T sqrt_rho_l = std::sqrt(qL[0]);
  T sqrt_rho_r = std::sqrt(qR[0]);
  T pL = compute_pressure(qL);
  T pR = compute_pressure(qR);
  T rho_sum = sqrt_rho_l + sqrt_rho_r;

  T rho_avg = std::sqrt(qL[0] * qR[0]);
  T u_avg = (qL[1] / sqrt_rho_l + qR[1] / sqrt_rho_r) / rho_sum;
  T v_avg = (qL[2] / sqrt_rho_l + qR[2] / sqrt_rho_r) / rho_sum;
  T H_avg = ((qL[3] + pL) / sqrt_rho_l + (qR[3] + pR) / sqrt_rho_r) / rho_sum;

  return RoeAvgState<T>{rho_avg, u_avg, v_avg, H_avg};
}

template <typename T>
constexpr RoeAvgState<T>
compute_roe_avg_jac(const Vec4<T> &qL, const Vec4<T> &qR,
                    Matrix<T, 4> &roe_avg_dotL, Matrix<T, 4> &roe_avg_dotR)
{
  T sqrt_rho_l = std::sqrt(qL[0]);
  T sqrt_rho_r = std::sqrt(qR[0]);
  T sqrt_rho_l_dot = 1.0 / (2 * sqrt_rho_l);
  T sqrt_rho_r_dot = 1.0 / (2 * sqrt_rho_r);

  auto [pL, pL_dot] = compute_pressure_jac(qL);
  auto [pR, pR_dot] = compute_pressure_jac(qR);
  T rho_sum = sqrt_rho_l + sqrt_rho_r;
  T rho_sum_dotL0 = sqrt_rho_l_dot;
  T rho_sum_dotR0 = sqrt_rho_r_dot;

  T rho_avg = std::sqrt(qL[0] * qR[0]);
  T rho_avg_dotL0 = qR[0] / (2 * rho_avg);
  T rho_avg_dotR0 = qL[0] / (2 * rho_avg);

  T t1 = qL[1] / sqrt_rho_l;
  T t1_dotL0 = -qL[1] * sqrt_rho_l_dot / qL[0];
  T t1_dotL1 = 1.0/sqrt_rho_l;

  T t2 = qR[1] / sqrt_rho_r;
  T t2_dotR0 = -qR[1] * sqrt_rho_r_dot / (qR[0]);
  T t2_dotR1 = 1.0/sqrt_rho_r;
  // T u_avg = (qL[1]/sqrt_rho_l + qR[1]/sqrt_rho_r)/rho_sum;

  T u_avg = (t1 + t2) / rho_sum;
  T u_avg_dotL0 = t1_dotL0/rho_sum - (t1 + t2) * rho_sum_dotL0/(rho_sum*rho_sum);
  T u_avg_dotR0 = t2_dotR0/rho_sum - (t1 + t2) * rho_sum_dotR0/(rho_sum*rho_sum);
  T u_avg_dotL1 = t1_dotL1/rho_sum;
  T u_avg_dotR1 = t2_dotR1/rho_sum;

  T t3 = qL[2] / sqrt_rho_l;
  T t3_dotL0 = -qL[2]*sqrt_rho_l_dot/(qL[0]);
  T t3_dotL2 = 1.0/sqrt_rho_l;

  T t4 = qR[2] / sqrt_rho_r;
  T t4_dotR0 = -qR[2]*sqrt_rho_r_dot/(qR[0]);
  T t4_dotR2 =  1.0/sqrt_rho_r;
  //T v_avg = (qL[2] / sqrt_rho_l + qR[2] / sqrt_rho_r) / rho_sum;
  T v_avg = (t3 + t4)/rho_sum;
  T v_avg_dotL0 = t3_dotL0/rho_sum - (t3 + t4)*rho_sum_dotL0/(rho_sum*rho_sum);
  T v_avg_dotR0 = t4_dotR0/rho_sum - (t3 + t4)*rho_sum_dotR0/(rho_sum*rho_sum);
  T v_avg_dotL2 = t3_dotL2/rho_sum;
  T v_avg_dotR2 = t4_dotR2/rho_sum;

  T t5 = (qL[3] + pL) / sqrt_rho_l;
  Vec4<T> t5_dotL = pL_dot / sqrt_rho_l;
  t5_dotL[3] += 1.0/sqrt_rho_l;
  t5_dotL[0] -= (qL[3] + pL) * sqrt_rho_l_dot / qL[0];

  T t6 = (qR[3] + pR) / sqrt_rho_r;
  Vec4<T> t6_dotR = pR_dot / sqrt_rho_r;
  t6_dotR[3] += 1.0/sqrt_rho_r;
  t6_dotR[0] -= (qR[3] + pR) * sqrt_rho_r_dot / qR[0];

  //T H_avg = ((qL[3] + pL) / sqrt_rho_l + (qR[3] + pR) / sqrt_rho_r) / rho_sum;
  T H_avg = (t5 + t6) / rho_sum;
  Vec4<T> H_avg_dotL = t5_dotL / rho_sum;
  H_avg_dotL[0] -= (t5 + t6) * rho_sum_dotL0 / (rho_sum*rho_sum);
  Vec4<T> H_avg_dotR = t6_dotR / rho_sum;
  H_avg_dotR[0] -= (t5 + t6) * rho_sum_dotR0 / (rho_sum*rho_sum);


  RoeAvgState<T> avg_state{rho_avg, u_avg, v_avg, H_avg};
  roe_avg_dotL = {rho_avg_dotL0, 0, 0, 0,
                  u_avg_dotL0, u_avg_dotL1, 0, 0,
                  v_avg_dotL0, 0, v_avg_dotL2, 0,
                  H_avg_dotL[0], H_avg_dotL[1], H_avg_dotL[2], H_avg_dotL[3]};

  roe_avg_dotR = {rho_avg_dotR0, 0, 0, 0,
                  u_avg_dotR0, u_avg_dotR1, 0, 0,
                  v_avg_dotR0, 0, v_avg_dotR2, 0,
                  H_avg_dotR[0], H_avg_dotR[1], H_avg_dotR[2], H_avg_dotR[3]}; 

  return avg_state;                 
}

template <typename T> constexpr T compute_sos2(const RoeAvgState<T> &state)
{
  T ke = 0.5 * (state.u * state.u + state.v * state.v);
  return Gamma_m1 * (state.H - ke);
}

template <typename T>
constexpr ScalarVectorPair<T> compute_sos2_jac(const RoeAvgState<T> &state)
{
  T ke = 0.5 * (state.u * state.u + state.v * state.v);
  return {Gamma_m1 * (state.H - ke),
          {0, -Gamma_m1 * state.u, -Gamma_m1 * state.v, Gamma_m1}};
}

template <typename T> 
constexpr T compute_sos(const RoeAvgState<T> &state)
{
  return std::sqrt(compute_sos2(state));
}

template <typename T>
constexpr ScalarVectorPair<T> compute_sos_jac(const RoeAvgState<T> &state)
{
  auto [sos2, sos2_jac] = compute_sos2_jac(state);
  T sos = std::sqrt(sos2);

  for (UInt i = 0; i < 4; ++i)
    sos2_jac[i] /= 2 * sos;

  return {sos, sos2_jac};
}

template <typename T>
constexpr T compute_pressure(const RoeAvgState<T> &avg_state)
{
  return avg_state.rho * compute_sos2(avg_state) / Gamma;
}

template <typename T>
constexpr ScalarVectorPair<T>
compute_pressure_jac(const RoeAvgState<T> &avg_state)
{
  auto [sos2, sos2_dot] = compute_sos2_jac(avg_state);

  for (UInt i = 0; i < 4; ++i)
    sos2_dot[i] *= avg_state.rho / Gamma;
  sos2_dot[0] += sos2 / Gamma;

  return {avg_state.rho * compute_sos2(avg_state) / Gamma, sos2_dot};
}

template <typename T>
constexpr T compute_un(const RoeAvgState<T> &avg_state,
                       const Vec2<Real> &normal)
{
  return (avg_state.u * normal[0] + avg_state.v * normal[1]);
}

template <typename T>
constexpr ScalarVectorPair<T> compute_un_jac(const RoeAvgState<T> &avg_state,
                                             const Vec2<Real> &normal)
{
  return {(avg_state.u * normal[0] + avg_state.v * normal[1]),
          {0, normal[0], normal[1], 0}};
}

struct RoeStateTag {};

template <typename T>
constexpr Vec4<T>
compute_conservative_variables(const RoeAvgState<T> &avg_state, RoeStateTag)
{
  T pressure = compute_pressure(avg_state);
  T E = avg_state.H * avg_state.rho - pressure;
  return {avg_state.rho, avg_state.rho * avg_state.u,
          avg_state.rho * avg_state.v, E};
}

template <typename T>
constexpr VectorMatrixPair<T>
compute_conservative_variables_jac(const RoeAvgState<T> &avg_state,
                                   RoeStateTag) {
  auto [pressure, pressure_dot] = compute_pressure_jac(avg_state);
  T E = avg_state.H * avg_state.rho - pressure;
  Vec4<T> E_dot = -pressure_dot;
  E_dot[0] += avg_state.H;
  E_dot[3] += avg_state.rho;

  Vec4<T> q = {avg_state.rho, avg_state.rho * avg_state.u,
               avg_state.rho * avg_state.v, E};
  Matrix<T, 4> q_dot({1, 0, 0, 0,
                      avg_state.u, avg_state.rho, 0, 0,
                      avg_state.v, 0, avg_state.rho, 0,
                      E_dot[0], E_dot[1], E_dot[2], E_dot[3]});
  return {q, q_dot};
}

} // namespace euler
} // namespace structured_fv

#endif