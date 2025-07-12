#include "gtest/gtest.h"
#include <stdexcept>
#include "physics/common/slope_limiters.h"
#include "physics/euler/euler_flux.h"
#include "physics/euler/reconstruction.h"
#include "jacobian.h"

namespace {
using namespace structured_fv;
using namespace structured_fv::euler;
}
// all slope limiters exactly reproduce linear functions

TEST(Reconstruction, Conservative)
{
  Vec4<Real> prim_vars = {2.0, 10.0, 20.0, 300.0};
  Vec4<Real> q_im1 = compute_conservative_variables(prim_vars, PrimitiveVarTag{});
  Vec4<Real> delta_q{0.1, 0.2, 0.3, 0.4};
  Vec4<Real> q_i = q_im1 + delta_q;
  Vec4<Real> q_ip1 = q_im1 + 2*delta_q;

  using Limiter = common::SlopeLimiterMinMod;
  Limiter limiter;
  ReconstructionConservative<Limiter> recon(limiter);

  Vec4<Real> qR = recon(q_im1, q_i, q_ip1, 1);
  Vec4<Real> qL = recon(q_im1, q_i, q_ip1, -1);

  for (UInt i=0; i < 4; ++i)
  {
    EXPECT_DOUBLE_EQ(qR[i], q_i[i] + delta_q[i]/2);
    EXPECT_DOUBLE_EQ(qL[i], q_i[i] - delta_q[i]/2);
  }
}

TEST(Reconstruction, ConservativeJac)
{
  Vec4<Real> prim_vars = {2.0, 10.0, 20.0, 300.0};
  Vec4<Real> q_im1 = compute_conservative_variables(prim_vars, PrimitiveVarTag{});
  Vec4<Real> delta_q{0.1, 0.2, 0.3, 0.4};
  Vec4<Real> q_i = q_im1 + delta_q;
  Vec4<Real> q_ip1 = q_im1 + 2*delta_q;

  using Limiter = common::SlopeLimiterVanAlba;
  Limiter limiter;
  ReconstructionConservative<Limiter> recon(limiter);

  UInt argument = 0;

  auto func = [&](auto q)
  {
    using VecType = decltype(q);
    using ValueType = typename VecType::value_type;
    VecType q_im1_local = copy<ValueType>(q_im1);
    VecType q_i_local = copy<ValueType>(q_i);
    VecType q_ip1_local = copy<ValueType>(q_ip1);

    if (argument == 0)
      q_im1_local = q;
    else if (argument == 1)
      q_i_local = q;
    else if (argument == 2)
      q_ip1_local = q;
    else
      throw std::runtime_error("unhandled value");
      

    return recon(q_im1_local, q_i_local, q_ip1_local, -1);
  };

  auto jac = [&](auto q)
  { 
    Vec4<Real> q_im1_local = q_im1;
    Vec4<Real> q_i_local = q_i;
    Vec4<Real> q_ip1_local = q_ip1;

    if (argument == 0)
      q_im1_local = q;
    else if (argument == 1)
      q_i_local = q;
    else if (argument == 2)
      q_ip1_local = q;
    else
      throw std::runtime_error("unhandled value");    

    Matrix<Real, 4> q_im1_jac, q_i_jac, q_ip1_jac;
    auto qL = recon(q_im1_local, q_im1_jac, q_i_local, q_i_jac, q_ip1_local, q_ip1_jac, -1);


    Matrix<Real, 4> jac_out;
    if (argument == 0)
      jac_out = q_im1_jac;
    else if (argument == 1)
      jac_out = q_i_jac;
    else if (argument == 2)
      jac_out = q_ip1_jac;
    else
      throw std::runtime_error("unhandled value");

    std::cout << "jac_out = \n" << jac_out << std::endl;
    return std::make_pair(qL, jac_out);
  };

  test_utils::checkJacobianVector(q_im1, func, jac);

  argument = 1;
  test_utils::checkJacobianVector(q_i, func, jac);

  argument = 2;
  test_utils::checkJacobianVector(q_ip1, func, jac);
}

TEST(Reconstruction, ConservativeLeftSlope)
{
  Vec4<Real> prim_vars = {2.0, 10.0, 20.0, 300.0};
  Vec4<Real> q_i = compute_conservative_variables(prim_vars, PrimitiveVarTag{});
  Vec4<Real> delta_q{0.1, 0.2, 0.3, 0.4};
  Vec4<Real> q_im1 = q_i - 0.5*delta_q;
  Vec4<Real> q_ip1 = q_i + delta_q;

  using Limiter = common::SlopeLimiterMinMod;
  Limiter limiter;
  ReconstructionConservative<Limiter> recon(limiter);

  Vec4<Real> qR = recon(q_im1, q_i, q_ip1, 1);
  Vec4<Real> qL = recon(q_im1, q_i, q_ip1, -1);

  for (UInt i=0; i < 4; ++i)
  {
    EXPECT_DOUBLE_EQ(qR[i], q_i[i] + 0.5*delta_q[i]/2);
    EXPECT_DOUBLE_EQ(qL[i], q_i[i] - 0.5*delta_q[i]/2);
  }
}

TEST(Reconstruction, ConservativeRightSlope)
{
  Vec4<Real> prim_vars = {2.0, 10.0, 20.0, 300.0};
  Vec4<Real> q_i = compute_conservative_variables(prim_vars, PrimitiveVarTag{});
  Vec4<Real> delta_q{0.1, 0.2, 0.3, 0.4};
  Vec4<Real> q_im1 = q_i - delta_q;
  Vec4<Real> q_ip1 = q_i + 0.5*delta_q;

  using Limiter = common::SlopeLimiterMinMod;
  Limiter limiter;
  ReconstructionConservative<Limiter> recon(limiter);

  Vec4<Real> qR = recon(q_im1, q_i, q_ip1, 1);
  Vec4<Real> qL = recon(q_im1, q_i, q_ip1, -1);

  for (UInt i=0; i < 4; ++i)
  {
    EXPECT_DOUBLE_EQ(qR[i], q_i[i] + 0.5*delta_q[i]/2);
    EXPECT_DOUBLE_EQ(qL[i], q_i[i] - 0.5*delta_q[i]/2);
  }
}

TEST(Reconstruction, Primitive)
{
  Vec4<Real> prim_i = {2.0, 10.0, 20.0, 300.0};
  Vec4<Real> delta_p{0.1, 0.2, 0.3, 0.4};

  Vec4<Real> q_i = compute_conservative_variables(prim_i);
  Vec4<Real> q_im1 = compute_conservative_variables(prim_i - delta_p);
  Vec4<Real> q_ip1 = compute_conservative_variables(prim_i + delta_p);

  using Limiter = common::SlopeLimiterMinMod;
  Limiter limiter;
  ReconstructionPrimitive<Limiter> recon(limiter);

  Vec4<Real> qR = recon(q_im1, q_i, q_ip1, 1);
  Vec4<Real> qL = recon(q_im1, q_i, q_ip1, -1);

  Vec4<Real> primL = compute_primitive_variables(qL);
  Vec4<Real> primR = compute_primitive_variables(qR);


  for (UInt i=0; i < 4; ++i)
  {
    EXPECT_DOUBLE_EQ(primR[i], prim_i[i] + delta_p[i]/2);
    EXPECT_DOUBLE_EQ(primL[i], prim_i[i] - delta_p[i]/2);
  }
}

TEST(Reconstruction, PrimitiveJac)
{
  Vec4<Real> prim_i = {2.0, 10.0, 20.0, 300.0};
  Vec4<Real> delta_p{0.1, 0.2, 0.3, 0.4};

  Vec4<Real> q_i = compute_conservative_variables(prim_i);
  Vec4<Real> q_im1 = compute_conservative_variables(prim_i - delta_p);
  Vec4<Real> q_ip1 = compute_conservative_variables(prim_i + delta_p);

  using Limiter = common::SlopeLimiterVanAlba;
  Limiter limiter;
  ReconstructionPrimitive<Limiter> recon(limiter);

  UInt argument = 0;

  auto func = [&](auto q)
  {
    using VecType = decltype(q);
    using ValueType = typename VecType::value_type;
    VecType q_im1_local = copy<ValueType>(q_im1);
    VecType q_i_local = copy<ValueType>(q_i);
    VecType q_ip1_local = copy<ValueType>(q_ip1);

    if (argument == 0)
      q_im1_local = q;
    else if (argument == 1)
      q_i_local = q;
    else if (argument == 2)
      q_ip1_local = q;
    else
      throw std::runtime_error("unhandled value");
      

    return recon(q_im1_local, q_i_local, q_ip1_local, -1);
  };

  auto jac = [&](auto q)
  { 
    Vec4<Real> q_im1_local = q_im1;
    Vec4<Real> q_i_local = q_i;
    Vec4<Real> q_ip1_local = q_ip1;

    if (argument == 0)
      q_im1_local = q;
    else if (argument == 1)
      q_i_local = q;
    else if (argument == 2)
      q_ip1_local = q;
    else
      throw std::runtime_error("unhandled value");    

    Matrix<Real, 4> q_im1_jac, q_i_jac, q_ip1_jac;
    auto qL = recon(q_im1_local, q_im1_jac, q_i_local, q_i_jac, q_ip1_local, q_ip1_jac, -1);


    Matrix<Real, 4> jac_out;
    if (argument == 0)
      jac_out = q_im1_jac;
    else if (argument == 1)
      jac_out = q_i_jac;
    else if (argument == 2)
      jac_out = q_ip1_jac;
    else
      throw std::runtime_error("unhandled value");

    return std::make_pair(qL, jac_out);
  };

  test_utils::checkJacobianVector(q_im1, func, jac);

  argument = 1;
  test_utils::checkJacobianVector(q_i, func, jac);

  argument = 2;
  test_utils::checkJacobianVector(q_ip1, func, jac);
}

TEST(Reconstruction, PrimitiveLeftSlope)
{
  Vec4<Real> prim_i = {2.0, 10.0, 20.0, 300.0};
  Vec4<Real> delta_p{0.1, 0.2, 0.3, 0.4};

  Vec4<Real> q_i = compute_conservative_variables(prim_i);
  Vec4<Real> q_im1 = compute_conservative_variables(prim_i - delta_p);
  Vec4<Real> q_ip1 = compute_conservative_variables(prim_i + 0.5*delta_p);

  using Limiter = common::SlopeLimiterMinMod;
  Limiter limiter;
  ReconstructionPrimitive<Limiter> recon(limiter);

  Vec4<Real> qR = recon(q_im1, q_i, q_ip1, 1);
  Vec4<Real> qL = recon(q_im1, q_i, q_ip1, -1);

  Vec4<Real> primL = compute_primitive_variables(qL);
  Vec4<Real> primR = compute_primitive_variables(qR);


  for (UInt i=0; i < 4; ++i)
  {
    EXPECT_DOUBLE_EQ(primR[i], prim_i[i] + 0.5*delta_p[i]/2);
    EXPECT_DOUBLE_EQ(primL[i], prim_i[i] - 0.5*delta_p[i]/2);
  }
}

TEST(Reconstruction, PrimitiveRightSlope)
{
  Vec4<Real> prim_i = {2.0, 10.0, 20.0, 300.0};
  Vec4<Real> delta_p{0.1, 0.2, 0.3, 0.4};

  Vec4<Real> q_i = compute_conservative_variables(prim_i);
  Vec4<Real> q_im1 = compute_conservative_variables(prim_i - 0.5*delta_p);
  Vec4<Real> q_ip1 = compute_conservative_variables(prim_i + delta_p);

  using Limiter = common::SlopeLimiterMinMod;
  Limiter limiter;
  ReconstructionPrimitive<Limiter> recon(limiter);

  Vec4<Real> qR = recon(q_im1, q_i, q_ip1, 1);
  Vec4<Real> qL = recon(q_im1, q_i, q_ip1, -1);

  Vec4<Real> primL = compute_primitive_variables(qL);
  Vec4<Real> primR = compute_primitive_variables(qR);


  for (UInt i=0; i < 4; ++i)
  {
    EXPECT_DOUBLE_EQ(primR[i], prim_i[i] + 0.5*delta_p[i]/2);
    EXPECT_DOUBLE_EQ(primL[i], prim_i[i] - 0.5*delta_p[i]/2);
  }
}

TEST(Reconstruction, OutputOperator)
{
  using Limiter = common::SlopeLimiterMinMod;
  Limiter limiter;
  ReconstructionConservative<Limiter> recon_conservative(limiter);
  ReconstructionPrimitive<Limiter> recon_primitive(limiter);

  {
    std::stringstream ss;
    ss << recon_conservative;
    EXPECT_EQ(ss.str(), "Conservative");
  }

  {
    std::stringstream ss;
    ss << recon_primitive;
    EXPECT_EQ(ss.str(), "Primitive");
  }
}
