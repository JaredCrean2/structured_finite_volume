#include "gtest/gtest.h"
#include "physics/common/slope_limiters.h"
#include "physics/euler/euler_flux.h"
#include "physics/euler/reconstruction.h"

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
