#include "gtest/gtest.h"
#include "utils/math.h"

namespace {

using namespace structured_fv;
//template <typename T>
//using Vec3 = structured_fv::Vec3<T>;

}


TEST(Math, ComplexIntOps)
{
  Complex a(2, 3);
  Int b = 4;
  UInt c = 5;
  Real br = b;
  Real cr = c;
  Complex d(4, 3);

  EXPECT_EQ(a + b, a  + br);
  EXPECT_EQ(b + a, br + a);
  EXPECT_EQ(a + c, a  + cr);
  EXPECT_EQ(c + a, cr + a);


  EXPECT_EQ(a - b, a  - br);
  EXPECT_EQ(b - a, br - a);
  EXPECT_EQ(a - c, a  - cr);
  EXPECT_EQ(c - a, cr - a);

  EXPECT_EQ(a * b, a  * br);
  EXPECT_EQ(b * a, br * a);
  EXPECT_EQ(a * c, a  * cr);
  EXPECT_EQ(c * a, cr * a);

  EXPECT_EQ(a / b, a  / br);
  EXPECT_EQ(b / a, br / a);
  EXPECT_EQ(a / c, a  / cr);
  EXPECT_EQ(c / a, cr / a);  

  EXPECT_EQ(std::min(a, d), a);
  EXPECT_EQ(std::max(a, d), d);
}

TEST(Math, ComplexCompare)
{
  Complex one(1.0, 2.0), two(2.0, 1.0);
  EXPECT_TRUE(one <   two);
  EXPECT_TRUE(one <=  two);
  EXPECT_FALSE(one >  two);
  EXPECT_FALSE(one >= two);

  EXPECT_TRUE(one <= one);
  EXPECT_TRUE(one >= one);

  EXPECT_FALSE(two <   one);
  EXPECT_FALSE(two <=  one);
  EXPECT_TRUE(two >  one);
  EXPECT_TRUE(two >= one);  
} 



TEST(Math, SmoothAbs)
{
  double delta = 0.1;
  EXPECT_EQ(smoothAbs(5.0, delta), std::abs(5.0));
  EXPECT_EQ(smoothAbs(-5.0, delta), std::abs(-5.0));
  EXPECT_NEAR(smoothAbs(0.05, delta), 0.05*0.05/delta, 1e-13);
  EXPECT_NEAR(smoothAbs(-0.05, delta), 0.05*0.05/delta, 1e-13);


  EXPECT_EQ(smoothAbs(Complex(5, 10), delta), Complex(5, 10));
  EXPECT_EQ(smoothAbs(Complex(-5, 10), delta), Complex(5, -10));
}

TEST(Math, SmoothAbsDot)
{
  Real h = 1e-40;
  Complex pert(0, h);

  double delta = 0.1;
  for (double val : {-2.0, -0.05, 0.05, 2.0})
  {

    Complex x = val + pert;
    Complex abs_xc = smoothAbs(x, delta);
    Real x_dotc = abs_xc.imag()/h;

    auto [abs_x, abs_x_dot] = smoothAbs_dot(val, 2, delta);
    EXPECT_NEAR(2.0*x_dotc, abs_x_dot, 1e-13);
  }
}

TEST(Math, RadianConversion)
{
  EXPECT_NEAR(degreesToRadians(30), PI/6, 1e-13);
  EXPECT_NEAR(radiansToDegrees(PI/6), 30, 1e-13);
}

TEST(Math, Sgn)
{
  EXPECT_EQ(sgn(2), 1);
  EXPECT_EQ(sgn(0), 0);
  EXPECT_EQ(sgn(-2), -1);

  EXPECT_EQ(sgn(2.0), 1);
  EXPECT_EQ(sgn(0.0), 0);
  EXPECT_EQ(sgn(-2.0), -1);  

  EXPECT_EQ(sgn(Complex(2, 5)), 1);
  EXPECT_EQ(sgn(Complex(0, 5)), 0);
  EXPECT_EQ(sgn(Complex(-2, 5)), -1);  
}

TEST(Math, QuadArea)
{
  FixedVec<FixedVec<Real, 2>, 4> coords;
  coords[0] = {0, 0};
  coords[1] = {1, 0};
  coords[2] = {1, 1};
  coords[3] = {0, 1};
  EXPECT_EQ(computeQuadArea(coords), 1);
}

TEST(Math, QuadAreaNegative)
{
  FixedVec<FixedVec<Real, 2>, 4> coords;
  coords[0] = {0, 0};
  coords[1] = {0, 1};
  coords[2] = {1, 1};
  coords[3] = {1, 0};

  EXPECT_EQ(computeQuadArea(coords), -1);
}