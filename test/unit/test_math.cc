#include "gtest/gtest.h"
#include "utils/math.h"

namespace {

using namespace structured_fv;
//template <typename T>
//using Vec3 = structured_fv::Vec3<T>;

void expect_near(const Vec3<double>& a, const Vec3<double>& b)
{
  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(a[i], b[i], 1e-13);
}
}


TEST(Math, ComplexIntOps)
{
  Complex a(2, 3);
  Int b = 4;
  UInt c = 5;
  Real br = b;
  Real cr = c;

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

TEST(Math, ArrayPlusArray)
{
  Vec3<double> a = {1, 2, 3};
  Vec3<double> b = {4, 5, 6};
  expect_near(a+b, {5, 7, 9});
}

TEST(Math, ArrayPlusEqualArray)
{
  Vec3<double> a = {1, 2, 3};
  Vec3<int> b = {4, 5, 6};
  a += b;
  expect_near(a, {5, 7, 9});
}

TEST(Math, ArrayPlusScalar)
{
  Vec3<double> a = {1, 2, 3};
  expect_near(a + 2, {3, 4, 5});
  expect_near(2 + a, {3, 4, 5});
}

TEST(Math, ArrayPlusEqualScalar)
{
  Vec3<double> a = {1, 2, 3};
  a += 2;
  expect_near(a, {3, 4, 5});
}

TEST(Math, ArrayMinusArray)
{
  Vec3<double> a = {1, 2, 3};
  Vec3<double> b = {2, 5, 7};
  expect_near(a-b, {-1, -3, -4});
}

TEST(Math, ArrayMinusEqualArray)
{
  Vec3<double> a = {1, 2, 3};
  Vec3<double> b = {2, 5, 7};
  a -= b;
  expect_near(a, {-1, -3, -4});
}

TEST(Math, ArrayMinusScalar)
{
  Vec3<double> a = {1, 2, 3};
  expect_near(a - 2, {-1, 0, 1});
  expect_near(2 - a, {1, 0, -1});
}

TEST(Math, ArrayUnaryMinus)
{
  Vec3<double> a = {1, 2, 3};
  Vec3<double> b = -a;
  expect_near(b, {-1, -2, -3});
}

TEST(Math, ArrayMinusEqualScalar)
{
  Vec3<double> a = {1, 2, 3};
  a -= 2;
  expect_near(a, {-1, 0, 1});
}

TEST(Math, ArrayMultiplyArray)
{
  Vec3<double> a = {1, 2, 3};
  Vec3<double> b = {4, 5, 6};
  expect_near(a*b, {4, 10, 18});
}

TEST(Math, ArrayMultiplyEqualArray)
{
  Vec3<double> a = {1, 2, 3};
  Vec3<double> b = {4, 5, 6};
  a *= b;
  expect_near(a, {4, 10, 18});
}

TEST(Math, ArrayMultiplyScalar)
{
  Vec3<double> a = {1, 2, 3};
  expect_near(a * 2, {2, 4, 6});
  expect_near(2 * a, {2, 4, 6});
}

TEST(Math, ArrayMultiplyEqualScalar)
{
  Vec3<double> a = {1, 2, 3};
  a *= 2;
  expect_near(a, {2, 4, 6});
}

TEST(Math, ArrayDivideArray)
{
  Vec3<double> a = {1, 2, 3};
  Vec3<double> b = {4, 5, 6};
  expect_near(a/b, {1.0/4, 2.0/5, 3.0/6});
}

TEST(Math, ArrayDivideEqualArray)
{
  Vec3<double> a = {1, 2, 3};
  Vec3<double> b = {4, 5, 6};
  a /= b;
  expect_near(a, {1.0/4, 2.0/5, 3.0/6});
}

TEST(Math, ArrayDivideScalar)
{
  Vec3<double> a = {1, 2, 3};
  expect_near(a / 2, {1.0/2, 1, 3.0/2});
  expect_near(2 / a, {2, 1, 2.0/3});
}

TEST(Math, ArrayDivideEqualScalar)
{
  Vec3<double> a = {1, 2, 3};
  a /= 2;
  expect_near(a, {1.0/2, 1, 3.0/2});
}

TEST(Math, DotProduct)
{
  Vec3<double> a = {1, 2, 3};
  Vec3<double> b = {4, 5, 6};
  EXPECT_NEAR(dot(a, b), 4 + 10 + 18, 1e-13);
}

TEST(Math, CrossProduct)
{
  Vec3<double> a = {1, 2, 3};
  Vec3<double> b = {4, 5, 6};
  expect_near(cross(a, b), {-3, 6, -3});
}

TEST(Math, RadianConversion)
{
  EXPECT_NEAR(degreesToRadians(30), PI/6, 1e-13);
  EXPECT_NEAR(radiansToDegrees(PI/6), 30, 1e-13);
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