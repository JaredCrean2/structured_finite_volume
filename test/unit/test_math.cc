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
  EXPECT_EQ(smoothAbs(5, delta), std::abs(5));
  EXPECT_EQ(smoothAbs(-5, delta), std::abs(-5));
  EXPECT_NEAR(smoothAbs(0.05, delta), 0.05*0.05/delta, 1e-13);
  EXPECT_NEAR(smoothAbs(-0.05, delta), 0.05*0.05/delta, 1e-13);
}

TEST(Math, SmoothAbsDeriv)
{
  double delta = 0.1;
  EXPECT_EQ(smoothAbsDeriv(5, delta), 1); 
  EXPECT_EQ(smoothAbsDeriv(-5, delta), -1);
  EXPECT_NEAR(smoothAbsDeriv(0.05, delta), 2*0.05/delta, 1e-13);
  EXPECT_NEAR(smoothAbsDeriv(-0.05, delta), -2*0.05/delta, 1e-13);
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