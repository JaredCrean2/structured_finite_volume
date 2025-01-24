#include "gtest/gtest.h"
#include "utils/math.h"

namespace {

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

TEST(Math, ArrayPlusScalar)
{
  Vec3<double> a = {1, 2, 3};
  expect_near(a + 2, {3, 4, 5});
  expect_near(2 + a, {3, 4, 5});

}

TEST(Math, ArrayMinusArray)
{
  Vec3<double> a = {1, 2, 3};
  Vec3<double> b = {2, 5, 7};
  expect_near(a-b, {-1, -3, -4});
}

TEST(Math, ArrayMinusScalar)
{
  Vec3<double> a = {1, 2, 3};
  expect_near(a - 2, {-1, 0, 1});
  expect_near(2 - a, {1, 0, -1});

}

TEST(Math, ArrayMultiplyArray)
{
  Vec3<double> a = {1, 2, 3};
  Vec3<double> b = {4, 5, 6};
  expect_near(a*b, {4, 10, 18});
}

TEST(Math, ArrayMultiplyScalar)
{
  Vec3<double> a = {1, 2, 3};
  expect_near(a * 2, {2, 4, 6});
  expect_near(2 * a, {2, 4, 6});
}

TEST(Math, ArrayDivideArray)
{
  Vec3<double> a = {1, 2, 3};
  Vec3<double> b = {4, 5, 6};
  expect_near(a/b, {1.0/4, 2.0/5, 3.0/6});
}

TEST(Math, ArrayDivideScalar)
{
  Vec3<double> a = {1, 2, 3};
  expect_near(a / 2, {1.0/2, 1, 3.0/2});
  expect_near(2 / a, {2, 1, 2.0/3});
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