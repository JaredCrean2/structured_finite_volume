#include "gtest/gtest.h"
#include "utils/matrix.h"

using namespace structured_fv;

TEST(Matrix, DefaultConstructor)
{
  Matrix<Int, 3, 4> mat;
  EXPECT_EQ(mat.extent0(), 3);
  EXPECT_EQ(mat.extent1(), 4);
  EXPECT_EQ(mat.size().first, 3);
  EXPECT_EQ(mat.size().second, 4);
  for (UInt i=0; i < 3; ++i)
    for (UInt j=0; j < 4; ++j)
      EXPECT_EQ(mat(i, j), 0);

  const Matrix<Int, 3, 4>& mat_const = mat;
  EXPECT_EQ(mat_const.extent0(), 3);
  EXPECT_EQ(mat_const.extent1(), 4);
  EXPECT_EQ(mat_const.size().first, 3);
  EXPECT_EQ(mat_const.size().second, 4);  
  for (UInt i=0; i < 3; ++i)
    for (UInt j=0; j < 4; ++j)
      EXPECT_EQ(mat_const(i, j), 0);  
}

TEST(Matrix, ValueConstructor)
{
  Matrix<Int, 3, 4> mat({1,  2,  3,  4,
                         5,  6,  7,  8,
                         9, 10, 11, 12});
  
  Int val = 1;
  for (UInt i=0; i < 3; ++i)
    for (UInt j=0; j < 4; ++j)
      EXPECT_EQ(mat(i, j), val++);

  const Matrix<Int, 3, 4>& mat_const = mat;
  val = 1;
  for (UInt i=0; i < 3; ++i)
    for (UInt j=0; j < 4; ++j)
      EXPECT_EQ(mat_const(i, j), val++);
}

TEST(Matrix, DefaultTemplateArgument)
{
  Matrix<Int, 3> mat;
  EXPECT_EQ(mat.extent0(), 3);
  EXPECT_EQ(mat.extent1(), 3);
  EXPECT_EQ(mat.size().first, 3);
  EXPECT_EQ(mat.size().second, 3);
}