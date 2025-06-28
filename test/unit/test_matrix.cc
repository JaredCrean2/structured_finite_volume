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

TEST(Matrix, AssignmentOperator)
{
  Matrix<Int, 3, 4> mat;
  mat = {1,  2,  3,  4,
         5,  6,  7,  8,
         9, 10, 11, 12};
  
  Int val = 1;
  for (UInt i=0; i < 3; ++i)
    for (UInt j=0; j < 4; ++j)
      EXPECT_EQ(mat(i, j), val++);
}

TEST(Matrix, DefaultTemplateArgument)
{
  Matrix<Int, 3> mat;
  EXPECT_EQ(mat.extent0(), 3);
  EXPECT_EQ(mat.extent1(), 3);
  EXPECT_EQ(mat.size().first, 3);
  EXPECT_EQ(mat.size().second, 3);
}

TEST(Matrix, MatVec)
{
  Matrix<Int, 3, 2> mat({1, 2, 3, 4, 5, 6});
  FixedVec<Int, 2> x{2, 3};
  FixedVec<Int, 3> b = mat * x;

  EXPECT_EQ(b[0], 1*2 + 2*3);
  EXPECT_EQ(b[1], 3*2 + 4*3);
  EXPECT_EQ(b[2], 5*2 + 6*3);
}

TEST(Matrix, RowIndex)
{
  Matrix<Int, 3, 2> mat({1, 2, 3, 4, 5, 6});

  auto x = mat(Row{1});
  EXPECT_EQ(x.size(), 2);
  EXPECT_EQ(x[0], 3);
  EXPECT_EQ(x[1], 4);
}

TEST(Matrix, ColIndex)
{
  Matrix<Int, 3, 2> mat({1, 2, 3, 4, 5, 6});

  auto x = mat(Column{1});
  EXPECT_EQ(x.size(), 3);
  EXPECT_EQ(x[0], 2);
  EXPECT_EQ(x[1], 4);
  EXPECT_EQ(x[2], 6);
}

TEST(Matrix, RowMajor)
{
  Matrix<Int, 3, 2> mat({1, 2, 3, 4, 5, 6});
  Int* data = mat.getData();
  EXPECT_EQ(&(mat(0, 0)), data);
  EXPECT_EQ(&(mat(0, 1)), data+1);
  EXPECT_EQ(&(mat(1, 0)), data+2);
  EXPECT_EQ(&(mat(1, 1)), data+3);
}