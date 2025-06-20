#include "gtest/gtest.h"
#include "utils/vec.h"

namespace {
using Vec3 = structured_fv::FixedVec<int, 3>;

void expect_eq(const Vec3& arr1, const Vec3& arr2)
{
  EXPECT_EQ(arr1[0], arr2[0]);
  EXPECT_EQ(arr1[1], arr2[1]);
  EXPECT_EQ(arr1[2], arr2[2]);  
}

}

TEST(Vec, Constructor)
{
  Vec3 v1{1, 2, 3};
  EXPECT_EQ(v1[0], 1);
  EXPECT_EQ(v1[1], 2);
  EXPECT_EQ(v1[2], 3);  

  Vec3 v2 = {1, 2, 3};
  expect_eq(v2, {1, 2, 3});

  Vec3 v3(v2);
  expect_eq(v3, {1, 2, 3});

  Vec3 v4(std::move(v2));
  expect_eq(v4, {1, 2, 3});
}


TEST(Vec, Assignment)
{
  Vec3 v1{1, 2, 3};

  Vec3 v2 = v1;
  expect_eq(v2, {1, 2, 3});

  Vec3 v3 = std::move(v1);
  expect_eq(v3, {1, 2, 3});
}

TEST(Vec, ConstAssignment)
{
  Vec3 v1{1, 2, 3};
  const Vec3& v2_const = v1;

  Vec3 v2 = v2_const;
  expect_eq(v2, {1, 2, 3});

  Vec3 v3 = std::move(v2_const);
  expect_eq(v3, {1, 2, 3});
}

TEST(Vec, Indexing)
{
  Vec3 v1{1, 2, 3};
  expect_eq(v1, {1, 2, 3});

  v1[1] = 4;
  expect_eq(v1, {1, 4, 3});

  v1.at(1) = 5;
  expect_eq(v1, {1, 5, 3});

  EXPECT_ANY_THROW(v1.at(3));
}

TEST(Vec, ConstIndexing)
{
  Vec3 v1{1, 2, 3};
  expect_eq(v1, {1, 2, 3});

  const Vec3& v2_const = v1;
  EXPECT_EQ(v2_const[1], 2);
  EXPECT_EQ(v2_const.at(1), 2);
  EXPECT_ANY_THROW(v2_const.at(3));
}

TEST(Vec, Front)
{
  Vec3 v1{1, 2, 3};
  const Vec3& v2_const = v1;

  EXPECT_EQ(v1.front(), 1);
  EXPECT_EQ(v2_const.front(), 1);

  v1.front() = 4;
  EXPECT_EQ(v1.front(), 4);
  EXPECT_EQ(v2_const.front(), 4);
}

TEST(Vec, Back)
{
  Vec3 v1{1, 2, 3};
  const Vec3& v2_const = v1;

  EXPECT_EQ(v1.back(), 3);
  EXPECT_EQ(v2_const.back(), 3);

  v1.back() = 4;
  EXPECT_EQ(v1.back(), 4);
  EXPECT_EQ(v2_const.back(), 4);
}

TEST(Vec, Data)
{
  Vec3 v1{1, 2, 3};
  const Vec3& v2_const = v1;

  EXPECT_EQ(v1.data(), &(v1[0]));
  EXPECT_EQ(v2_const.data(), &(v2_const[0]));
}

TEST(Vec, Iterators)
{
  std::vector<int> copy1, copy2, copy3;

  Vec3 v1{1, 2, 3};
  const Vec3& v2_const = v1;
  for (auto val : v1)
    copy1.push_back(val);

  for (auto val : v2_const)
    copy2.push_back(val);

  for (auto it = v2_const.cbegin(); it != v2_const.cend(); ++it)
    copy3.push_back(*it);

  EXPECT_EQ(copy1.size(), 3U);
  EXPECT_EQ(copy2.size(), 3U);
  EXPECT_EQ(copy3.size(), 3U);
  for (int i=0; i < 3; ++i)
  {
    EXPECT_EQ(copy1[i], v1[i]);
    EXPECT_EQ(copy2[i], v1[i]);
    EXPECT_EQ(copy3[i], v1[i]);
  }
}

TEST(Vec, Empty)
{
  Vec3 v1{1, 2, 3};

  EXPECT_FALSE(v1.empty());
  structured_fv::FixedVec<int, 0> v2{};
  EXPECT_TRUE(v2.empty());
}

TEST(Vec, Size)
{
  Vec3 v1{1, 2, 3};
  EXPECT_EQ(v1.size(), 3U);
  EXPECT_EQ(v1.max_size(), 3U);
}

TEST(Vec, Fill)
{
  Vec3 v1{1, 2, 3};
  v1.fill(4);
  expect_eq(v1, {4, 4, 4});
}

TEST(Vec, Swap)
{
  Vec3 v1{1, 2, 3}, v2{4, 5, 6};
  auto it1 = v1.begin();
  auto it2 = v2.begin();
  EXPECT_EQ(*it1, 1);
  EXPECT_EQ(*it2, 4);

  v1.swap(v2);
  expect_eq(v1, {4, 5, 6});
  expect_eq(v2, {1, 2, 3});

  EXPECT_EQ(*it1, 4);
  EXPECT_EQ(*it2, 1);
}

TEST(Vec, Equality)
{
  Vec3 v1{1, 2, 3}, v2{4, 5, 6}, v3{1, 2, 3};
  EXPECT_FALSE(v1 == v2);
  EXPECT_TRUE(v1 == v3);

  EXPECT_TRUE(v1 != v2);
  EXPECT_FALSE(v1 != v3);
}

TEST(Vec, LessThan)
{
  Vec3 v1{1, 2, 3}, v2{0, 1, 2}, v3{1, 2, 3}, v4{2, 3, 4};
  EXPECT_FALSE(v1 < v2);
  EXPECT_FALSE(v1 < v3);
  EXPECT_TRUE(v1 < v4);
}

TEST(Vec, LessThanOrEqul)
{
  Vec3 v1{1, 2, 3}, v2{0, 1, 2}, v3{1, 2, 3}, v4{2, 3, 4};
  EXPECT_FALSE(v1 < v2);
  EXPECT_TRUE(v1 <= v3);
  EXPECT_TRUE(v1 < v4);
}

TEST(Vec, GreaterThan)
{
  Vec3 v1{1, 2, 3}, v2{0, 1, 2}, v3{1, 2, 3}, v4{2, 3, 4};
  EXPECT_TRUE(v1 > v2);
  EXPECT_FALSE(v1 > v3);
  EXPECT_FALSE(v1 > v4);
}

TEST(Vec, GreaterThanOrEqual)
{
  Vec3 v1{1, 2, 3}, v2{0, 1, 2}, v3{1, 2, 3}, v4{2, 3, 4};
  EXPECT_TRUE(v1 >= v2);
  EXPECT_TRUE(v1 >= v3);
  EXPECT_FALSE(v1 >= v4);
}

TEST(Vec, Get)
{
  Vec3 v1{1, 2, 3};
  const Vec3& v2_const = v1;

  EXPECT_EQ(structured_fv::get<0>(v1), 1);
  EXPECT_EQ(structured_fv::get<1>(v1), 2);
  EXPECT_EQ(structured_fv::get<2>(v1), 3);


  EXPECT_EQ(structured_fv::get<0>(v2_const), 1);
  EXPECT_EQ(structured_fv::get<1>(v2_const), 2);
  EXPECT_EQ(structured_fv::get<2>(v2_const), 3);

  structured_fv::get<0>(v1) = 4;
  EXPECT_EQ(structured_fv::get<0>(v1), 4);
}

TEST(Vec, OutputOperator)
{
  std::stringstream ss;
  Vec3 v1{1, 2, 3};
  ss << v1;
  EXPECT_EQ(ss.str(), "1, 2, 3");
}