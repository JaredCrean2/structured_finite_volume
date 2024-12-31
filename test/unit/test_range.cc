#include "gtest/gtest.h"
#include "utils/range.h"
#include <vector>

using namespace structured_fv;

TEST(RangeIter, Dereference)
{
  RangeIter iter(2);
  EXPECT_EQ(*iter, 2);
}

TEST(RangeIter, PrefixIncrement)
{
  RangeIter iter(2);
  EXPECT_EQ(*(++iter), 3);
  EXPECT_EQ(*iter, 3);
}

TEST(RangeIter, PostfixIncrement)
{
  RangeIter iter(2);
  EXPECT_EQ(*(iter++), 2);
  EXPECT_EQ(*iter, 3);
}

TEST(RangeIter, PrefixDecrement)
{
  RangeIter iter(2);
  EXPECT_EQ(*(--iter), 1);
  EXPECT_EQ(*iter, 1);
}

TEST(RangeIter, PostfixDecrement)
{
  RangeIter iter(2);
  EXPECT_EQ(*(iter--), 2);
  EXPECT_EQ(*iter, 1);
}

TEST(RangeIter, PlusEqual)
{
  RangeIter iter(2);
  EXPECT_EQ(*(iter += 2), 4);
  EXPECT_EQ(*iter, 4);
}

TEST(RangeIter, MinusEqual)
{
  RangeIter iter(2);
  EXPECT_EQ(*(iter -= 2), 0);
  EXPECT_EQ(*iter, 0);
}

TEST(RangeIter, Plus)
{
  RangeIter iter(2);
  EXPECT_EQ(*(iter + 2), 4);
  EXPECT_EQ(*iter, 2);
}

TEST(RangeIter, Minus)
{
  RangeIter iter(2);
  EXPECT_EQ(*(iter - 2), 0);
  EXPECT_EQ(*iter, 2);
}

TEST(RangeIter, MinusIterator)
{
  RangeIter iter(2);
  EXPECT_EQ(iter - RangeIter(2), 0);
  EXPECT_EQ(*iter, 2);
}

TEST(RangeIter, Equality)
{
  EXPECT_TRUE(RangeIter(2) == RangeIter(2));
  EXPECT_FALSE(RangeIter(2) != RangeIter(2));

  EXPECT_FALSE(RangeIter(2) == RangeIter(3));
  EXPECT_TRUE(RangeIter(2) != RangeIter(3));
}

TEST(RangeIter, RelationalOperators)
{
  EXPECT_TRUE(RangeIter(2)  < RangeIter(3));
  EXPECT_TRUE(RangeIter(2)  <= RangeIter(3));
  EXPECT_FALSE(RangeIter(2) > RangeIter(3));
  EXPECT_FALSE(RangeIter(2) >= RangeIter(3));

  EXPECT_FALSE(RangeIter(3)  < RangeIter(2));
  EXPECT_FALSE(RangeIter(3)  <= RangeIter(2));
  EXPECT_TRUE(RangeIter(3) > RangeIter(2));
  EXPECT_TRUE(RangeIter(3) >= RangeIter(2));

  EXPECT_FALSE(RangeIter(2)  < RangeIter(2));
  EXPECT_TRUE(RangeIter(2)  <= RangeIter(2));
  EXPECT_FALSE(RangeIter(2) > RangeIter(2));
  EXPECT_TRUE(RangeIter(2) >= RangeIter(2));    
}

TEST(Range, CallOperator)
{
  Range range(2, 4);
  EXPECT_EQ(range(0), 2);
  EXPECT_EQ(range(1), 3);
  EXPECT_EQ(range(2), 4);
}

TEST(Range, RangeFor)
{
  Range range(2, 4);
  std::vector<UInt> vals;

  for (UInt i : range)
    vals.push_back(i);


  std::vector<UInt> vals_expected = {2, 3, 4};
  EXPECT_EQ(vals, vals_expected);
}

TEST(Range2D, RangeFor)
{
  Range2D range(2, 4, 5, 6);
  std::vector<std::pair<UInt, UInt>> vals;

  for (UInt i : range.getXRange())
    for (UInt j : range.getYRange())
      vals.emplace_back(i, j);

  std::vector<std::pair<UInt, UInt>> vals_expected = { {2, 5}, {2, 6},
                                                       {3, 5}, {3, 6},
                                                       {4, 5}, {4, 6}};
  EXPECT_EQ(vals, vals_expected);
}