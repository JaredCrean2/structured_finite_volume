#ifndef STRUCTURED_FINITE_VOLUME_UTILS_RANGE_H
#define STRUCTURED_FINITE_VOLUME_UTILS_RANGE_H

#include "project_defs.h"
#include <iterator>
#include <limits>
#include <iostream>

namespace structured_fv {

class RangeIter
{
  public:
    using value_type = UInt;
    using difference_type = std::ptrdiff_t;
    //using pointer = value_type*;
    using reference = value_type&;
    using const_reference = const value_type&;
    using iterator_catagory = std::random_access_iterator_tag;

    constexpr RangeIter() :
      m_idx(std::numeric_limits<UInt>::max())
    {}

    constexpr RangeIter(UInt idx) :
      m_idx(idx)
    {}

    // prefix
    constexpr RangeIter& operator++()
    {
      ++m_idx;
      return *this;
    }

    constexpr RangeIter& operator--()
    {
      --m_idx;
      return *this;
    }

    // postfix
    constexpr RangeIter operator++(int)
    {
      RangeIter tmp = *this;
      m_idx++;
      return tmp;
    }

    constexpr RangeIter operator--(int)
    {
      RangeIter tmp = *this;
      m_idx--;
      return tmp;
    }

    constexpr reference operator*()
    {
      return m_idx;
    }

    constexpr const_reference operator*() const
    {
      return m_idx;
    }

    constexpr RangeIter& operator+=(const difference_type rhs)
    {
      m_idx += rhs;
      return *this;
    }

    constexpr RangeIter& operator-=(const difference_type rhs)
    {
      m_idx -= rhs;
      return *this;
    }

    constexpr RangeIter operator+(const difference_type rhs)
    {
      RangeIter tmp = *this;
      return tmp += rhs;
    }

    constexpr RangeIter operator-(const difference_type rhs)
    {
      RangeIter tmp = *this;
      return tmp -= rhs;
    }

    constexpr difference_type operator-(const RangeIter& rhs)
    {
      return difference_type(m_idx) - difference_type(rhs.m_idx);
    }

    //TODO: should have operator[n] here, but I'm not sure how to implement it
    //      because it returns a reference

    constexpr bool operator<(const RangeIter& rhs) const { return m_idx < rhs.m_idx; }

    constexpr bool operator<=(const RangeIter& rhs) const { return m_idx <= rhs.m_idx; }

    constexpr bool operator>(const RangeIter& rhs) const { return m_idx > rhs.m_idx; }

    constexpr bool operator>=(const RangeIter& rhs) const { return m_idx >= rhs.m_idx; }

    constexpr bool operator==(const RangeIter& rhs) const { return m_idx == rhs.m_idx; }

    constexpr bool operator!=(const RangeIter& rhs) const { return m_idx != rhs.m_idx; }

  private:
    UInt m_idx;
};


class Range
{
  public:
    using iterator = RangeIter;

    constexpr Range(UInt start, UInt past_the_end) :
      m_start(start),
      m_past_the_end(std::max(start, past_the_end))
    {}

    constexpr UInt operator()(UInt i) const 
    {
      assert (i < size());
      return m_start + i;
    }

    constexpr RangeIter begin() const { return RangeIter(m_start); }

    constexpr RangeIter end() const { return RangeIter(m_past_the_end); }

    constexpr UInt size() const { return m_past_the_end - m_start; }

    constexpr bool operator==(const Range& rhs) const
    {
      return m_start == rhs.m_start && m_past_the_end == rhs.m_past_the_end;
    }

    constexpr bool operator!=(const Range& rhs) const
    {
      return !(*this == rhs);
    }    

  private:
    UInt m_start;
    UInt m_past_the_end;
};

inline std::ostream& operator<<(std::ostream& os, const Range& range)
{
  os << "[" << *(range.begin()) << ", " << *(range.end()) << ")";
  return os;
}

constexpr bool in(const Range& range, UInt val)
{
  return val >= *range.begin() && val < *range.end();

}


class Range2D
{
  public:
    constexpr Range2D(UInt xstart=0, UInt x_past_the_end=0, UInt ystart=0, UInt y_past_the_end=0) :
      m_xrange(xstart, x_past_the_end),
      m_yrange(ystart, y_past_the_end)
    {}

    constexpr const Range& getXRange() const { return m_xrange; }

    constexpr const Range& getYRange() const { return m_yrange; }

    constexpr bool operator==(const Range2D& rhs) const
    {
      return m_xrange == rhs.m_xrange &&
             m_yrange == rhs.m_yrange;
    }

    constexpr bool operator!=(const Range2D& rhs) const
    {
      return !(*this == rhs);
    }

  private:
    Range m_xrange;
    Range m_yrange;
};

inline std::ostream& operator<<(std::ostream& os, const Range2D& range)
{
  os << range.getXRange() << " x " << range.getYRange();
  return os;
}

constexpr bool in(const Range2D range, UInt x, UInt y)
{
  return in(range.getXRange(), x) && in(range.getYRange(), y);
}

}

#endif