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

    RangeIter() :
      m_idx(std::numeric_limits<UInt>::max())
    {}

    RangeIter(UInt idx) :
      m_idx(idx)
    {}

    // prefix
    RangeIter& operator++()
    {
      ++m_idx;
      return *this;
    }

    RangeIter& operator--()
    {
      --m_idx;
      return *this;
    }

    // postfix
    RangeIter operator++(int)
    {
      RangeIter tmp = *this;
      m_idx++;
      return tmp;
    }

    RangeIter operator--(int)
    {
      RangeIter tmp = *this;
      m_idx--;
      return tmp;
    }

    reference operator*()
    {
      return m_idx;
    }

    const_reference operator*() const
    {
      return m_idx;
    }

    RangeIter& operator+=(const difference_type rhs)
    {
      m_idx += rhs;
      return *this;
    }

    RangeIter& operator-=(const difference_type rhs)
    {
      m_idx -= rhs;
      return *this;
    }

    RangeIter operator+(const difference_type rhs)
    {
      RangeIter tmp = *this;
      return tmp += rhs;
    }

    RangeIter operator-(const difference_type rhs)
    {
      RangeIter tmp = *this;
      return tmp -= rhs;
    }

    difference_type operator-(const RangeIter& rhs)
    {
      return difference_type(m_idx) - difference_type(rhs.m_idx);
    }

    //TODO: should have operator[n] here, but I'm not sure how to implement it
    //      because it returns a reference

    bool operator<(const RangeIter& rhs) const { return m_idx < rhs.m_idx; }

    bool operator<=(const RangeIter& rhs) const { return m_idx <= rhs.m_idx; }

    bool operator>(const RangeIter& rhs) const { return m_idx > rhs.m_idx; }

    bool operator>=(const RangeIter& rhs) const { return m_idx >= rhs.m_idx; }

    bool operator==(const RangeIter& rhs) const { return m_idx == rhs.m_idx; }

    bool operator!=(const RangeIter& rhs) const { return m_idx != rhs.m_idx; }

  private:
    UInt m_idx;
};


class Range
{
  public:
    using iterator = RangeIter;

    Range(UInt start, UInt past_the_end) :
      m_start(start),
      m_past_the_end(std::max(start, past_the_end))
    {}

    UInt operator()(UInt i) const 
    {
      assert (i < size());
      return m_start + i;
    }

    RangeIter begin() const { return RangeIter(m_start); }

    RangeIter end() const { return RangeIter(m_past_the_end); }

    UInt size() const { return m_past_the_end - m_start; }

    bool operator==(const Range& rhs) const
    {
      return m_start == rhs.m_start && m_past_the_end == rhs.m_past_the_end;
    }

    bool operator!=(const Range& rhs) const
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


class Range2D
{
  public:
    Range2D(UInt xstart=0, UInt x_past_the_end=0, UInt ystart=0, UInt y_past_the_end=0) :
      m_xrange(xstart, x_past_the_end),
      m_yrange(ystart, y_past_the_end)
    {}

    const Range& getXRange() const { return m_xrange; }

    const Range& getYRange() const { return m_yrange; }

    bool operator==(const Range2D& rhs) const
    {
      return m_xrange == rhs.m_xrange &&
             m_yrange == rhs.m_yrange;
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

}

#endif