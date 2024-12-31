#ifndef STRUCTURED_FINITE_VOLUME_UTILS_RANGE_H
#define STRUCTURED_FINITE_VOLUME_UTILS_RANGE_H

#include "project_defs.h"
#include <iterator>
#include <limits>

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

    Range(UInt start, UInt end) :
      m_start(start),
      m_end(end)
    {}

    UInt operator()(UInt i) const { return m_start + i; }

    RangeIter begin() const { return RangeIter(m_start); }

    RangeIter end() const { return RangeIter(m_end + 1); }

  private:
    UInt m_start;
    UInt m_end;
};

class Range2D
{
  public:
    Range2D(UInt xstart, UInt xend, UInt ystart, UInt yend) :
      m_xrange(xstart, xend),
      m_yrange(ystart, yend)
    {}

    const Range& getXRange() const { return m_xrange; }

    const Range& getYRange() const { return m_yrange; }

  private:
    Range m_xrange;
    Range m_yrange;
};
}

#endif