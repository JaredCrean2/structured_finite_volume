#include "range.h"
#include <iostream>

namespace structured_fv {

std::ostream& operator<<(std::ostream& os, const Range& range)
{
  os << "[" << *(range.begin()) << ", " << *(range.end()) << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, const Range2D& range)
{
  os << range.getXRange() << " x " << range.getYRange();
  return os;
}

}