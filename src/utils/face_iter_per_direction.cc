#include "face_iter_per_direction.h"
#include <iostream>

namespace structured_fv {


std::ostream& operator<<(std::ostream& os, const XDirTag& tag)
{
  os << "XDirTag";
  return os;
}

std::ostream& operator<<(std::ostream& os, const YDirTag& tag)
{
  os << "YDirTag";
  return os;
}

}