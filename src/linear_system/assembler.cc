#include "assembler.h"
#include <iostream>

namespace structured_fv {
namespace linear_system {

std::ostream& operator<<(std::ostream& os, const Indices& ind)
{
  os << ind.i << ", " << ind.j;
  return os;
}

}
}