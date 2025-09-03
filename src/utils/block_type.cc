#include "block_type.h"
#include <iostream>

namespace structured_fv {

std::ostream& operator<<(std::ostream& os, BlockType type)
{
  switch (type)
  {
    case BlockType::Regular: { os << "Regular"; break; }
    case BlockType::GhostBC: { os << "GhostBC"; break; }
    case BlockType::Invalid: { os << "Invalid"; break; }
    default:
      throw std::runtime_error("unhandled BlockType enum value");
  }

  return os;
}

}