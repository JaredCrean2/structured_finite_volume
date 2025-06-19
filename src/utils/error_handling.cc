#include <stdexcept>
#include "utils/error_handling.h"

namespace structured_fv {

void assertAlways(const bool condition, const std::string& msg)
{
  if (!condition)
    throw std::runtime_error(msg);
}

}
