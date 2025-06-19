#ifndef STRUCTURED_FINITE_VOLUME_UTILS_ERROR_HANDLING_H
#define STRUCTURED_FINITE_VOLUME_UTILS_ERROR_HANDLING_H

#include <string>

namespace structured_fv {

void assertAlways(const bool condition, const std::string& msg);

}

#endif