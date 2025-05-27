#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_RECONSTRUCTION_ENUM_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_RECONSTRUCTION_ENUM_H

#include <iosfwd>

namespace structured_fv {
namespace euler {

// a reconstruction is the combination of two choices: a slope limiter
// and which variables to apply the slope limiter to
enum class Reconstruction
{
  Conservative,
  Primitive,
};

std::string get_name(Reconstruction recon);

std::ostream& operator<<(std::ostream& os, Reconstruction recon);

}
}

#endif