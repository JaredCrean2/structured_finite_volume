#include "reconstruction_enum.h"

#include <iostream>

namespace structured_fv {
namespace euler {

std::string get_name(Reconstruction recon)
{
  switch (recon)
  {
    case Reconstruction::Conservative: { return "Conservative"; }
    case Reconstruction::Primitive:    { return "Primitive"; }
    default:
      throw std::runtime_error("Unhandled Reconstruction enum");
  }
}

std::ostream& operator<<(std::ostream& os, Reconstruction recon)
{
  return os << get_name(recon);
  return os;
}

}
}