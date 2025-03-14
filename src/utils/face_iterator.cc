#include "face_iterator.h"
#include <iostream>

namespace structured_fv {

std::ostream& operator<<(std::ostream& os, const FaceId& faceid)
{
  os << "cellL: (" << faceid.cell_i_left << ", " << faceid.cell_j_left << "), dir " << faceid.dirL << ", "
     << "cellR: (" << faceid.cell_i_right << ", " << faceid.cell_j_right << "), dir " << faceid.dirR;

  return os;
}
}