#ifndef STRUCTURED_FINITE_VOLUME_TEST_UNIT_ROTATION_UTILS_H
#define STRUCTURED_FINITE_VOLUME_TEST_UNIT_ROTATION_UTILS_H

#include "utils/project_defs.h"
#include "mesh/block_spec.h"
#include <stdexcept>

#include <iostream>

namespace test_utils {

using structured_fv::UInt;
using structured_fv::Real;
// utilities for rotating the coordinate system of blocks while
// keeping the shape of the block the same

constexpr std::pair<UInt, UInt> compute_block_size(UInt nx, UInt ny, UInt rotation)
{
  bool is_reversed = rotation % 2 == 1;
  UInt nx_prime = is_reversed ? ny : nx;
  UInt ny_prime = is_reversed ? nx : ny;

  return {nx_prime, ny_prime};
}

// given the indices (i, j) of a cell for a rotation=0 block, and the size of that block,
// returns the indices of the same cell for a rotated block
constexpr std::pair<UInt, UInt> compute_indices(UInt i, UInt j, UInt nx, UInt ny, UInt rotation)
{
  switch (rotation)
  {
    case 0 : { return {i, j}; }
    case 1 : { return {j, nx - i - 1}; }
    case 2 : { return {nx - i - 1, ny - j - 1}; }
    case 3 : { return {ny - j - 1, i}; }
    default :
      throw std::runtime_error("invalid rotation value");
  }
}

//TODO: it would be better if BlockSpec used a coord_func that is in the given rotation, and
//      then we handle the rotation here.
/*
inline structured_fv::mesh::Fxy rotate_mapping(structured_fv::mesh::Fxy func, UInt rotation)
{
  // x, y are in rotation=0 coordinate system
  // xprime, yprime are in rotated coordinate system

  auto f1 = [=](Real xprime, Real yprime)
  { 
    std::cout << "xprime, yprime = " << xprime << ", " << yprime << std::endl;
    Real x = 1 - yprime;
    Real y = xprime;
    //std::cout << "x, y = " << x<< ", " << y << std::endl;    
    //return func(x, y);
    return func(xprime, yprime);

  };

  auto f2 = [=](Real xprime, Real yprime)
  { 
    Real x = 1 - xprime;
    Real y = 1 - yprime;
    return func(x, y);
  };

  auto f3 = [=](Real xprime, Real yprime)
  { 
    Real x = yprime;
    Real y = 1 - xprime;
    return func(x, y);
  };    

  switch (rotation)
  {
    case 0: { return func; }
    case 1: { return f1; }
    case 2: { return f2; }
    case 3: { return f3; }
    default:
      throw std::runtime_error("invalid rotation");
  }
}
*/

// rotates the coordinate system of the block while keeping it the same
// shape in the global coordinate system.
// For example, if a block is 3 x 4 in the global coordinate system (rotation=0),
// then the block size is 4 x 3 in rotation=1 (90 degrees ccw rotation).
inline structured_fv::mesh::MeshBlockSpec rotate_block(const structured_fv::mesh::MeshBlockSpec& block, UInt rotation)
{
  auto [nx, ny] = compute_block_size(block.num_cells_x, block.num_cells_y, rotation);
  //auto func = rotate_mapping(block.coord_func, rotation);
  return structured_fv::mesh::MeshBlockSpec(nx, ny, rotation, block.coord_func);
}


}


#endif