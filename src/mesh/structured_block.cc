#include "structured_block.h"

namespace structured_fv {
namespace mesh {

namespace {
// given an index (i, j) in the standard orientation (ie rotation=0), compute
// the location of the corresponding point under the rotation
inline std::pair<double, double> applyRotation(UInt rotation, double x, double y)
{
  switch (rotation)
  {
    case 0: { return {x,   y}; }
    case 1: { return {1 - y,  x}; }
    case 2: { return {1 - x, 1 - y}; }
    case 3: { return {y, 1-x}; }
    default: { throw std::runtime_error("invalid rotation"); }
  }
}
}

StructuredBlock::StructuredBlock(const MeshBlockSpec& spec, UInt block_id, BlockType block_type) :
  m_block_id(block_id),
  m_block_type(block_type),
  m_owned_cell_range(0,  spec.num_cells_x, 0, spec.num_cells_y),
  m_owned_vert_range(0, spec.num_cells_x+1, 0, spec.num_cells_y+1),
  m_owned_vert_coords("block_coords", spec.num_cells_x+1, spec.num_cells_y+1),
  m_offset_into_block{0, 0},
  m_all_block_size{spec.num_cells_x, spec.num_cells_y}
{
  UInt dimx = spec.num_cells_x + 1;
  UInt dimy = spec.num_cells_y + 1;

  for (UInt i : Range(0, dimx))
    for (UInt j : Range(0, dimy))
    {
      auto [x, y] = applyRotation(spec.rotation, double(i) / (dimx - 1), double(j) / (dimy - 1));
      auto [xprime, yprime] = spec.coord_func(x, y);

      m_owned_vert_coords(i, j, 0) = xprime;
      m_owned_vert_coords(i, j, 1) = yprime;
    }
}

StructuredBlock::StructuredBlock(Kokkos::View<Real**[2], HostMemorySpace> ghost_coords, UInt block_id) :
  m_block_id(block_id),
  m_block_type(BlockType::GhostBC),
  m_owned_cell_range(0, ghost_coords.extent(0)-1, 0, ghost_coords.extent(1)-1),
  m_owned_vert_range(0, ghost_coords.extent(0),    0, ghost_coords.extent(1)),
  m_owned_vert_coords(ghost_coords),
  m_offset_into_block{0, 0},
  m_all_block_size{static_cast<UInt>(ghost_coords.extent(0)-1), static_cast<UInt>(ghost_coords.extent(1)-1)}
{}

UInt StructuredBlock::getBlockId() const { return m_block_id;};  // block IDs are global

BlockType StructuredBlock::getBlockType() const { return m_block_type; };

// these are the verts of owned cells
// slightly misleading, because we dont define ownership of verts
Range2D StructuredBlock::getOwnedVerts() const { return m_owned_vert_range; }

Range2D StructuredBlock::getOwnedCells() const { return m_owned_cell_range; };

const StructuredBlock::CoordsHostView& StructuredBlock::getOwnedVertCoords() const { return m_owned_vert_coords; };

std::array<UInt, 2> StructuredBlock::getOffsetIntoBlock() const { return m_offset_into_block; };  // owned idx + offset = block idx

// returns size of the entire block (including non-owned cells)
std::array<UInt, 2> StructuredBlock::getAllBlockSize() const { return m_all_block_size; };

UInt getNumOwnedCells(const StructuredBlock& block, NeighborDirection dir)
{
  return block.getOwnedCells().getDimensions()[to_int(dir) % 2];
}

}
}