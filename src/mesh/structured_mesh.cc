#include "structured_mesh.h"
#include "mesh/adjacent_block_indexer.h"
#include "mesh/block_spec.h"
#include "mesh/structured_block.h"
#include "mesh/structured_block_interface.h"
#include "utils/project_defs.h"
#include "utils/math.h"

namespace structured_fv {
namespace mesh {

UInt getConstantIndexAlongBoundary(const StructuredBlock& block, NeighborDirection dir)
{
  switch (dir)
  {
    case NeighborDirection::North: { return *(block.getOwnedCells().getYRange().end())-1; }
    case NeighborDirection::East:  { return *(block.getOwnedCells().getXRange().end())-1; }
    case NeighborDirection::South: { return *(block.getOwnedCells().getXRange().begin()); }
    case NeighborDirection::West:  { return *(block.getOwnedCells().getYRange().begin()); }
  }
}

Range getVariableIndexAlongBoundary(const StructuredBlock& block, NeighborDirection dir)
{
  if (static_cast<int>(dir) % 2 == 0)
  {
    return block.getOwnedCells().getXRange();
  } else
  {
    return block.getOwnedCells().getYRange();
  }
}

std::array<UInt, 2> getMinCellOnBoundary(const StructuredBlock& right_block, NeighborDirection left_block_dir, int right_block_rotation)
{
  UInt nx = right_block.getOwnedCells().getXRange().size();
  UInt ny = right_block.getOwnedCells().getYRange().size();
  std::array<std::array<UInt, 2>, 4> min_cells;
  if (left_block_dir == NeighborDirection::North || left_block_dir == NeighborDirection::East)
  {
    min_cells[0] = {0, 0};
    min_cells[1] = {0, ny-1};
    min_cells[2] = {nx-1, ny-1};
    min_cells[3] = {nx-1, 0};
  } else if (left_block_dir == NeighborDirection::South)
  {
    min_cells[0] = {0, ny-1};
    min_cells[1] = {nx-1, ny-1};
    min_cells[2] = {nx-1, 0};
    min_cells[3] = {0, 0};
  } else if (left_block_dir == NeighborDirection::West)
  {
    min_cells[0] = {nx-1, 0};
    min_cells[1] = {0, 0};
    min_cells[2] = {0, ny-1};
    min_cells[3] = {nx-1, ny-1};    
  }

  return min_cells[right_block_rotation];
}

std::array<Int, 2> getTransform(Int rotation_left, Int rotation_right)
{
  Int delta_r = (rotation_right - rotation_left + 4) % 4;
  std::array<std::array<Int, 2>, 4> transforms;
  transforms[0] = {1, 2};
  transforms[1] = {-2, 1};
  transforms[2] = {-1, -2};
  transforms[3] = {2, -1};

  return transforms[delta_r];
}

StructuredMesh::StructuredMesh(const MeshSpec& spec)
{
  assert(spec.blocks.extent(0) > 0 && spec.blocks.extent(1) > 0);
  Range2D block_range(0, spec.blocks.extent(0), 0, spec.blocks.extent(1));
  Range x_iface_range(1, spec.blocks.extent(0)-1);
  Range y_iface_range(1, spec.blocks.extent(1)-1);

  UInt block_id = 0;
  for (UInt j : block_range.getYRange())
    for (UInt i : block_range.getXRange())
      m_blocks.emplace_back(spec.blocks(i, j), block_id++, BlockType::Regular);

  m_block_counts[0] = block_range.getXRange().size() * block_range.getYRange().size();


  for (UInt j : block_range.getYRange())
    for (UInt i : x_iface_range)
    {
      UInt left_block_id = (i-1) + j * spec.blocks.extent(0);
      UInt right_block_id = left_block_id + 1;

      NeighborDirection dir = rotate(NeighborDirection::East, spec.blocks(i-1, j).rotation);
      UInt left_block_constant_index = getConstantIndexAlongBoundary(m_blocks[left_block_id], dir);
      Range left_block_variable_index = getVariableIndexAlongBoundary(m_blocks[left_block_id], dir);
      std::array<Int, 2> transform = getTransform(spec.blocks(i-1, j).rotation, spec.blocks(i, j).rotation);

      std::array<UInt, 2> right_block_min_cell = getMinCellOnBoundary(m_blocks[right_block_id], dir, spec.blocks(i, j).rotation);

      m_block_interfaces.emplace_back(left_block_id, left_block_constant_index, 
                                      left_block_variable_index, dir, transform,
                                      right_block_id, right_block_min_cell);

      m_block_iface_counts[0]++;
    }

  for (UInt j : y_iface_range)
    for (UInt i : block_range.getXRange())
    {
      UInt bottom_block_id = i + (j-1) * spec.blocks.extent(0);
      UInt top_block_id = bottom_block_id + spec.blocks.extent(0);

      NeighborDirection dir = rotate(NeighborDirection::North, spec.blocks(i, j-1).rotation);
      UInt bottom_block_constant_index = getConstantIndexAlongBoundary(m_blocks[bottom_block_id], dir);
      Range bottom_block_variable_index = getVariableIndexAlongBoundary(m_blocks[bottom_block_id], dir);
      std::array<Int, 2> transform = getTransform(spec.blocks(i, j-1).rotation, spec.blocks(i, j).rotation);

      std::array<UInt, 2> top_block_min_cell = getMinCellOnBoundary(m_blocks[top_block_id], dir, spec.blocks(i, j).rotation);

      m_block_interfaces.emplace_back(bottom_block_id, bottom_block_constant_index, 
                                      bottom_block_variable_index, dir, transform,
                                      top_block_id, top_block_min_cell);
      m_block_iface_counts[0]++;
    }

  createBCGhosts(spec);

  m_block_interface_connectivity.resize(getNumBlocks(), {-1, -1, -1, -1});
  for (UInt i=0; i < getNumBlockInterfaces(); ++i)
  {
    const StructuredBlockInterface& iface = getBlockInterface(i);
    m_block_interface_connectivity[iface.getLeftBlockId()][to_int(iface.getNeighborDirection())] = i;
    m_block_interface_connectivity[iface.getRightBlockId()][to_int(iface.getOtherBlockNeighborDirection())] = i;
  }
}

// dir is the direction of the interface from the ghost block perspective
Kokkos::View<double**[2], HostMemorySpace> computeGhostBCCoords(const StructuredBlock& block, NeighborDirection dir, UInt num_ghost_cells)
{
  //TODO: assert block size

  Range xrange = static_cast<int>(dir) % 2 == 0 ? block.getOwnedVerts().getXRange() : Range(0, num_ghost_cells+1);
  Range yrange = static_cast<int>(dir) % 2 == 0 ? Range(0, num_ghost_cells+1) : block.getOwnedVerts().getYRange();
  Kokkos::View<double**[2], HostMemorySpace> ghost_coords("ghost_bc_coords", xrange.size(), yrange.size());
  auto block_coords = block.getOwnedVertCoords();
  UInt block_nx = block.getOwnedVerts().getXRange().size();
  UInt block_ny = block.getOwnedVerts().getYRange().size();

  for (UInt i : xrange)
    for (UInt j : yrange)
    {
      std::array<double, 2> x1;
      std::array<double, 2> x0;

      switch(dir)
      {
        case NeighborDirection::North: 
        {
          UInt j2 = num_ghost_cells - j;
          x0 = {block_coords(i, 0, 0), block_coords(i, 0, 1)};
          x1 = {block_coords(i, 1, 0), block_coords(i, 1, 1)};
          auto normal = x0 - x1;

          ghost_coords(i, j, 0) = j2*normal[0] + x0[0];          
          ghost_coords(i, j, 1) = j2*normal[1] + x0[1];
          break;
        }

        case NeighborDirection::East: 
        {
          UInt i2 = num_ghost_cells - i;
          x0 = {block_coords(0, j, 0), block_coords(0, j, 1)};
          x1 = {block_coords(1, j, 0), block_coords(1, j, 1)};
          auto normal = x0 - x1;

          ghost_coords(i, j, 0) = i2*normal[0] + x0[0];
          ghost_coords(i, j, 1) = i2*normal[1] + x0[1];
          break;
        }

        case NeighborDirection::South: 
        {
          x0 = {block_coords(i, block_ny-1, 0), block_coords(i, block_ny-1, 1)};
          x1 = {block_coords(i, block_ny-2, 0), block_coords(i, block_ny-2, 1)};
          auto normal = x0 - x1;


          ghost_coords(i, j, 0) = j*normal[0] + x0[0];
          ghost_coords(i, j, 1) = j*normal[1] + x0[1];
          break;
        }

        case NeighborDirection::West: 
        {          
          x0 = {block_coords(block_nx-1, j, 0), block_coords(block_nx-1, j, 1)};
          x1 = {block_coords(block_nx-2, j, 0), block_coords(block_nx-2, j, 1)};
          auto normal = x0 - x1;

          ghost_coords(i, j, 0) = i*normal[0] + x0[0];   
          ghost_coords(i, j, 1) = i*normal[1] + x0[1];
          break;
        }
      }
    }

  return ghost_coords;
}

void StructuredMesh::createBCGhosts(const MeshSpec& spec)
{
  std::array<NeighborDirection, 4> directions = {NeighborDirection::North,
                                                 NeighborDirection::East,
                                                 NeighborDirection::South,
                                                 NeighborDirection::West};
  std::array<Range2D, 4> ranges = {Range2D(0, spec.blocks.extent(0), spec.blocks.extent(1)-1, spec.blocks.extent(1)),
                                   Range2D(spec.blocks.extent(0)-1, spec.blocks.extent(0), 0, spec.blocks.extent(1)),
                                   Range2D(0, spec.blocks.extent(0), 0, 1),
                                   Range2D(0, 1, 0, spec.blocks.extent(1))};

  for (int dir=0; dir < 4; ++dir)
  {
    UInt ghost_block_id_start = m_blocks.size();
    UInt ghost_iface_start = m_block_interfaces.size();    
    for (UInt i : ranges[dir].getXRange())
      for (UInt j : ranges[dir].getYRange())
      {
        UInt regular_block_id = i + j * spec.blocks.extent(0);
        const MeshBlockSpec& block_spec = spec.blocks(i, j);
        createBCGhost(block_spec, regular_block_id, directions[dir]);
      }

    UInt ghost_block_id_end = m_blocks.size();
    UInt ghost_iface_end = m_block_interfaces.size();
    m_bc_block_ranges.push_back(Range(ghost_block_id_start, ghost_block_id_end));
    m_bc_iface_ranges.push_back(Range(ghost_iface_start, ghost_iface_end));
  }
}

/*
void StructuredMesh::createBCGhosts(const MeshSpec& spec)
{
  Range block_x_range(0, spec.blocks.extent(0));
  Range block_y_range(0, spec.blocks.extent(1));

  // North
  UInt ghost_block_id_start = m_blocks.size();
  UInt ghost_iface_start = m_block_interfaces.size();
  for (UInt i : block_x_range)
  {
    UInt j = spec.blocks.extent(1)-1;
    UInt regular_block_id = i + j * spec.blocks.extent(0);
    const MeshBlockSpec& block_spec = spec.blocks(i, j);
    createBCGhost(block_spec, regular_block_id, NeighborDirection::North);
  }
  UInt ghost_block_id_end = m_blocks.size();
  UInt ghost_iface_end = m_block_interfaces.size();
  m_bc_block_ranges.push_back(Range(ghost_block_id_start, ghost_block_id_end));
  m_bc_iface_ranges.push_back(Range(ghost_iface_start, ghost_iface_end));

  // East
  ghost_block_id_start = m_blocks.size();
  ghost_iface_start = m_block_interfaces.size();  
  for (UInt j : block_y_range)
  {
    UInt i = spec.blocks.extent(0) - 1;
    UInt regular_block_id = i + j*spec.blocks.extent(0);
    const MeshBlockSpec& block_spec = spec.blocks(i, j);
    createBCGhost(block_spec, regular_block_id, NeighborDirection::East);
  }  
  ghost_block_id_end = m_blocks.size();
  ghost_iface_end = m_block_interfaces.size();
  m_bc_block_ranges.push_back(Range(ghost_block_id_start, ghost_block_id_end));
  m_bc_iface_ranges.push_back(Range(ghost_iface_start, ghost_iface_end));

  // South
  ghost_block_id_start = m_blocks.size();
  ghost_iface_start = m_block_interfaces.size();   
  for (UInt i : block_x_range)
  {
    UInt j = 0 ;
    UInt regular_block_id = i;
    const MeshBlockSpec& block_spec = spec.blocks(i, j);
    createBCGhost(block_spec, regular_block_id, NeighborDirection::South);
  }
  ghost_block_id_end = m_blocks.size();
  ghost_iface_end = m_block_interfaces.size();
  m_bc_block_ranges.push_back(Range(ghost_block_id_start, ghost_block_id_end));
  m_bc_iface_ranges.push_back(Range(ghost_iface_start, ghost_iface_end));

  // West
  ghost_block_id_start = m_blocks.size();
  ghost_iface_start = m_block_interfaces.size();   
  for (UInt j : block_y_range)
  {
    UInt i = 0;
    UInt regular_block_id = i + j*spec.blocks.extent(0);
    const MeshBlockSpec& block_spec = spec.blocks(i, j);
    createBCGhost(block_spec, regular_block_id, NeighborDirection::West);
  }
  ghost_block_id_end = m_blocks.size();
  ghost_iface_end = m_block_interfaces.size();
  m_bc_block_ranges.push_back(Range(ghost_block_id_start, ghost_block_id_end));
  m_bc_iface_ranges.push_back(Range(ghost_iface_start, ghost_iface_end));       
}
*/

void StructuredMesh::createBCGhost(const MeshBlockSpec& spec, UInt regular_block_id, NeighborDirection domain_boundary)
{
  UInt ghost_block_id = m_blocks.size();
  UInt regular_block_rotation = spec.rotation;
  // the ghost block always has the same coordinate system as the regular block
  UInt ghost_block_rotation = regular_block_rotation;
  std::array<Int, 2> transform = {1, 2};

  NeighborDirection regular_block_dir = rotate(domain_boundary, regular_block_rotation);
  NeighborDirection ghost_block_dir = static_cast<NeighborDirection>((static_cast<int>(domain_boundary) + 2) % 4);
  ghost_block_dir = rotate(ghost_block_dir, regular_block_rotation);


  auto ghost_coords = computeGhostBCCoords(m_blocks[regular_block_id], ghost_block_dir, 2);

  m_blocks.emplace_back(ghost_coords, ghost_block_id);
  const StructuredBlock& ghost_block = m_blocks[ghost_block_id];
  m_block_counts[1]++;

  const StructuredBlock& regular_block = m_blocks[regular_block_id];
  UInt left_block_constant_index = getConstantIndexAlongBoundary(regular_block,  regular_block_dir);
  Range left_block_variable_index = getVariableIndexAlongBoundary(regular_block, regular_block_dir);
  std::array<UInt, 2> min_cell_on_boundary = getMinCellOnBoundary(ghost_block, regular_block_dir, ghost_block_rotation);

  m_block_interfaces.emplace_back(regular_block_id, left_block_constant_index, left_block_variable_index, regular_block_dir,
                                  transform, ghost_block_id, min_cell_on_boundary);
  m_block_iface_counts[1]++;
}
  




}
}
