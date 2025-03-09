#include "structured_mesh.h"
#include "mesh/adjacent_block_indexer.h"
#include "mesh/block_spec.h"
#include "mesh/structured_block.h"
#include "mesh/structured_block_interface.h"
#include "utils/project_defs.h"
#include "utils/math.h"

namespace structured_fv {
namespace mesh {



StructuredMesh::StructuredMesh(const MeshSpec& spec)
{
  assert(spec.blocks.extent(0) > 0 && spec.blocks.extent(1) > 0);
  Range2D block_range(0, spec.blocks.extent(0), 0, spec.blocks.extent(1));
  Range x_iface_range(0, spec.blocks.extent(0)-1);
  Range y_iface_range(0, spec.blocks.extent(1)-1);

  UInt block_id = 0;
  for (UInt j : block_range.getYRange())
    for (UInt i : block_range.getXRange())
      m_blocks.emplace_back(spec.blocks(i, j), block_id++, BlockType::Regular);

  m_block_counts[0] = block_range.getXRange().size() * block_range.getYRange().size();


  for (UInt j : block_range.getYRange())
    for (UInt i : x_iface_range)
    {
      UInt left_block_id = i + j * spec.blocks.extent(0);
      UInt right_block_id = left_block_id + 1;

      const StructuredBlock& blockL = m_blocks[left_block_id];
      const StructuredBlock& blockR = m_blocks[right_block_id];

      NeighborDirection dirL = rotate(NeighborDirection::East, spec.blocks(i, j).rotation);
      Range rangeL = to_int(dirL) % 2 == 0 ? blockL.getOwnedCells().getXRange() : 
                                             blockL.getOwnedCells().getYRange();
      std::array<Int, 2> transformL = getTransform(spec.blocks(i, j).rotation, spec.blocks(i+1, j).rotation);

      NeighborDirection dirR = getNeighborImage(dirL, transformL);
      Range rangeR = to_int(dirR) % 2 == 0 ?  m_blocks[right_block_id].getOwnedCells().getXRange() : 
                                              m_blocks[right_block_id].getOwnedCells().getYRange();

      m_block_interfaces.emplace_back(blockL, dirL, rangeL, transformL, blockR, rangeR);
      m_block_iface_counts[0]++;
    }

  for (UInt j : y_iface_range)
    for (UInt i : block_range.getXRange())
    {
      UInt bottom_block_id = i + j * spec.blocks.extent(0);
      UInt top_block_id = bottom_block_id + spec.blocks.extent(0);

      const StructuredBlock& blockB = m_blocks[bottom_block_id];
      const StructuredBlock& blockT = m_blocks[top_block_id];

      NeighborDirection dirB = rotate(NeighborDirection::North, spec.blocks(i, j).rotation);
      Range rangeB = to_int(dirB) % 2 == 0 ? blockB.getOwnedCells().getXRange() : 
                                             blockB.getOwnedCells().getYRange();      
      std::array<Int, 2> transformB = getTransform(spec.blocks(i, j).rotation, spec.blocks(i, j+1).rotation);

      NeighborDirection dirT = getNeighborImage(dirB, transformB);
      Range rangeT = to_int(dirT) % 2 == 0 ?  m_blocks[top_block_id].getOwnedCells().getXRange() : 
                                              m_blocks[top_block_id].getOwnedCells().getYRange();

      m_block_interfaces.emplace_back(blockB, dirB, rangeB, transformB, blockT, rangeT);
      m_block_iface_counts[0]++;
    }

  createBCGhosts(spec);

  m_block_interface_connectivity.resize(getNumBlocks(), {-1, -1, -1, -1});
  for (UInt i=0; i < getNumBlockInterfaces(); ++i)
  {
    //TODO: this does not work for T junctions
    const StructuredBlockInterface& iface = getBlockInterface(i);
    m_block_interface_connectivity[iface.getBlockIdL()][to_int(iface.getNeighborDirectionL())] = i;
    m_block_interface_connectivity[iface.getBlockIdR()][to_int(iface.getNeighborDirectionR())] = i;
  }
}


UInt StructuredMesh::getNumBlocks() const { return m_block_counts[0] + m_block_counts[1]; }

UInt StructuredMesh::getNumRegularBlocks() const { return m_block_counts[0]; }

UInt StructuredMesh::getNumGhostBCBlocks() const { return m_block_counts[1]; }

const StructuredBlock& StructuredMesh::getBlock(UInt block) const
{
  return m_blocks.at(block);
}

const StructuredBlock& StructuredMesh::getRegularBlock(UInt block) const
{
  return m_blocks.at(block);
}

const StructuredBlock& StructuredMesh::getGhostBCBlock(UInt bc_block) const
{
  return m_blocks.at(getNumRegularBlocks() + bc_block);
}    

UInt StructuredMesh::getNumBlockInterfaces() const { return m_block_iface_counts[0] + m_block_iface_counts[1]; }

UInt StructuredMesh::getNumRegularBlockInterfaces() const { return m_block_iface_counts[0]; }

UInt StructuredMesh::getNumGhostBCBlockInterfaces() const { return m_block_iface_counts[1]; }

const StructuredBlockInterface& StructuredMesh::getBlockInterface(UInt iface) const
{
  return m_block_interfaces.at(iface);
}

const StructuredBlockInterface& StructuredMesh::getRegularBlockInterface(UInt iface) const
{
  return m_block_interfaces.at(iface);
}

const StructuredBlockInterface& StructuredMesh::getGhostBCBlockInterface(UInt iface) const
{
  return m_block_interfaces.at(getNumRegularBlockInterfaces() + iface);
}

UInt StructuredMesh::getNumBCs() const { return m_bc_block_ranges.size(); }

// gives the range of ghost blocks that are part of the given BC
Range StructuredMesh::getBCBlockRange(UInt bc) const { return m_bc_block_ranges.at(bc); }

// gives the range of ghost block interfaces that are part of the BC
// (ie. one of the blocks in a ghost BC block)
Range StructuredMesh::getBCInterfaceRange(UInt bc) const { return m_bc_iface_ranges.at(bc); }

std::array<Int, 4> StructuredMesh::getBlockInterfaces(UInt block) const { return m_block_interface_connectivity.at(block); }


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

void StructuredMesh::createBCGhost(const MeshBlockSpec& spec, UInt regular_block_id, NeighborDirection domain_boundary)
{
  UInt ghost_block_id = m_blocks.size();
  UInt regular_block_rotation = spec.rotation;
  // the ghost block always has the same coordinate system as the regular block
  std::array<Int, 2> transform = {1, 2};

  NeighborDirection regular_block_dir = rotate(domain_boundary, regular_block_rotation);
  NeighborDirection ghost_block_dir = getNeighborImage(regular_block_dir, transform);

  auto ghost_coords = computeGhostBCCoords(m_blocks[regular_block_id], ghost_block_dir, 2);

  m_blocks.emplace_back(ghost_coords, ghost_block_id);
  m_block_counts[1]++;

  const StructuredBlock& regular_block = m_blocks[regular_block_id];
  const StructuredBlock& ghost_block = m_blocks[ghost_block_id];

  Range regular_block_range, ghost_block_range;
  if (to_int(regular_block_dir) % 2 == 0)
  {
    regular_block_range = regular_block.getOwnedCells().getXRange();
    ghost_block_range   = ghost_block.getOwnedCells().getXRange();
  } else
  {
    regular_block_range = regular_block.getOwnedCells().getYRange();
    ghost_block_range   = ghost_block.getOwnedCells().getYRange();
  }

  m_block_interfaces.emplace_back(regular_block, regular_block_dir, regular_block_range, transform,
                                  ghost_block, ghost_block_range);
  m_block_iface_counts[1]++;
}
  




}
}
