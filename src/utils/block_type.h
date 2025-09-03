#ifndef STRUCTURED_FINITE_VOLUME_BLOCK_TYPE_H
#define STRUCTURED_FINITE_VOLUME_BLOCK_TYPE_H

#include "utils/neighbor_direction.h"
#include "vec.h"
#include <iosfwd>

namespace structured_fv {

enum class BlockType
{
  Regular,
  GhostBC,
  Invalid
};

std::ostream& operator<<(std::ostream& os, BlockType type);


// stores the type of all the neighboring blocks.
// Note: this still works for meshes that have T intersections,
//       because even if a block has two neighboring blocks in the
//       same direction, those blocks still have to be the same type
class NeighborBlockTypes
{
  public:
    NeighborBlockTypes() :
      m_neighbor_block_types{}
    {}

    NeighborBlockTypes(const Vec4<BlockType>& neighbor_block_types) :
      m_neighbor_block_types(neighbor_block_types)
    {}

    BlockType getBlockType(NeighborDirection dir) const
    {
      return m_neighbor_block_types[to_int(dir)];
    }

  private:
    Vec4<BlockType> m_neighbor_block_types;
};

}

#endif