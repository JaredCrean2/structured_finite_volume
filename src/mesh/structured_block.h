#ifndef STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_BLOCK_H
#define STRUCTURED_FINITE_VOLUME_MESH_STRUCTURED_BLOCK_H

#include "utils/project_defs.h"
#include "utils/range.h"

namespace structured_fv {
namespace mesh {


struct MeshBlockSpec
{
  UInt numCellsX;
  UInt numCellsY;

  UInt rotation;  // 0 = no rotation, 1 = rotate 90 degrees counter-clockwise, 2 = 180 degrees, 3 = 270 degrees
                  // this allows testing the transform between block interfaces
  std::function<std::array<double, 2>(UInt i, UInt j, UInt numVertsX, UInt numVertsY)> m_coordFunc;
};

enum class BlockType
{
  Regular,
  GhostBC
};

class StructuredBlock
{
  public:
    BlockType getBlockType() const;

    UInt getBlockId() const;  // block IDs are global

    Range2D getOwnedVerts() const;  // these are the verts of owned cells
                                    // slightly misleading, because we dont define ownership of
                                    // verts

    Range2D getOwnedCells() const;

    const Kokkos::View<Real**, HostMemorySpace>& getOwnedVertCoords() const;

    std::array<int, 2> getOffsetIntoBlock();  // owned idx + offset = block idx

    std::array<int, 2> getBlockSize() const;

  private:

};


class StructuredBlockInterface
{
  public:

  private:
    UInt m_left_block_id;
    UInt m_right_block_id;
    std::array<Int, 2> m_transform;  // for each direction in the left block, gives the
                                     // corresponding direction in the right block
                                     // ex. if the +x direction in the left block is
                                     // the negative y direction in right block, then
                                     // m_transform[0] = -1;
}


}
}

