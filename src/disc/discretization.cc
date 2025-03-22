#include "discretization.h"
#include "disc/disc_block.h"
#include "disc/face_field.h"
#include "disc/vert_field.h"
#include "disc/elem_field.h"
#include "utils/face_iter_per_direction.h"
#include "utils/math.h"

#include <iostream>  //TODO: DEBUGGING

namespace structured_fv {
namespace disc {

StructuredDisc::StructuredDisc(std::shared_ptr<mesh::StructuredMesh> mesh, UInt num_ghost_cells) :
  m_mesh(mesh)
  //m_num_ghost_cells(num_ghost_cells)
{
  for (UInt i=0; i < mesh->getNumBlocks(); ++i)
  {
    m_blocks.emplace_back(*mesh, mesh->getBlock(i), num_ghost_cells);
  }

  for (UInt i=0; i < mesh->getNumBlockInterfaces(); ++i)
  {
    const mesh::StructuredBlockInterface& mesh_iface = mesh->getBlockInterface(i);
    m_ifaces.emplace_back(m_blocks[mesh_iface.getBlockIdL()], m_blocks[mesh_iface.getBlockIdR()], mesh_iface);
  }

  m_coordField = createCoordField();
  m_normalField = createFaceNormalField();
  m_invCellVolumeField = createInvCellVolumeField();
}


VertexFieldPtr<Real> StructuredDisc::createCoordField()
{
  auto coordField = std::make_shared<VertexField<Real>>(*this, 2);
  for (UInt i=0; i < getNumBlocks(); ++i)
  {
    const StructuredBlock& block = getBlock(i);
    const mesh::StructuredBlock& mesh_block = m_mesh->getBlock(i);
    const auto& meshVertCoords = mesh_block.getOwnedVertCoords();
    auto& blockVertCoords = coordField->getData(i);
    for (UInt i : block.getOwnedVerts().getXRange())
      for (UInt j : block.getOwnedVerts().getYRange())
      {
        auto [imesh, jmesh] = block.blockVertToMeshVert(i, j);
        for (UInt d=0; d < 2; ++d)
        {
          blockVertCoords(i, j, d) = meshVertCoords(imesh, jmesh, d);
        }
      }
  }

  coordField->updateGhostValues();

  return coordField;
}

ElementFieldPtr<Real> StructuredDisc::createInvCellVolumeField()
{
  auto inv_cell_volume_field = std::make_shared<ElementField<Real>>(*this, 1);
  inv_cell_volume_field->set(0);

  for (UInt block_id=0; block_id < getNumBlocks(); ++block_id)
  {
    const StructuredBlock& block = getBlock(block_id);
    const auto& coords = getCoordField()->getData(block_id);
    const auto& inv_volumes = inv_cell_volume_field->getData(block_id);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
      {
        std::array<std::array<Real, 2>, 4> cell_coords;
        for (UInt d=0; d < 2; ++d)
        {
          cell_coords[0][d] = coords(i, j, d);
          cell_coords[1][d] = coords(i+1, j, d);
          cell_coords[2][d] = coords(i+1, j+1, d);
          cell_coords[3][d] = coords(i, j+1, d);

          Real volume = computeQuadArea(cell_coords);
          if (volume < 0)
          {
            throw std::runtime_error(std::string("found negative volume for block ") + std::to_string(block_id) +
                                     ", cell " + std::to_string(i) + ", " + std::to_string(j));
          }

          inv_volumes(i, j, 0) = 1.0/volume;
        }
      }
  }

  inv_cell_volume_field->updateGhostValues();

  return inv_cell_volume_field;
}

FaceFieldPtr<Real> StructuredDisc::createFaceNormalField()
{
  auto normal_field = std::make_shared<FaceField<Real>>(*this, 2);
  normal_field->set(0);

  constexpr  Vec3<Real> unit_z{0, 0, 1};
  for (UInt block_id=0; block_id < getNumBlocks(); ++block_id)
  {
    const StructuredBlock& block = getBlock(block_id);
    const auto& coords = getCoordField()->getData(block_id);
    const auto& east_normals = normal_field->getData(block_id, NeighborDirection::East);
    const auto& north_normals = normal_field->getData(block_id, NeighborDirection::North);

    auto compute_normals = [&](UInt i, UInt j)
    {
      Vec3<Real> bottom_edge{coords(i+1, j, 0) - coords(i, j, 0),
                              coords(i+1, j, 1) - coords(i, j, 1),
                              0};

      Vec3<Real> right_edge{coords(i+1, j+1, 0) - coords(i+1, j, 0),
                            coords(i+1, j+1, 1) - coords(i+1, j, 1),
                            0};

      Vec3<Real> top_edge{coords(i, j+1, 0) - coords(i+1, j+1, 0),
                          coords(i, j+1, 1) - coords(i+1, j+1, 1),
                          0};

      Vec3<Real> left_edge{coords(i, j, 0) - coords(i, j+1, 0),
                           coords(i, j, 1) - coords(i, j+1, 1),
                           0};
      
      // normal is outward from the left/bottom cell
      Vec3<Real> bottom_normal = cross(unit_z, bottom_edge);
      Vec3<Real> right_normal  = cross(right_edge, unit_z);
      Vec3<Real> top_normal    = cross(top_edge, unit_z);
      Vec3<Real> left_normal   = cross(unit_z, left_edge);

      for (UInt d=0; d < 2; ++d)
      {
        north_normals(i,   j,   d) = bottom_normal[d];
        east_normals( i+1, j,   d) = right_normal[d];
        north_normals(i,   j+1, d) = top_normal[d];
        east_normals( i,   j,   d) = left_normal[d];
      }
    };

    for (UInt i : block.getOwnedAndGhostCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        compute_normals(i, j);

    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedAndGhostCells().getYRange())
        compute_normals(i, j);     
  }

  return normal_field;
}



}
}