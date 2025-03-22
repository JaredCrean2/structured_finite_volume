#ifndef STRUCTURED_FINITE_VOLUME_DISC_VERT_FIELD_H
#define STRUCTURED_FINITE_VOLUME_DISC_VERT_FIELD_H

#include "discretization.h"
#include "utils/traits.h"

namespace structured_fv {
namespace disc {


// defines a field that stores n values in each element of each block (including ghosts)
template <typename T>
class VertexField
{
  public:
    using FieldData = Kokkos::View<T***, HostMemorySpace>;
    using ConstFieldData = Kokkos::View<const T***, HostMemorySpace>;

    VertexField(const StructuredDisc& disc, Int nvals_per_element) :
      m_disc(disc),
      m_nvals_per_element(nvals_per_element)
    {
      for (UInt i=0; i < disc.getNumBlocks(); ++i)
      {
        const StructuredBlock& block = disc.getBlock(i);
        auto dimensions = block.getVertDimensions();
        m_data.emplace_back("field_data", dimensions[0], dimensions[1], nvals_per_element);
      }
    }

    VertexField(const VertexField& other) = delete;

    VertexField& operator=(const VertexField& other) = delete;

    UInt getNumBlocks() const { return m_data.size(); }

    FieldData& getData(UInt block) { return m_data[block]; }

    const ConstFieldData& getData(UInt block) const { return m_data[block]; }

    T& operator()(UInt block, UInt i, UInt j, UInt v) { return m_data[block](i, j, v); }      

    const T& operator()(UInt block, UInt i, UInt j, UInt v) const { return m_data[block](i, j, v); }

    void set(const T& val);

    // Func is a callable object (Real x, Real y) -> std::array<T, num_vals_per_element>
    // return type can be anything of the correct length that supports operator[]
    template <typename Func, IsFuncXY_t<Func> = true>
    void set(Func func);

    void updateGhostValues();

  private:
    std::vector<FieldData> m_data;
    const StructuredDisc& m_disc;
    const UInt m_nvals_per_element;
};


template <typename T>
void VertexField<T>::set(const T& val)
{
  for (UInt b=0; b < m_disc.getNumBlocks(); ++b)
  {
    const StructuredBlock& block = m_disc.getBlock(b);
    auto data = getData(b);
    for (UInt i : block.getOwnedAndGhostVerts().getXRange())
      for (UInt j : block.getOwnedAndGhostVerts().getYRange())
        for (UInt d=0; d < m_nvals_per_element; ++d)
          data(i, j, d) = val;
  }
}

// Func is a callable object (Real x, Real y) -> std::array<T, num_vals_per_element>
// return type can be anything of the correct length that supports operator[]
template <typename T>
template <typename Func, IsFuncXY_t<Func>>
void VertexField<T>::set(Func func)
{
  for (UInt b=0; b < m_disc.getNumBlocks(); ++b)
  {
    const StructuredBlock& block = m_disc.getBlock(b);
    auto data = getData(b);
    auto coords = m_disc.getCoordField()->getData(b);
    for (UInt i : block.getOwnedAndGhostCells().getXRange())
      for (UInt j : block.getOwnedAndGhostCells().getYRange())
      {
        const auto& vals = func(coords(i, j, 0), coords(i, j, 1));
        for (UInt d=0; d < m_nvals_per_element; ++d)
          data(i, j, d) = vals[d];
      }
  }      
}


template <typename T>
void VertexField<T>::updateGhostValues()
{
  for (UInt i=0; i < m_disc.getNumBlockInterfaces(); ++i)
  {
    const StructuredBlockInterface& iface = m_disc.getBlockInterface(i);
    const StructuredBlock& blockL = m_disc.getBlock(iface.getBlockIdL());
    int num_ghost_cells = blockL.getNumGhostCellsPerDirection()[to_int(iface.getNeighborDirectionL())];
    const auto& indexerL = iface.getAdjBlockVertIndexerL();
    const auto& indexerR = iface.getAdjBlockVertIndexerR();
    auto fieldL = getData(iface.getBlockIdL());
    auto fieldR = getData(iface.getBlockIdR());
    
    for (UInt i : iface.getOwnedBoundaryVertsL().getXRange())
      for (UInt j : iface.getOwnedBoundaryVertsL().getYRange())
        for (int v=1; v <= num_ghost_cells; ++v)
        {
          auto [iprime, jprime] = computeIndices(iface.getNeighborDirectionL(), v, i, j);
          auto [ineighbor, jneighbor] = indexerL(iprime, jprime);
          auto [ineighbor2, jneighbor2] = computeIndices(iface.getNeighborDirectionR(), -1, ineighbor, jneighbor);

          for (UInt d=0; d < fieldL.extent(2); ++d)
          {
            fieldL(iprime, jprime, d) = fieldR(ineighbor2, jneighbor2, d);
          } 
        }

    //TODO: it would have better temporal locality to merge this
    //      into the above loops
    for (UInt i : iface.getOwnedBoundaryVertsR().getXRange())
      for (UInt j : iface.getOwnedBoundaryVertsR().getYRange())
        for (int v=1; v <= num_ghost_cells; ++v)
        {
          auto [iprime, jprime] = computeIndices(iface.getNeighborDirectionR(), v, i, j);
          auto [ineighbor, jneighbor] = indexerR(iprime, jprime);
          auto [ineighbor2, jneighbor2] = computeIndices(iface.getNeighborDirectionL(), -1, ineighbor, jneighbor);

          for (UInt d=0; d < fieldL.extent(2); ++d)
          {
            fieldR(iprime, jprime, d) = fieldL(ineighbor2, jneighbor2, d);
          }      
        }        
  }
}

}
}

#endif