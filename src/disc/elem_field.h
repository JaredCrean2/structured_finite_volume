#ifndef STRUCTURED_FINITE_VOLUME_DISC_ELEM_FIELD_H
#define STRUCTURED_FINITE_VOLUME_DISC_ELEM_FIELD_H

#include "discretization.h"
#include "vert_field.h"

namespace structured_fv {
namespace disc {

// defines a field that stores n values in each element of each block (including ghosts)
template <typename T>
class ElementField
{
  public:
    using FieldData = Kokkos::View<T***, HostMemorySpace>;
    using ConstFieldData = Kokkos::View<const T***, HostMemorySpace>;

    ElementField(const StructuredDisc& disc, Int nvals_per_element) :
      m_disc(disc),
      m_nvals_per_element(nvals_per_element)
    {
      for (UInt i=0; i < disc.getNumBlocks(); ++i)
      {
        const StructuredBlock& block = disc.getBlock(i);
        auto dimensions = block.getCellDimensions();
        m_data.emplace_back("field_data", dimensions[0], dimensions[1], nvals_per_element);
      }
    }

    ElementField(const ElementField& other) = delete;

    ElementField& operator=(const ElementField& other) = delete;

    UInt getNumBlocks() const { return m_data.size(); }

    FieldData& getData(UInt block) { return m_data[block]; }

    const ConstFieldData& getData(UInt block) const { return m_data[block]; }

    T& operator()(UInt block, UInt i, UInt j, UInt v) { return m_data[block](i, j, v); }      

    const T& operator()(UInt block, UInt i, UInt j, UInt v) const { return m_data[block](i, j, v); }

    void set(const T& val);

    // Func is a callable object (Real x, Real y) -> std::array<T, num_vals_per_element>
    // return type can be anything of the correct length that supports operator[]
    template <typename Func>
    void set(Func func);

    // overwrite ghost values with the owner values
    void updateGhostValues();

  private:
    const StructuredDisc& m_disc;
    const UInt m_nvals_per_element;
    std::vector<FieldData> m_data;
};


template <typename T>
void ElementField<T>::set(const T& val)
{
  for (UInt b=0; b < m_disc.getNumBlocks(); ++b)
  {
    const StructuredBlock& block = m_disc.getBlock(b);
    auto data = getData(b);
    for (UInt i : block.getOwnedAndGhostCells().getXRange())
      for (UInt j : block.getOwnedAndGhostCells().getYRange())
        for (UInt d=0; d < m_nvals_per_element; ++d)
          data(i, j, d) = val;
  }
}

// Func is a callable object (Real x, Real y) -> std::array<T, num_vals_per_element>
// return type can be anything of the correct length that supports operator[]
template <typename T>
template <typename Func>
void ElementField<T>::set(Func func)
{
  for (UInt b=0; b < m_disc.getNumBlocks(); ++b)
  {
    const StructuredBlock& block = m_disc.getBlock(b);
    auto data = getData(b);
    auto coords = m_disc.getCoordField()->getData(b);
    for (UInt i : block.getOwnedAndGhostCells().getXRange())
      for (UInt j : block.getOwnedAndGhostCells().getYRange())
      {
        auto [x, y] = computeCellCentroid(coords, i, j);
        const auto& vals = func(x, y);
        for (UInt d=0; d < m_nvals_per_element; ++d)
          data(i, j, d) = vals[d];
      }
  }      
}


template <typename T>
void ElementField<T>::updateGhostValues()
{
  for (UInt i=0; i < m_disc.getNumBlockInterfaces(); ++i)
  {
    const StructuredBlockInterface& iface = m_disc.getBlockInterface(i);
    const StructuredBlock& blockL = m_disc.getBlock(iface.getBlockIdL());
    int num_ghost_cells = blockL.getNumGhostCellsPerDirection()[to_int(iface.getNeighborDirectionL())];
    const auto& indexerL = iface.getAdjBlockCellIndexerL();
    const auto& indexerR = iface.getAdjBlockCellIndexerR();
    auto fieldL = getData(iface.getBlockIdL());
    auto fieldR = getData(iface.getBlockIdR());
    
    for (UInt i : iface.getOwnedBoundaryCellsL().getXRange())
      for (UInt j : iface.getOwnedBoundaryCellsL().getYRange())
        for (int v=1; v <= num_ghost_cells; ++v)
        {
          auto [iprime, jprime] = computeIndices(iface.getNeighborDirectionL(), v, i, j);
          auto [ineighbor, jneighbor] = indexerL(iprime, jprime);

          for (UInt d=0; d < fieldL.extent(2); ++d)
          {
            fieldL(iprime, jprime, d) = fieldR(ineighbor, jneighbor, d);
          } 
        }

    //TODO: it would have better temporal locality to merge this
    //      into the above loops
    for (UInt i : iface.getOwnedBoundaryCellsR().getXRange())
      for (UInt j : iface.getOwnedBoundaryCellsR().getYRange())
        for (int v=1; v <= num_ghost_cells; ++v)
        {
          auto [iprime, jprime] = computeIndices(iface.getNeighborDirectionR(), v, i, j);
          auto [ineighbor, jneighbor] = indexerR(iprime, jprime);

          for (UInt d=0; d < fieldL.extent(2); ++d)
          {
            fieldR(iprime, jprime, d) = fieldL(ineighbor, jneighbor, d);
          }      
        }        
  }
}

}
}

#endif