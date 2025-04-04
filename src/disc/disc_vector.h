#ifndef STRUCTURED_FINITE_VOLUME_DISC_VECTOR_H
#define STRUCTURED_FINITE_VOLUME_DISC_VECTOR_H

#include "utils/project_defs.h"
#include "utils/traits.h"
#include "discretization.h"
#include "elem_field.h"
#include "utils/math.h"

namespace structured_fv {
namespace disc {

template <typename T>
class DiscVector
{
  public:
    explicit DiscVector(StructuredDiscPtr disc, const std::string& name) :
      m_disc(*disc),
      m_data(name, disc->getNumDofs())
    {}

    GlobalDof size() const { return m_data.extent(0); }

    T& operator()(GlobalDof dof) { return m_data(dof); }

    const T& operator()(GlobalDof dof) const { return m_data(dof); }

    void set(const T& val)
    {
      for (GlobalDof i=0; i < size(); ++i)
        m_data(i) = val;
    }

    template <typename Fxy, IsFuncXY_t<Fxy> = true>
    void set(Fxy func)
    {
      for (UInt block_id : m_disc.getRegularBlocksIds())
      {
        const StructuredBlock& block = m_disc.getBlock(block_id);
        const auto& dof_nums = m_disc.getDofNumbering()->getData(block_id);
        const auto& coords   = m_disc.getCoordField()->getData(block_id);

        for (UInt i : block.getOwnedCells().getXRange())
          for (UInt j : block.getOwnedCells().getYRange())
          {
             Vec2<Real> x = computeCellCentroid(coords, i, j);
             auto vals = func(x[0], x[1]);

             for (int k=0; k < m_disc.getNumDofsPerNode(); ++k)
               m_data(dof_nums(i, j, k)) = vals[k];
          }
      }
    }

  private:
    const StructuredDisc& m_disc;
    Kokkos::View<T*, HostMemorySpace> m_data;
};

template <typename T>
using DiscVectorPtr = std::shared_ptr<DiscVector<T>>;

}
}

#endif