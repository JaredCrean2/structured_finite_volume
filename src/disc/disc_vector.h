#ifndef STRUCTURED_FINITE_VOLUME_DISC_VECTOR_H
#define STRUCTURED_FINITE_VOLUME_DISC_VECTOR_H

#include "utils/project_defs.h"

namespace structured_fv {
namespace disc {

template <typename T>
class DiscVector
{
  public:
    explicit DiscVector(GlobalDof num_dofs, const std::string& name) :
      m_data(name, num_dofs)
    {}

    GlobalDof size() const { return m_data.extent(0); }

    T& operator()(GlobalDof dof) { return m_data(dof); }

    const T& operator()(GlobalDof dof) const { return m_data(dof); }

    void set(const T& val)
    {
      for (GlobalDof i=0; i < size(); ++i)
        m_data(i) = val;
    }

  private:
    Kokkos::View<T*, HostMemorySpace> m_data;
};

}
}

#endif