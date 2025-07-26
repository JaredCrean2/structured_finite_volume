#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_MODEL_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_MODEL_H

#include "disc/discretization.h"
#include "physics/common/slope_limiters.h"  // TODO: separate out enums from defs
#include "utils/project_defs.h"
#include "physics/physics_model.h"
#include "flux_function_enums.h"
#include "reconstruction_enum.h"
#include "typedefs.h"

namespace structured_fv {
namespace euler {

//TODO: disable source term if it is zero
struct EulerOpts
{
  common::SlopeLimiter limiter = common::SlopeLimiter::FirstOrder;
  Reconstruction recon = Reconstruction::Conservative;
  FluxFunction flux = FluxFunction::Roe;
  Real roe_efix_delta = 0.3;
};

template <typename T>
struct Fields
{
  Fields(const disc::StructuredDisc& disc) :
    solution(std::make_shared<disc::ElementField<T>>(disc, 4)),
    residual(std::make_shared<disc::ElementField<T>>(disc, 4))
  {}
  
  ElementFieldPtr<T> solution;
  ElementFieldPtr<T> residual;
};

class EulerModel : public PhysicsModel
{
  public:
    EulerModel(const EulerOpts& opts, StructuredDiscPtr disc,
               const std::vector<Fxyt>& bc_functions,
               Fxyt source_func) :
      m_opts(opts),
      m_disc(disc),
      m_fields_real(*disc),
      m_fields_complex(*disc),
      m_bc_functions(bc_functions),
      m_source_func(source_func)
      //m_solution(std::make_shared<disc::ElementField<Real>>(*disc, DofsPerCell)),
      //m_residual(std::make_shared<disc::ElementField<Real>>(*disc, DofsPerCell))
    {
      if (disc->getNumDofsPerNode() != DofsPerCell)
        throw std::runtime_error("Euler can only have 4 dof per node");

      if (m_bc_functions.size() != disc->getNumGhostBCBlocks())
        throw std::runtime_error("number of BC functions must be equal to number of BC blocks");
    }

    void evaluateRhs(DiscVectorPtr<Real> q, Real t, DiscVectorPtr<Real> residual) override;

    void evaluateJacobian(disc::DiscVectorPtr<Real> q, Real t, disc::DiscVectorPtr<Real> residual, linear_system::AssemblerBasePtr assembler) override;

    // computes h = dR/dq * v without forming dR/dq
    void computeJacVecProduct(disc::DiscVectorPtr<Real> q, Real t, disc::DiscVectorPtr<Real> v, disc::DiscVectorPtr<Real> h) override;

    std::shared_ptr<linear_system::SparsityPattern> getSparsityPattern() const override;

    const disc::StructuredDiscPtr& getDisc() const override { return m_disc; }


  private:

    template <typename T>
    void evaluateRhsT(Fields<T>& fields, Real t);

    Fxyt& getBCFunction(UInt block_id);

    template <typename T>
    void setBCValues(ElementFieldPtr<T> solution, Real t);

    template <typename T>
    void evaluateSourceTerm(Real t, Fields<T>& fields);

    template <typename T>
    void checkPositivity(const ElementFieldPtr<T>& solution);

    template <typename T>
    T computeMaxWaveSpeed(ElementFieldPtr<T> solution);


    EulerOpts m_opts;
    StructuredDiscPtr m_disc;
    Fields<Real> m_fields_real;
    Fields<Complex> m_fields_complex;    
    std::vector<Fxyt> m_bc_functions;
    Fxyt m_source_func;

    //ElementFieldPtr<Real> m_solution;
    //ElementFieldPtr<Real> m_residual;
};

//template <typename T>
//inline Vec4<T> getValues(const typename disc::ElementField<T>::ConstFieldData& data, UInt i, UInt j)

template <typename ViewType>
inline Vec4<typename ViewType::value_type> getValues(const ViewType& data, UInt i, UInt j)
{
  static_assert(Kokkos::is_view_v<ViewType>);
  using T = typename ViewType::value_type;
  assert(data.extent(2) == DofsPerCell);
  return Vec4<T>{data(i, j, 0), data(i, j, 1), data(i, j, 2), data(i, j, 3)};
}

}
}

#endif