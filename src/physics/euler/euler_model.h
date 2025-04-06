#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_MODEL_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_MODEL_H

#include "disc/discretization.h"
#include "physics/euler/euler_flux.h"
#include "utils/project_defs.h"
#include "physics/physics_model.h"

#include "typedefs.h"
#include "hlle_flux.h"
#include "lax_friedrich_flux.h"
#include "disc/face_field.h"
#include "disc/elem_field.h"

#include <iostream>

namespace structured_fv {
namespace euler {

struct EulerOpts {};

class EulerModel : public PhysicsModel
{
  public:
    EulerModel(const EulerOpts& opts, StructuredDiscPtr disc,
               const std::vector<Fxyt>& bc_functions,
               Fxyt source_func) :
      //m_opts(opts),
      m_disc(disc),
      m_bc_functions(bc_functions),
      m_source_func(source_func),
      m_solution(std::make_shared<disc::ElementField<Real>>(*disc, DofsPerCell)),
      m_residual(std::make_shared<disc::ElementField<Real>>(*disc, DofsPerCell))
    {
      if (disc->getNumDofsPerNode() != DofsPerCell)
        throw std::runtime_error("Euler can only have 4 dof per node");

      if (m_bc_functions.size() != disc->getNumGhostBCBlocks())
        throw std::runtime_error("number of BC functions must be equal to number of BC blocks");
    }

    void evaluateRhs(DiscVectorPtr<Real> q, Real t, DiscVectorPtr<Real> residual) override;

    const disc::StructuredDiscPtr& getDisc() const override { return m_disc; }


  private:
    Fxyt& getBCFunction(UInt block_id);

    void setBCValues(ElementFieldPtr<Real> solution, Real t);

    template <typename Flux, typename Tag>
    void evaluateInterfaceTerms(const ElementFieldPtr<Real>& solution, Real t,
                                Flux& flux_func, Tag dir_tag, ElementFieldPtr<Real> residual);

    void evaluateSourceTerm(Real t, ElementFieldPtr<Real> residual);

    void checkPositivity(const ElementFieldPtr<Real>& solution);


    //EulerOpts m_opts;
    StructuredDiscPtr m_disc;
    std::vector<Fxyt> m_bc_functions;
    Fxyt m_source_func;

    ElementFieldPtr<Real> m_solution;
    ElementFieldPtr<Real> m_residual;
};

inline Vec4<Real> getValues(const disc::ElementField<Real>::ConstFieldData& data, UInt i, UInt j)
{
  assert(data.extent(2) == DofsPerCell);
  return {data(i, j, 0), data(i, j, 1), data(i, j, 2), data(i, j, 3)};
}

}
}

#endif