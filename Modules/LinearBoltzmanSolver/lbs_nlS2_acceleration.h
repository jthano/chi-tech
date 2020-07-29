#ifndef _lbs_power_iteration_h
#define _lbs_power_iteration_h

#include "lbs_linear_boltzman_solver.h"
#include "SweepChunks/lbs_sweepchunk_pwl.h"

#include "ChiMath/Quadratures/product_quadrature.h"

#include <functional>
#include <map>

namespace LinearBoltzman
{
class NlS2Acceleration : public Solver
{
private:
  const double TINY=1.e-8;
  std::map< LBSGroupset *, std::map<unsigned int,unsigned int> > angular_quad_s2_map;
  std::shared_ptr<chi_math::ProductQuadrature> quad_uniform_ref;
  std::vector<std::shared_ptr<chi_math::AngularQuadrature>> quad_store;

public:
  std::vector<LBSGroupset *> group_sets_nlS2;

 // Vec  phi_old_k;
  std::vector<double> phi_old_k_local;
  std::vector<double> delta_phi_local_k;

  class MomentCallBack{
  public:
    LBSGroupset* my_group_reference;
    std::vector<std::vector<double>> phi_nlS2;
    std::vector<std::vector<chi_mesh::Vector3>> M_nlS2;
    std::map<int,int> angle_octant_map;
    std::vector<chi_mesh::Vector3> omegas;
    LBSGroupset::MomentCallbackF s2_integrals;
    void update_funcitonals(int dof_index, int m, int angle_num, double psi);
    MomentCallBack(LBSGroupset* my_group_reference):my_group_reference(my_group_reference),phi_nlS2(8),M_nlS2(8){};
  };

  std::vector<MomentCallBack> nlS2_moment_data;

  std::vector<chi_mesh::sweep_management::SPDS*> sweep_orderings_store;
  std::vector<chi_mesh::sweep_management::SPDS*> uniform_s2_orderings_store;

  NlS2Acceleration();

  virtual void Initialize() override;

  virtual void Execute() override;

//  void SetSource(int group_set_num,
//      bool apply_mat_src = false,
//      bool suppress_phi_old = false) override;

  void InitializeParrays() override;

  void ClassicRichardsonNLS2(int group_set_num);

  void ComputeSweepOrderings_nlS2(LBSGroupset *groupset);

  void InitFluxDataStructures_nlS2(LBSGroupset *groupset);

};

}




#endif
