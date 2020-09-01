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
  chi_math::ProductQuadrature quad_uniform_ref;

public:
  std::vector<LBSGroupset *> group_sets_nlS2;
  std::vector<chi_mesh::sweep_management::SPDS*> sweep_orderings_nlS2;
  std::vector<SweepBndry*>                      sweep_boundaries_nlS2;
  std::vector<std::vector<double>>              incident_P0_mg_boundaries_nlS2;

  class MomentCallBack{
  public:
    LBSGroupset* my_group_reference;
    std::vector<std::vector<double>> phi_nlS2;
    std::vector<std::vector<chi_mesh::Vector3>> M_nlS2;
    std::map<int,int> angle_octant_map;
    std::vector<chi_mesh::Vector3> omegas;
    LBSGroupset::MomentCallbackF s2_integrals;
    void update_funcitonals(const LinearBoltzman::CellViewFull *, int dof_index, int group, int m, int angle_num, double psi);
    MomentCallBack(LBSGroupset* my_group_reference):my_group_reference(my_group_reference),phi_nlS2(8),M_nlS2(8){};
  };

  std::vector<MomentCallBack> nlS2_moment_data;

  NlS2Acceleration();

  virtual void Initialize() override;

  virtual void Execute() override;

//  void SetSource(int group_set_num,
//      bool apply_mat_src = false,
//      bool suppress_phi_old = false) override;

  void InitializeParraysNLS2();

  void ClassicRichardsonNLS2(int group_set_num, SweepChunk* sweep_chunk, MainSweepScheduler & sweepScheduler, bool log_info=true);

  void ComputeSweepOrderings_nlS2(LBSGroupset *groupset);

  SweepChunk *SetSweepChunk_nlS2(int group_set_num);

  void InitAngleAggSingle_nlS2(LBSGroupset *groupset);

};

}




#endif
