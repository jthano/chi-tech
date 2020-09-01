#include "lbs_nlS2_acceleration.h"


void LinearBoltzman::NlS2Acceleration::InitializeParraysNLS2()
{

  auto pwl_discretization = (SpatialDiscretization_PWL*)discretization;

  auto domain_ownership = pwl_discretization->OrderNodesDFEM(grid);
  local_dof_count = domain_ownership.first;
  glob_dof_count  = domain_ownership.second;
  //================================================== Compute num of unknowns
  int num_grps = groups.size();
  int M = num_moments;
  unsigned long long local_unknown_count = local_dof_count * num_grps;
  //================================================== Size local vectors
  //
  for (auto & nlS2 : nlS2_moment_data)
  {
    for (unsigned int i=0; i<8; i++){
      nlS2.phi_nlS2[i].resize(local_unknown_count,0.0);
      nlS2.M_nlS2[i].resize(local_unknown_count);
    }
  }
}
