#include "lbs_nlS2_acceleration.h"
#include "IterativeMethods/lbs_iterativemethods.h"

#include <chi_mpi.h>
#include <chi_log.h>
#include <ChiConsole/chi_console.h>

extern ChiMPI&      chi_mpi;
extern ChiLog&     chi_log;
extern ChiConsole&  chi_console;

#include <iomanip>


namespace LinearBoltzman
{

NlS2Acceleration::NlS2Acceleration()
{

}

void NlS2Acceleration::Initialize(){
  Solver::Initialize();

  using namespace std::placeholders;

  for (auto group_set : group_sets)
  {
    nlS2_moment_data.push_back(MomentCallBack(group_set));

    for (unsigned int i_dir=0; i_dir<group_set->quadrature->omegas.size(); i_dir++)
    {
      const auto & omega = group_set->quadrature->omegas[i_dir];

      unsigned int z_offset = 0;
      if (omega[3] < 0.0)
        z_offset = 3;

      if (omega[0]==0 || omega[1]==0){
        chi_log.Log(LOG_ALLERROR)
          << "NlS2Acceleration::Initialize - A direction in the quadrature set has a value of 0"<<
         " for the x or y component which was not yet considered.";
        exit(EXIT_FAILURE);
      }

      unsigned int octant;
      if (omega[0]>0 && omega[1]>0)
        octant=0;
      else if (omega[0]<0 && omega[1]>0)
        octant=1;
      else if (omega[0]<0 && omega[1]<0)
        octant=2;
      else
        octant=3;

      octant += z_offset;

      nlS2_moment_data.back().angle_octant_map.insert( std::pair<unsigned int,unsigned int>(i_dir,octant) );
    }
    nlS2_moment_data.back().omegas = group_set->quadrature->omegas;

    group_set->moment_callbacks.push_back(
        std::bind(&NlS2Acceleration::MomentCallBack::update_funcitonals,&nlS2_moment_data.back(),_1,_2,_3,_4));

    LBSGroupset* newgs = new LBSGroupset;
    group_sets_nlS2.push_back(newgs);
    newgs->groups = group_set->groups;
    //TODO: aggregation type shouldn't matter here
    newgs->angleagg_method = LinearBoltzman::AngleAggregationType::SINGLE;
    //
    newgs->master_num_grp_subsets = group_set->master_num_grp_subsets;
    //



  }

  // TODO: if uniform mesh option set
  quad_uniform_ref->InitializeWithGLC(1,1,false);


}

void NlS2Acceleration::Execute(){


  std::vector<SweepChunk*> sweep_chunks;
  std::vector<MainSweepScheduler> main_sweep_schedulers;

  for (int gs=0; gs<group_sets.size(); gs++)
  {
    group_sets[gs]->BuildDiscMomOperator(options.scattering_order,
                                         options.geometry_type);
    group_sets[gs]->BuildMomDiscOperator(options.scattering_order,
                                         options.geometry_type);
    group_sets[gs]->BuildSubsets();

    ComputeSweepOrderings(group_sets[gs]);
    InitFluxDataStructures(group_sets[gs]);

    LBSGroupset* group_set = group_sets[gs];

    sweep_chunks.push_back(SetSweepChunk(gs));
    main_sweep_schedulers.push_back(
        MainSweepScheduler(SchedulingAlgorithm::DEPTH_OF_GRAPH,
        group_set->angle_agg));

  }

//  if (groupset->angleagg_method == LinearBoltzman::AngleAggregationType::SINGLE)
//  {
//    for (auto& angle : groupset->quadrature->abscissae)
//    {
//      chi_mesh::sweep_management::SPDS* new_swp_order =
//        chi_mesh::sweep_management::
//        CreateSweepOrder(angle.theta,
//                         angle.phi,
//                         this->grid,
//                         groupset->allow_cycles);
//      this->sweep_orderings.push_back(new_swp_order);
//    }
//  }


  for (auto& angle : quad_uniform_ref->abscissae)
  {
    chi_mesh::sweep_management::SPDS* new_swp_order =
      chi_mesh::sweep_management::
      CreateSweepOrder(angle.theta,
                       angle.phi,
                       this->grid,
                             true);
  }

  for (unsigned int iter = 0; iter<2; iter ++)
  {
    for (int gs=0; gs<group_sets.size(); gs++)
    {
      //
      LBSGroupset* group_set = group_sets[gs];
      //
//      if (group_set->iterative_method == NPT_CLASSICRICHARDSON)
//      {
        Solver::ClassicRichardson(gs, sweep_chunks[gs], main_sweep_schedulers[gs]);
//      }

    }

    for (auto nlS2 : nlS2_moment_data)
    {
        for (unsigned int oct_index=0; oct_index<8; oct_index++){
          for (unsigned int i=0; i<nlS2.M_nlS2[oct_index].size(); i++){
            if (nlS2.phi_nlS2[oct_index][i] > TINY)
              nlS2.M_nlS2[oct_index][i] /= nlS2.phi_nlS2[oct_index][i];
          }
        }
    }


  }




//  for (unsigned int oct_index=0; oct_index<8; oct_index++){
//    for (unsigned int i=0; i<M_nls2[oct_index].size(); i++){
//      if (phi_nls2[oct_index][i] > TINY)
//        M_nls2[oct_index][i] /= phi_nls2[oct_index][i];
//    }
//  }

//  sweep_orderings_store = std::move(sweep_orderings);
//
//  for (unsigned int gs = 0; gs<group_sets.size(); gs++){
//    group_sets[gs]->quadrature = uniform_s2_orderings_store;
//  }

}

void NlS2Acceleration::MomentCallBack::update_funcitonals(int dof_index, int m, int angle_num, double psi){

  if (m == 0){
    const double moment_wieght = my_group_reference->d2m_op[m][angle_num];
    const int i_octant = angle_octant_map.at(angle_num);
    phi_nlS2[i_octant][dof_index] += moment_wieght * psi;
    M_nlS2[i_octant][dof_index] += moment_wieght * psi * omegas[angle_num];
  }

}

}
