#include "lbs_nlS2_acceleration.h"
#include "IterativeMethods/lbs_iterativemethods.h"

#include "SweepChunks/lbs_sweepchunk_pwl_nlS2.h"

#include <chi_mpi.h>
#include <chi_log.h>
#include <ChiConsole/chi_console.h>

#include <math.h>

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
      if (omega[2] < 0.0)
        z_offset = 4;

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
    nlS2_moment_data.back().my_group_reference = group_set;
    nlS2_moment_data.back().omegas = group_set->quadrature->omegas;

    group_set->moment_callbacks.push_back(
        std::bind(&NlS2Acceleration::MomentCallBack::update_funcitonals,&nlS2_moment_data.back(),_1,_2,_3,_4,_5,_6));

    LBSGroupset* newgs = new LBSGroupset;
    group_sets_nlS2.push_back(newgs);
    //
    newgs->groups = group_set->groups;
    //TODO: aggregation type shouldn't matter here
    newgs->angleagg_method = LinearBoltzman::AngleAggregationType::SINGLE;
    //
    newgs->master_num_grp_subsets = group_set->master_num_grp_subsets;
    //
    // TODO: only if some uniform mesh option is set?
    //
    newgs->quadrature = std::make_shared<chi_math::AngularQuadrature>();
    //
    //InitializeWithCustom
//    void InitializeWithCustom(std::vector<double>& azimuthal,
//                              std::vector<double>& polar,
//                              std::vector<double>& in_weights,
//                              bool verbose=false) override;
    std::vector<double> azmuthal = {   M_PI/4., 3.*M_PI/4., 5.*M_PI/4., 7.*M_PI/4., M_PI/4., 3.*M_PI/4., 5.*M_PI/4., 7.*M_PI/4.};
    std::vector<double> polar = {M_PI/4.,    M_PI/4.,    M_PI/4.,    M_PI/4., 3.*M_PI/4., 3.*M_PI/4., 3.*M_PI/4., 3.*M_PI/4.};
    std::vector<double> weights(8, 1/8.);

   // static_cast<chi_math::ProductQuadrature*>(newgs->quadrature.get())->InitializeWithCustom(azmuthal,polar,weights,true);
    newgs->quadrature->InitializeWithCustom(azmuthal,polar,weights,true);
    //

    for (unsigned int i_dir=0; i_dir<newgs->quadrature->omegas.size(); i_dir++)
    {
      const auto & omega = newgs->quadrature->omegas[i_dir];
      std::cout<<"omega[0] = "<<omega[0]<<", omega[1] = "<<omega[1]<<", omega[2] = "<<omega[2]<<std::endl;
    }
    //
    newgs->allow_cycles = true;
    //
    newgs->iterative_method = NPT_CLASSICRICHARDSON;
    //
    //for (const auto & bndry: newgs->angle_agg->sim_boundaries);

//  double* Psi = &ref_boundaries[bndry_map]->boundary_flux[g];


  }

  InitializeParraysNLS2();


}

void NlS2Acceleration::Execute(){


  std::vector<SweepChunk*> sweep_chunks;
  std::vector<MainSweepScheduler> main_sweep_schedulers;

  std::vector<SweepChunk*> sweep_chunks_nlS2;
  std::vector<MainSweepScheduler> main_sweep_schedulers_nlS2;

  for (int gs=0; gs<group_sets.size(); gs++)
  {
    group_sets[gs]->BuildDiscMomOperator(options.scattering_order,
                                         options.geometry_type);
    group_sets[gs]->BuildMomDiscOperator(options.scattering_order,
                                         options.geometry_type);
    group_sets[gs]->BuildSubsets();

    ComputeSweepOrderings(group_sets[gs]);
    InitFluxDataStructures(group_sets[gs]);

    group_sets[gs]->max_iterations = 1;

    sweep_chunks.push_back(SetSweepChunk(gs));
    main_sweep_schedulers.push_back(
        MainSweepScheduler(SchedulingAlgorithm::DEPTH_OF_GRAPH,
            group_sets[gs]->angle_agg));

    //
    // nlS2 groups now
    //
  //  group_sets_nlS2[gs]->BuildDiscMomOperator(0, options.geometry_type);
  //  group_sets_nlS2[gs]->BuildMomDiscOperator(0, options.geometry_type);

    std::vector<double> d2m_op(8, 1/8.);
    std::vector<double> m2d_op(8, 1.);

    group_sets_nlS2[gs]->d2m_op.push_back(d2m_op);
    group_sets_nlS2[gs]->m2d_op.push_back(m2d_op);

    group_sets_nlS2[gs]->BuildSubsets();

    group_sets_nlS2[gs]->max_iterations = 100;

    ComputeSweepOrderings_nlS2(group_sets_nlS2[gs]);
    InitAngleAggSingle_nlS2(group_sets_nlS2[gs]);

    sweep_chunks_nlS2.push_back(SetSweepChunk_nlS2(gs));
    main_sweep_schedulers_nlS2.push_back(
        MainSweepScheduler(SchedulingAlgorithm::DEPTH_OF_GRAPH,
            group_sets_nlS2[gs]->angle_agg));
  }

  for (unsigned int iter = 0; iter<10; iter ++)
  {

    for (auto nlS2 : nlS2_moment_data)
    {
        for (unsigned int i_octant=0; i_octant<8; i_octant++){
          for (unsigned int i=0; i<nlS2.M_nlS2[i_octant].size(); i++){
             nlS2.phi_nlS2[i_octant][i] = 0.0;
             nlS2.M_nlS2[i_octant][i] = chi_mesh::Vector3(0.0,0.0,0.0);
          }
        }
    }

    //TODO clear the M vectors
    std::cout<<"========= iteration "<<iter<<" ========="<<std::endl;
    for (int gs=0; gs<group_sets.size(); gs++)
    {
      //
      LBSGroupset* group_set = group_sets[gs];
      //
//      if (group_set->iterative_method == NPT_CLASSICRICHARDSON)
//      {
        Solver::ClassicRichardson(gs, sweep_chunks[gs], main_sweep_schedulers[gs], true);
//      }
    }

    for (auto nlS2 : nlS2_moment_data)
    {
        for (unsigned int i_octant=0; i_octant<8; i_octant++){
          for (unsigned int i=0; i<nlS2.M_nlS2[i_octant].size(); i++){
            if (nlS2.phi_nlS2[i_octant][i] > TINY)
              nlS2.M_nlS2[i_octant][i] /= nlS2.phi_nlS2[i_octant][i];
          }
        }
    }

    bool verbose=false;
    if (verbose)
    {
      for (auto nlS2 : nlS2_moment_data)
      {
          for (unsigned int oct_index=0; oct_index<8; oct_index++){
            std::cout<<"================ Octant "<<oct_index<<"================\n";
            for (unsigned int i=0; i<nlS2.M_nlS2[oct_index].size(); i++){
              double norm = nlS2.M_nlS2[oct_index][i].x * nlS2.M_nlS2[oct_index][i].x;
              norm += nlS2.M_nlS2[oct_index][i].y * nlS2.M_nlS2[oct_index][i].y;
              norm += nlS2.M_nlS2[oct_index][i].z * nlS2.M_nlS2[oct_index][i].z;
              norm = sqrt(norm);

              std::cout<<"M_nlS2["<<i<<"] = "<<"("<<nlS2.M_nlS2[oct_index][i].x<<","<<nlS2.M_nlS2[oct_index][i].y<<","
                  <<nlS2.M_nlS2[oct_index][i].z<<")"<<", norm = "<<norm<<'\n';
            }
          }
      }
    }

    const auto scattering_order_store = options.scattering_order;
    options.scattering_order = 0;

    for (int gs=0; gs<group_sets_nlS2.size(); gs++)
    {
      ClassicRichardsonNLS2(gs, sweep_chunks_nlS2[gs], main_sweep_schedulers_nlS2[gs], true);
    }

    options.scattering_order = scattering_order_store;

  }

}

SweepChunk * NlS2Acceleration::SetSweepChunk_nlS2(int group_set_num){
  LBSGroupset* groupset = group_sets_nlS2[group_set_num];

  SweepChunk* sweep_chunk = new LBSSweepChunkPWL_nlS2(
        grid,                                    //Spatial grid of cells
        (SpatialDiscretization_PWL*)discretization, //Spatial discretization
        &cell_transport_views,                   //Cell transport views
        &phi_new_local,                          //Destination phi
        &q_moments_local,                        //Source moments
        groupset,                                //Reference groupset
        &material_xs,                            //Material cross-sections
        num_moments,max_cell_dof_count, nlS2_moment_data[group_set_num].M_nlS2);

  return sweep_chunk;

}

void NlS2Acceleration::MomentCallBack::update_funcitonals(const LinearBoltzman::CellViewFull * cell_view, int dof_index,int group, int m, int angle_num, double psi){

  if (m == 0){
    const int dof = dof_index*cell_view->num_grps + cell_view->dof_phi_map_start + group;
    const double moment_wieght = my_group_reference->d2m_op[m][angle_num];
    const int i_octant = angle_octant_map.at(angle_num);
    if (psi > 0)
    {
      phi_nlS2[i_octant][dof] += moment_wieght * psi;
      M_nlS2[i_octant][dof] += moment_wieght * psi * omegas[angle_num];
    }
  //  std::cout<<"i_octant = "<<i_octant<<", psi = "<<psi<<std::endl;
  //  std::cout<<"omega[0] = "<<omegas[angle_num][0]<<", omega[1] = "<<omegas[angle_num][1]<<", omega[2] = "<<omegas[angle_num][2]<<std::endl;

  }

}

void NlS2Acceleration::ComputeSweepOrderings_nlS2(LBSGroupset *groupset)
{
//  chi_log.Log(LOG_0)
//    << chi_program_timer.GetTimeString()
//    << " Computing Sweep ordering.\n";

  //============================================= Clear sweep ordering
  sweep_orderings_nlS2.clear();
  sweep_orderings_nlS2.shrink_to_fit();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Single angle aggr.
  for (auto& angle : groupset->quadrature->abscissae)
  {
    chi_mesh::sweep_management::SPDS* new_swp_order =
      chi_mesh::sweep_management::
      CreateSweepOrder(angle.theta,
                       angle.phi,
                       this->grid,
                       groupset->allow_cycles);
    this->sweep_orderings_nlS2.push_back(new_swp_order);
  }

//  chi_log.Log(LOG_0)
//    << chi_program_timer.GetTimeString()
//    << " Done computing nlS2 sweep orderings.           Process memory = "
//    << std::setprecision(3)
//    << chi_console.GetMemoryUsageInMB() << " MB";

}

void NlS2Acceleration::InitAngleAggSingle_nlS2(LBSGroupset *groupset)
{
  //
  // fix the isotropic boundary source for the S2 system
  //
  for (unsigned int ii = 0; ii<sweep_boundaries.size(); ++ii)
  {
    if (typeid(*sweep_boundaries[ii]) == typeid(chi_mesh::sweep_management::BoundaryIncidentHomogenous) )
    {
      std::vector<double> boundary_flux;
      boundary_flux = sweep_boundaries[ii]->boundary_flux;

      for (auto & flux : boundary_flux)
        flux *= M_PI/2.0;

      incident_P0_mg_boundaries_nlS2.push_back(boundary_flux);

      sweep_boundaries_nlS2.push_back( new chi_mesh::sweep_management::BoundaryIncidentHomogenous(incident_P0_mg_boundaries_nlS2.back()) );

    } else
    {
      sweep_boundaries_nlS2.push_back(sweep_boundaries[ii]);
    }
  }

  groupset->angle_agg->sim_boundaries  = sweep_boundaries_nlS2;

    groupset->angle_agg->number_of_groups        = groupset->groups.size();
    groupset->angle_agg->number_of_group_subsets = groupset->grp_subsets.size();
    groupset->angle_agg->quadrature              = groupset->quadrature;
    groupset->angle_agg->grid                    = grid;

    //=========================================== Set angle aggregation
    for (int q=0; q<1; q++)  //%%%%%%%%% Just a single group
    {
      auto angle_set_group = new TAngleSetGroup;
      groupset->angle_agg->angle_set_groups.push_back(angle_set_group);

      for (int n=0; n<groupset->quadrature->abscissae.size(); ++n)
      {
        bool make_primary = true;
        chi_mesh::sweep_management::PRIMARY_FLUDS* primary_fluds;

        for (int gs_ss=0; gs_ss<groupset->grp_subsets.size(); gs_ss++)
        {
          std::vector<int> angle_indices;

          angle_indices.push_back(n);

          chi_mesh::sweep_management::FLUDS* fluds;
          if (make_primary)
          {
            primary_fluds = new chi_mesh::sweep_management::
            PRIMARY_FLUDS(groupset->grp_subset_sizes[gs_ss]);

            primary_fluds->InitializeAlphaElements(sweep_orderings_nlS2[n]);
            primary_fluds->InitializeBetaElements(sweep_orderings_nlS2[n]);

            fluds = primary_fluds;
          }

          auto angleSet =
            new TAngleSet(groupset->grp_subset_sizes[gs_ss],
                          gs_ss,
                          sweep_orderings_nlS2[n],
                          fluds,
                          angle_indices,
                          sweep_boundaries_nlS2,
                          options.sweep_eager_limit,
                          &grid->GetCommunicator());

          angle_set_group->angle_sets.push_back(angleSet);
        }

      } //angle
    }//for q top


  }

}
