#include "ChiLua/chi_lua.h"

#include "../lbs_linear_boltzman_solver.h"
#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

//###################################################################
/**Executes the LBS solver.
\param SolverIndex int Handle to the solver.
 \ingroup LuaNPT
 */
int chiLBSExecute(lua_State *L)
{
  int solver_index = lua_tonumber(L,1);

  //============================================= Get pointer to solver
  chi_physics::Solver* psolver;
  LinearBoltzman::Solver* solver;
  try{
    psolver = chi_physics_handler.solver_stack.at(solver_index);
    solver = dynamic_cast<LinearBoltzman::Solver*>(psolver);
    if (solver == nullptr)
    {
      fprintf(stderr,"ERROR: Incorrect solver-type"
                     "in chiLBSExecute\n");
      exit(EXIT_FAILURE);
    }
  }
  catch(const std::out_of_range& o)
  {
    fprintf(stderr,"ERROR: Invalid handle to solver"
                   "in chiLBSExecute\n");
    exit(EXIT_FAILURE);
  }

  solver->Execute();

  return 0;
}
