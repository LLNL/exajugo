// the problem to be solved


#include "jl_NlpPriDec.hpp"
#include "jl_interface.hpp"

// the solver
#include "hiopAlgPrimalDecomp.hpp"


#include <iostream>
#include <julia.h>

JULIA_DEFINE_FAST_TLS // Required for thread-local storage in Julia


#ifdef HIOP_USE_MAGMA
#include "magma_v2.h"
#endif

#include <cstdlib>
#include <string>


/**t
 * Driver for PriDec Example 1 that illustrates the use of hiop::hiopAlgPrimalDecomposition
 *
 * @note This example is built only when HIOP_USE_MPI is enabled during cmake build
 * and require at least two MPI ranks in MPI_COMM_WORLD.
 *
 */



int main(int argc, char** argv)
{


  int rank = 0;
#ifdef HIOP_USE_MPI
  MPI_Init(&argc, &argv);
  int comm_size;
  int ierr = MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  assert(MPI_SUCCESS == ierr);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  assert(MPI_SUCCESS == ierr);
#endif

#ifdef HIOP_USE_MAGMA
  magma_init();
#endif

   jl_init();

  // include_jl_functions();

   JL_Interface prob_data;

//  int S = prob_data.number_of_contingencies(); //3; //100;
  int nc = prob_data.number_of_columns(); //6//20;

  int* list = new int[nc];
  for(int i = 0; i < nc; i++) list[i] = i;

  
  JL_PriDecMasterProblem pridec_problem(prob_data);
  hiop::hiopAlgPrimalDecomposition pridec_solver(&pridec_problem, nc, list, MPI_COMM_WORLD);

  auto status = pridec_solver.run();

  if(status != Solve_Success) {
    if(rank == 0) printf("Solve was NOT successfull.");
  } else {
    if(rank == 0) printf("Solve was successfull. Optimal value: %12.5e\n", pridec_solver.getObjective());
  }

  delete[] list;

 jl_atexit_hook(0);

#ifdef HIOP_USE_MAGMA
  magma_finalize();
#endif
#ifdef HIOP_USE_MPI
  MPI_Finalize();
#endif

  printf("Returned successfully from driver! Rank=%d\n", rank);

  return 0;
}


