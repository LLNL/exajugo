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
#include <climits>    // for INT_MAX

/**t
 * Driver for PriDec Example 1 that illustrates the use of hiop::hiopAlgPrimalDecomposition
 *
 * @note This example is built only when HIOP_USE_MPI is enabled during cmake build
 * and require at least two MPI ranks in MPI_COMM_WORLD.
 *
 */

#include <filesystem>   // for std::filesystem
//const char preferred_separator = '/';
//namespace fs = std::filesystem;

int main(int argc, char** argv)
{
   std::string instance;


  int rank = 0;
#ifdef HIOP_USE_MPI
  MPI_Init(&argc, &argv);
  int comm_size;
  int ierr = MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  assert(MPI_SUCCESS == ierr);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  assert(MPI_SUCCESS == ierr);
#endif

//std::string sep(1, fs::path::preferred_separator); // convert char to string
std::string sep(1, preferred_separator); // convert char to string

std::string defoutput = "output"+sep+"rank_" + std::to_string(rank)+sep;
// Assume 'rank' is already defined as an int
std::string outputDir = std::getenv("OUTPUT_ITER") ? std::getenv("OUTPUT_ITER"): defoutput;
if (!outputDir.empty() && outputDir.back() != preferred_separator  && outputDir.back() != '\\') {
    outputDir += preferred_separator;
}

// Create the directory (including parent directories if needed)
//fs::create_directories(outputDir);

#ifdef HIOP_USE_MAGMA
  magma_init();
#endif

    if (argc > 1) {
       instance = argv[1];
       if (rank==0)
          std::cout << " Instance: " << instance << std::endl;

    } else {
       if (rank==0)
          std::cout << " Instance not provided! Execution aborted!\n\n" << std::endl;
      exit(0);
    }
   
  jl_init();

  int max_iter = INT_MAX;
  const char* env_max_iter = std::getenv("MAX_ITER");
  if (env_max_iter) 
  {
    max_iter = std::atoi(env_max_iter);
    }

  // JL_Interface constructor: base system and maximum number of iterations
  JL_Interface prob_data(outputDir, instance, max_iter);

  int ncont = prob_data.number_of_contingencies(); //6//20;
  if (rank==0)
     std::cout<<" # of contingencies: "<<ncont<<" comm_size: "<<comm_size<<"\n\n";
  if ((comm_size < 2) || (comm_size > ncont+1))
  {
       if (rank==0)
          std::cout << " Total number of processes must be >=2 and <= "<<ncont+1<<"! Execution aborted!\n\n" << std::endl;
      exit(0);
  }

  int nc = prob_data.number_of_columns(); //6//20;

  int* list = new int[nc];
  for(int i = 0; i < nc; i++) list[i] = i;

  JL_PriDecMasterProblem pridec_problem(prob_data);

 std::cout<<"\n # of NC: "<<nc<<"\n\n";

  hiop::hiopAlgPrimalDecomposition pridec_solver(&pridec_problem, nc, list, MPI_COMM_WORLD);

  pridec_solver.set_max_iteration(prob_data.get_max_iter()); // Set maximum iterations

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


