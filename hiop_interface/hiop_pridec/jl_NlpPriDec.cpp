
#include "jl_NlpPriDec.hpp"


#include <cassert>
#include <cstring>  //for memcpy
#include <cstdio>
#include <cmath>


hiopSolveStatus JL_PriDecMasterProblem::solve_master(hiopVector& x,
                                                     const bool& include_r,
                                                     const double& rval /* = 0*/,
                                                     const double* grad /*=0*/,
                                                     const double* hess /*=0*/,
                                                     const char* master_options_file /*=nullptr*/)
{

  obj_ = -1e+20;
  hiopSolveStatus status;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // needs to fix to get the solver status
  opt_data.solve_base(this->get_recourse_gradient(), this->get_recourse_hessian()); //JL_solve_base_case(opt_data);

  double* x_vec = x.local_data();

  opt_data.getSolution(x_vec);
  obj_ = opt_data.getObjective();

  status=Solve_Success;
  if(status < 0) {
    printf("solver returned negative solve status: %d (with objective is %18.12e)\n", status, obj_);
    return status;
  }

  if(sol_ == nullptr) {
    sol_ = new double[n_];
  }

  memcpy(sol_, x_vec, n_ * sizeof(double));

  // send full solution to the contingency problems
  opt_data.send_solution();

  return Solve_Success;
};

bool JL_PriDecMasterProblem::eval_f_rterm(size_type idx, const int& n, const double* x, double& rval)
{
   //solve recourse

   opt_data.solve_contingency_recourse(idx, rval); 
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   return true;
};

// x is handled by primalDecomp to be the correct coupled x
bool JL_PriDecMasterProblem::eval_grad_rterm(size_type idx, const int& n, double* x, hiopVector& grad)
{

  assert(static_cast<int>(nc_) == n);
  double* grad_vec = grad.local_data();

   opt_data.getGradient(grad_vec);

  return true;
};

bool JL_PriDecMasterProblem::set_recourse_approx_evaluator(const int n,
                                                           hiopInterfacePriDecProblem::RecourseApproxEvaluator* evaluator)
{
  evaluator_= evaluator;
  return true;
}

