
#include "jl_NlpPriDec.hpp"


#include <cassert>
#include <cstring>  //for memcpy
#include <cstdio>
#include <cmath>


JL_PriDec::JL_PriDec(int ns_, int S_)
   : ns(ns_), S(S_), nc(ns_),
      evaluator_(nullptr)
{}


JL_PriDec::JL_PriDec(int ns_, int S_, bool include_)
    : ns(ns_), S(S_), nc(ns_), include_r(include_)
{
  if(include_r) {
    evaluator_ = new hiopInterfacePriDecProblem::RecourseApproxEvaluator(nc, "default");
  }
}

JL_PriDec::JL_PriDec(int ns_, int S_, bool include, hiopInterfacePriDecProblem::RecourseApproxEvaluator* evaluator)
   : ns(ns_), S(S_), include_r(include), evaluator_(evaluator)
{}

JL_PriDec::~JL_PriDec() {}

bool JL_PriDec::get_prob_sizes(size_type& n, size_type& m)
{
  n = ns;
  m = 0;
  return true;
}


bool JL_PriDec::get_cons_info(const size_type& m, double* clow, double* cupp, NonlinearityType* type)
{
  assert(m == 0);
  return true;
};

bool JL_PriDec::eval_cons(const size_type& n,
                          const size_type& m,
                          const size_type& num_cons,
                          const index_type* idx_cons,
                          const double* x,
                          bool new_x,
                          double* cons)
{
  assert(num_cons == 0);
  return true;
};


bool JL_PriDec::get_starting_point(const size_type& global_n, double* x0_)
{
  assert(global_n == ns);
  for(int i = 0; i < global_n; i++) x0_[i] = 2.;
  return true;
};

bool JL_PriDec::get_starting_point(const size_type& n,
                                   const size_type& m,
                                   double* x0_,
                                   bool& duals_avail,
                                   double* z_bndL0,
                                   double* z_bndU0,
                                   double* lambda0,
                                   bool& slacks_avail,
                                   double* ineq_slack)
{
  duals_avail = false;
  slacks_avail = false;
  return false;
};

bool JL_PriDec::get_MPI_comm(MPI_Comm& comm_out)
{
  comm_out = MPI_COMM_SELF;
  return true;
};

bool JL_PriDec::get_vars_info(const size_type& n, double* xlow, double* xupp, NonlinearityType* type)
{
  // assert(n>=4 && "number of variables should be greater than 4 for this example");
  assert(n == ns);
  // define x bounds
  for(int i = 0; i < ns; ++i) xlow[i] = 0.;
  for(int i = 0; i < ns; ++i) xupp[i] = +1e+20;
  for(int i = 0; i < ns; ++i) type[i] = hiopNonlinear;
  // uncoupled x fixed
  // for testing
  if(nc < ns) {
    for(int i = nc + 1; i < ns; ++i) xlow[i] = 1.;
    for(int i = nc + 1; i < ns; ++i) xupp[i] = 1.;
    xupp[0] = 1.;
    xupp[0] = 1.;
  }
  return true;
};

bool JL_PriDec::eval_Jac_cons(const size_type& n,
                              const size_type& m,
                              const size_type& num_cons,
                              const index_type* idx_cons,
                              const double* x,
                              bool new_x,
                              double* Jac)
{
  assert(m == 0);
  return true;
};

bool JL_PriDec::quad_is_defined()
{
  if(evaluator_ != NULL)
    return true;
  else
    return false;
};

bool JL_PriDec::set_quadratic_terms(const int& n, hiopInterfacePriDecProblem::RecourseApproxEvaluator* evaluator)
{
  assert(nc == n);
  evaluator_ = evaluator;

  return true;
};

bool JL_PriDec::set_include(bool include)
{
  include_r = include;
  return true;
};

hiopSolveStatus JL_PriDecMasterProblem::solve_master(hiopVector& x,
                                                     const bool& include_r,
                                                     const double& rval /* = 0*/,
                                                     const double* grad /*=0*/,
                                                     const double* hess /*=0*/,
                                                     const char* master_options_file /*=nullptr*/)
{

  obj_ = -1e+20;
  hiopSolveStatus status;
  if(my_nlp == NULL) {
     // my_nlp = new JL_PriDec(n_, S_, opt_data);
      my_nlp = new JL_PriDec(n_, S_);
  }


//  [[maybe_unused]] 
//bool ierr = 
my_nlp->set_include(include_r);


  if(include_r) {
    assert(my_nlp->quad_is_defined());
  }

  // check to see if the resource value and gradient are correct
  // printf("recourse value: is %18.12e)\n", rval_);
  hiopNlpDenseConstraints nlp(*my_nlp, master_options_file);

  // any of the options below can be overwritten by specifying them in the 'hiop_pridec_master.options' file

  nlp.options->SetStringValue("duals_update_type", "linear");
  nlp.options->SetStringValue("duals_init", "zero");  // "lsq" or "zero"
  nlp.options->SetStringValue("compute_mode", "cpu");
  nlp.options->SetStringValue("KKTLinsys", "xdycyd");
  nlp.options->SetStringValue("fixed_var", "relax");
  /*
  nlp.options->SetStringValue("dualsInitialization", "zero");
  nlp.options->SetStringValue("Hessian", "analytical_exact");
  nlp.options->SetStringValue("KKTLinsys", "xdycyd");
  nlp.options->SetStringValue("compute_mode", "hybrid");
  nlp.options->SetNumericValue("mu0", 1e-1);
  nlp.options->SetNumericValue("tolerance", 1e-5);
  */

  nlp.options->SetIntegerValue("verbosity_level", 1);

int rank;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

 
  // needs to fix to get the solver status
  opt_data.solve_base(my_nlp->get_recourse_gradient(), my_nlp->get_recourse_hessian()); //JL_solve_base_case(opt_data);

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
  my_nlp->set_quadratic_terms(n, evaluator);
  return true;
}

