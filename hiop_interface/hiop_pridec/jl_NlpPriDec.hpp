

#ifndef HIOP_EXAMPLE_PRIDEC_EX1
#define HIOP_EXAMPLE_PRIDEC_EX1

#include <julia.h>

#include "jl_interface.hpp"

//JULIA_DEFINE_FAST_TLS // Required for thread-local storage in Julia


#include "hiopInterfacePrimalDecomp.hpp"

#include "hiopInterface.hpp"
#include "hiopNlpFormulation.hpp"
#include "hiopAlgFilterIPM.hpp"

#include <cassert>
#include <cstring>  //for memcpy
#include <cstdio>
#include <cmath>
#include <chrono>


using namespace hiop;

/**
 * Master problem class based on the base case problem, which is a Ex8 class.
 *
 */
class JL_PriDecMasterProblem : public hiopInterfacePriDecProblem
{

public:

JL_Interface opt_data;

  JL_PriDecMasterProblem(const JL_Interface& _opt_data)
      : 
        opt_data(_opt_data),
        obj_(-1e20),
        sol_(nullptr), evaluator_(nullptr)
  {
 
     std::cout<<"!11!!\n";

    n_ = opt_data.number_of_columns();
    S_ = opt_data.number_of_contingencies();
    nc_=n_;

  }

  virtual ~JL_PriDecMasterProblem()
  {
    delete[] sol_;

  }

  virtual hiopSolveStatus solve_master(hiopVector& x,
                                       const bool& include_r,
                                       const double& rval = 0,
                                       const double* grad = 0,
                                       const double* hess = 0,
                                       const char* master_options_file = nullptr);

  /**
   * This function returns the recourse objective, which is 0.5*(x+Se_i)(x+Se_i).
   */
  virtual bool eval_f_rterm(size_type idx, const int& n, const double* x, double& rval);

  /**
   * This function returns the recourse gradient.
   */
  virtual bool eval_grad_rterm(size_type idx, const int& n, double* x, hiopVector& grad);

  /**
   * This function sets up the approximation of the recourse objective based on the function value and gradient
   * returned by eval_f_rterm and eval_grad_rterm.
   * Implemented with alpha = 1 for now only.
   * This function is called only if quadratic regularization is included.
   */
  virtual bool set_recourse_approx_evaluator(const int n, hiopInterfacePriDecProblem::RecourseApproxEvaluator* evaluator);
  // Returns the number S of recourse terms.
  size_type get_num_rterms() const { return S_; }
  size_type get_num_vars() const { return n_; }
  // Returns the solution.
  void get_solution(double* x) const
  {
    for(int i = 0; i < static_cast<int>(n_); i++) x[i] = sol_[i];
  }
  double get_objective() { return obj_; }

private:
  int iter;
  size_type n_;
  size_type S_;
  size_type nc_;

  double obj_;
  double* sol_;

hiopInterfacePriDecProblem::RecourseApproxEvaluator* evaluator_;

double *get_recourse_gradient() const
{
   if (evaluator_ == nullptr) return nullptr;

    hiopVector* hograd = evaluator_->get_rgrad(); 
    return hograd->local_data();

}

double *get_recourse_hessian() const
{
   if (evaluator_ == nullptr) return nullptr;

    hiopVector* hograd = evaluator_->get_rhess(); 
    return hograd->local_data();
}


};

#endif
