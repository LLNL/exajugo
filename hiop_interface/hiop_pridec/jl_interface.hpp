
#ifndef JL_INTERFACE_HPP

#define JL_INTERFACE_HPP

#include <string>
#include <iostream>

#include <julia.h>


extern jl_function_t* jl_load_ACOPF;
extern jl_function_t* jl_copy_ACOPF;

extern jl_function_t* jl_number_of_contingencies;
extern jl_function_t* jl_number_of_columns;

extern jl_function_t* jl_solve_base_case ;

extern jl_function_t* jl_solve_base_case_recourse;
extern jl_function_t* jl_get_recourse_derivatives;

extern jl_function_t* jl_getModel;
extern jl_function_t* jl_getDim;
extern jl_function_t* jl_getObjective;
extern jl_function_t* jl_getSolution;

extern jl_function_t* jl_solve_contingency_pridec;
extern jl_function_t* jl_getCost;
extern jl_function_t* jl_getGradient;
    
extern jl_function_t*  deepcopy_func;
extern jl_function_t*  jl_serialize_obj;
extern jl_function_t*  jl_deserialize_obj;

extern jl_function_t*  jl_get_data_size;
extern jl_function_t*  jl_get_data_bytes;

extern jl_function_t* jl_debug_base_case;
extern jl_function_t* jl_save_solution;
extern jl_function_t* jl_save_cont_solution;
extern jl_function_t* jl_save_array;


void include_jl_functions();


class JL_Interface
{
 int source_process;

    int rank;
    int nproc;
    int size_buffer;


protected:

   jl_value_t* opt_data;
   jl_value_t* base_sol;
   jl_value_t* cont_sol;
    uint8_t* send_buffer;
    std::string instance;

    jl_value_t* read_data() 
    {

       std::string exajugo_path = std::getenv("PATH_TO_EXAJUGO");
       std::string example_path = exajugo_path+"/examples/"+instance+"/";
       std::string userCommand= "rm -Rf "+instance+"; mkdir -p "+instance;
       int result = system(userCommand.c_str());

       std::string raw_file = (example_path+"case.raw"); 
       std::string rop_file = (example_path+"case.rop"); 
       std::string con_file = (example_path+"case.con");
      
       jl_value_t* jl_raw_file = jl_cstr_to_string(raw_file.c_str()); 
       jl_value_t* jl_rop_file = jl_cstr_to_string(rop_file.c_str()); 
       jl_value_t* jl_con_file = jl_cstr_to_string(con_file.c_str());

       return jl_call3(jl_load_ACOPF, jl_raw_file, jl_rop_file, jl_con_file);
       }

    void release_buffer();
    uint8_t *alloc_buffer(int);
 

  public:

   ~JL_Interface() { release_buffer(); }
    

    int getRank(){ return rank;}

     jl_value_t* ptr() const { return opt_data; } 
    
    void getCost(double& rval) { rval =  jl_unbox_float64(jl_call1(jl_getCost, cont_sol)); }
   
    jl_value_t* solve_base_case_with_recourse(jl_value_t* prevsol, double *grad, double *hess) 
    {
       jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);
       jl_array_t* jl_grad = jl_ptr_to_array_1d(array_type, grad, getDim(), 0);

       // this line is not needed: can use the same type above 'array_type'
       jl_value_t* array_typeh = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);
       jl_array_t* jl_hess = jl_ptr_to_array_1d(array_typeh, hess, getDim(), 0);

       jl_value_t* ptr_rderivatives = jl_call2(jl_get_recourse_derivatives, (jl_value_t *)jl_grad, (jl_value_t *)jl_hess);

       jl_value_t* retptr= jl_call3(jl_solve_base_case_recourse, opt_data, prevsol, ptr_rderivatives);
 
       // save gradient
       std::string fname=instance+"/gradient.csv";
       jl_value_t*  jl_fname = jl_cstr_to_string(fname.c_str());
       jl_call3(jl_save_array, jl_fname, opt_data, (jl_value_t *)jl_grad);

       // save hessian
       std::string fnameh=instance+"/hessian.csv";
       jl_value_t*  jl_fnameh = jl_cstr_to_string(fnameh.c_str());
       jl_call3(jl_save_array, jl_fnameh, opt_data, (jl_value_t *)jl_hess);

       return retptr;
   }

    void solve_base() { jl_call1(jl_solve_base_case, opt_data); }


  jl_value_t* solve_base(jl_value_t* prevsol, double *grad, double *hess) 
  {

     base_sol = (prevsol == nullptr? jl_call1(jl_solve_base_case, opt_data): solve_base_case_with_recourse(prevsol, grad, hess)); 
   //  jl_call1(jl_debug_base_case, base_sol);

     std::string fname=instance+"/solution.csv";
     jl_value_t*  jl_fname = jl_cstr_to_string(fname.c_str());

    // save solution     
     jl_call3(jl_save_solution, jl_fname, opt_data, base_sol);

     return base_sol;
   }

   void solve_contingency_prob(int i, jl_value_t* bsol) 
   {  
      cont_sol = jl_call3(jl_solve_contingency_pridec, opt_data, jl_box_int64(i+1), bsol); 

     std::string fname=instance+"/contingency_"+std::to_string(i+1) +".csv";
     jl_value_t*  jl_fname = jl_cstr_to_string(fname.c_str());

    // save solution     
     jl_call3(jl_save_cont_solution, jl_fname, opt_data, cont_sol);

    }

   void solve_contingency_recourse(int i, double& rval) 
   {  
      receive_solution();
      solve_contingency_prob(i, base_sol);  //cont_sol
      getCost(rval);
    }

    void getGradient(double* x_vec) { getGradient(cont_sol, x_vec); }

    void getGradient(jl_value_t* _base_sol, double* x_vec) 
    {
       jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);
       jl_array_t* jl_x = jl_ptr_to_array_1d(array_type, x_vec, getDim(), 0);

       jl_value_t *res = jl_call2(jl_getGradient, _base_sol, (jl_value_t*)jl_x);
       }

    double getObjective(jl_value_t* _base_sol) { return jl_unbox_float64(jl_call1(jl_getObjective, _base_sol)); }

    void getSolution(jl_value_t* _base_sol, double* x_vec) 
    {
       jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);
       jl_array_t* jl_x = jl_ptr_to_array_1d(array_type, x_vec, getDim(), 0);

       jl_value_t *res = jl_call2(jl_getSolution, _base_sol, (jl_value_t*)jl_x);
       }

    JL_Interface(const JL_Interface& jlobj): 
      opt_data(jlobj.opt_data), cont_sol(nullptr), base_sol(nullptr), send_buffer(nullptr),
      size_buffer(0), instance(jlobj.instance)
      { 
           init_MPI();  

       }
  
    int64_t getDim() const {   return jl_unbox_int64(jl_call1(jl_getDim, opt_data)); }

    int64_t number_of_contingencies() const { return jl_unbox_int64(jl_call1(jl_number_of_contingencies, opt_data)); }
    int64_t number_of_columns() const { return jl_unbox_int64(jl_call1(jl_number_of_columns, opt_data)); }
 
  jl_value_t* jl_array(uint8_t *, int );

   void* serialize_data(int &_size, jl_value_t *_data_ptr)  
   { 

      jl_value_t* data_bytes= jl_call1(jl_serialize_obj, _data_ptr); 

     _size =jl_unbox_int64(jl_call1(jl_get_data_size, data_bytes));

    //  uint8_t *_ptr = new uint8_t[size]; 
      uint8_t *_ptr = alloc_buffer(_size); 
      jl_call2(jl_get_data_bytes, data_bytes,jl_array(_ptr, _size));
    
     return _ptr;
   }
    jl_value_t* deserialize_data(uint8_t *_ptr, int _size) //: cont_sol(nullptr)  
    { 
   
        return jl_call1(jl_deserialize_obj, jl_array(_ptr, _size)); 
    }

    JL_Interface();

  void init_MPI();

   void send_MPI_data(jl_value_t*, int tag=0, bool block=true);
   jl_value_t* receive_MPI_data(int tag=0, bool block=true);
   void send_instance()  {  send_MPI_data(opt_data,77); }
   void receive_instance() { opt_data =receive_MPI_data(77); }

   void send_solution()  {  send_MPI_data(base_sol,99, false); }
   void receive_solution() 
   { 
      base_sol =receive_MPI_data(99, false); 
    }

  };

#endif

