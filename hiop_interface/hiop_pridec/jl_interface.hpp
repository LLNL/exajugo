
#ifndef JL_INTERFACE_HPP

#define JL_INTERFACE_HPP

#pragma once 

#include <julia.h>

#include <string>
#include <iostream>
#include <vector>

#include <memory>


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
extern jl_function_t* jl_debug_array;

extern jl_function_t* jl_save_solution;
extern jl_function_t* jl_save_cont_solution;
extern jl_function_t* jl_save_array;

extern jl_function_t* jl_struct_to_array_generic;
extern jl_function_t* jl_array_to_struct;

extern jl_function_t* jl_full_solution_dim;
extern jl_function_t* jl_define_array_lengths;

extern jl_function_t* jl_get_data_ptr;
extern jl_function_t* jl_hold_pointer;


void include_jl_functions();

// Pointer wrapper class

class JL_Pointer
{
    jl_value_t* ptr;  

    int getId() const { return static_cast<int>(reinterpret_cast<uintptr_t>(this)); }

public:

    JL_Pointer(jl_value_t *_ptr=nullptr) { set(_ptr); }

    jl_value_t* get() const { return ptr; }

    void set(jl_value_t *_ptr) 
    { 
        if ((ptr = _ptr) == nullptr) return;

        jl_call2(jl_hold_pointer, get(), jl_box_int64(getId()));

    }

};


class JL_Interface
{
private:
    int max_iter;
    int source_process = 0; // Default source process for MPI communication
    int rank;               // MPI rank
    int nproc;              // Number of MPI processes
    int size_buffer;
    double* data_buffer;   // Buffer for MPI communication
    std::string instance;   // Instance name

protected:

    JL_Pointer opt_data;   // Julia object for optimization data
    JL_Pointer base_sol;   // Julia object for base solution
    JL_Pointer cont_sol;   // Julia object for contingency solution
    JL_Pointer fieldsizes; // Used to reconstruct Basesolution struct from master solution

    double* alloc_buffer(int size)
    {
       release_buffer();
        double* buffer = new double [size]; 
        if (!buffer)
        {
            std::cerr << "Failed to allocate buffer!" << std::endl;
            exit(EXIT_FAILURE);
        }
        return buffer;
    }
 
    // Helper method to release buffer
    void release_buffer()
    {
        if (data_buffer)
        {
            delete[]data_buffer; 
            data_buffer = nullptr;
        }
    }

    jl_value_t* read_data() 
    {

       std::string exajugo_path = std::getenv("PATH_TO_EXAJUGO");
       std::string example_path = exajugo_path+"/examples/"+instance+"/";
       std::string userCommand= "rm -Rf "+instance+"; mkdir -p "+instance;
       system(userCommand.c_str());

       std::string raw_file = (example_path+"case.raw"); 
       std::string rop_file = (example_path+"case.rop"); 
       std::string con_file = (example_path+"case.con");
      
       jl_value_t* jl_raw_file = jl_cstr_to_string(raw_file.c_str()); 
       jl_value_t* jl_rop_file = jl_cstr_to_string(rop_file.c_str()); 
       jl_value_t* jl_con_file = jl_cstr_to_string(con_file.c_str());

       return jl_call3(jl_load_ACOPF, jl_raw_file, jl_rop_file, jl_con_file);
       }

    jl_value_t* get_field_data()
    {
        return jl_call1(jl_define_array_lengths, opt_data.get());
    }

public:
  
    int get_max_iter() { return max_iter; }

    JL_Interface(const std::string& , const int _max_it=100);

    // Destructor
    ~JL_Interface() {  release_buffer(); 
     }

    // Initialize MPI
    void init_MPI();
    
    // Send Julia object via MPI
    void send_MPI_data(jl_value_t* _dt_ptr, int tag = 0, bool block = true);
     
    // Receive Julia object via MPI
    jl_value_t* receive_MPI_data(int tag = 0, bool block = true);

    // Send base solution
    void send_solution() {  send_MPI_data(base_sol.get(), 99, false);  }

    // Receive base solution
    void receive_solution()  { base_sol.set(receive_MPI_data(99));  }

    void getCost(double& rval) { rval =  jl_unbox_float64(jl_call1(jl_getCost, cont_sol.get())); }

    void solve_contingency_recourse(int i, double& rval) 
    {  
        receive_solution();
        solve_contingency_prob(i);  //cont_sol
        getCost(rval);
    }

    // Solve contingency problem
    void solve_contingency_prob(int i)
    {

       cont_sol.set(jl_call3(jl_solve_contingency_pridec, opt_data.get(), jl_box_int64(i+1), base_sol.get()));

       std::string fsolcont="contingency_"+std::to_string(i+1) +".csv";

       save_jl_array(jl_save_cont_solution, fsolcont, cont_sol.get());

    }
    
    jl_value_t* jl_array(double *_ptr, int _size);


    // Get gradient
    void getGradient(double* x_vec)
    {
        jl_call2(jl_getGradient, cont_sol.get(), jl_array(x_vec, getDim()));
    }

    // Get objective value
    double getObjective()
    {
        return jl_unbox_float64(jl_call1(jl_getObjective, base_sol.get()));
    }

    // Get solution
    void getSolution(double* x_vec)
    {
       jl_value_t* jl_x =  jl_array(x_vec, getDim());
       jl_call2(jl_getSolution, base_sol.get(), (jl_value_t*)jl_x);
    }

    int64_t number_of_contingencies() const { return jl_unbox_int64(jl_call1(jl_number_of_contingencies, opt_data.get())); }

    int64_t number_of_columns() const 
    {
         return jl_unbox_int64(jl_call1(jl_number_of_columns, opt_data.get()));
    }

    // Get optimization problem dimension
    int64_t getDim() const
    {
       return jl_unbox_int64(jl_call1(jl_getDim, opt_data.get()));
    }
   
   // Function to save a Julia array to a CSV file
    void save_jl_array(jl_function_t* jl_save,
                       const std::string& filename, jl_value_t* array_ptr) 
    {
       std::string fullpath = instance + "/" + filename;
       jl_value_t* jl_fname = jl_cstr_to_string(fullpath.c_str());
       jl_call3(jl_save, jl_fname, opt_data.get(), array_ptr);
    }

    // Solve base optimization problem
    void solve_base_case_with_recourse(double *grad, double *hess) 
    {
       jl_value_t* jl_grad = jl_array(grad, getDim());
       jl_value_t* jl_hess = jl_array(hess, getDim());

// this is necessary to root pointers tp protect from Julia GC
       JL_GC_PUSH2(&jl_grad, &jl_hess);

       jl_value_t* ptr_rderivatives = 
          jl_call3(jl_get_recourse_derivatives, jl_grad, jl_hess, jl_box_int64(getDim()));

       base_sol.set(jl_call3(jl_solve_base_case_recourse, opt_data.get(), base_sol.get(), ptr_rderivatives));

       save_jl_array(jl_save_array, "gradient.csv", (jl_value_t*)jl_grad);
       save_jl_array(jl_save_array, "hessian.csv", (jl_value_t*)jl_hess);

       JL_GC_POP();
     }

     void solve_base() 
     {
        base_sol.set(jl_call1(jl_solve_base_case, opt_data.get())); 
     }


     void solve_base(double *grad, double *hess) 
     {
        if (base_sol.get() == nullptr)
           solve_base();
        else
           solve_base_case_with_recourse(grad, hess);

       save_jl_array(jl_save_solution, "solution.csv", base_sol.get());
     }

};

#endif

