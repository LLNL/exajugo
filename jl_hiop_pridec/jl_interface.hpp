
#ifndef JL_INTERFACE_HPP

#define JL_INTERFACE_HPP

#pragma once 

#include <string>
#include <iostream>
#include <vector>

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
extern jl_function_t* jl_debug_array;

extern jl_function_t* jl_save_solution;
extern jl_function_t* jl_save_cont_solution;
extern jl_function_t* jl_save_array;

extern jl_function_t* jl_struct_to_array_generic;
extern jl_function_t* jl_array_to_struct;

extern jl_function_t* jl_full_solution_dim;
extern jl_function_t* jl_define_array_lengths;

extern jl_function_t* jl_get_data_ptr;

/*
extern jl_value_t* opt_data;   // Julia object for optimization data
extern jl_value_t* base_sol;   // Julia object for base solution
extern jl_value_t* cont_sol;   // Julia object for contingency solution
extern jl_value_t* fieldsizes;
*/

void include_jl_functions();

#include <julia.h>
#include <iostream>

class JuliaRooted {
public:
    jl_value_t* julia_obj;
    jl_array_t* root_array;

    JuliaRooted(jl_value_t* obj) : julia_obj(obj) 
    {
        julia_obj = obj;
        jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_any_type, 1);

        root_array = (jl_array_t*)jl_alloc_array_1d(array_type, 1);

        jl_array_ptr_set(root_array, 0, obj);
    }

    ~JuliaRooted() {
        // No manual cleanup needed for root_array
    }

    jl_value_t* get() const {
        return julia_obj;
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

    // global variables in Julia are used to keep the following pointers valid
    // erialization/deserialization can be used too
    // jl_gc_preserve/jl_gc_unpreserve can also be used but the are not available 
    //     in the current julia module used

    jl_value_t* opt_data;   // Julia object for optimization data
    jl_value_t* base_sol;   // Julia object for base solution
    jl_value_t* cont_sol;   // Julia object for contingency solution
    jl_value_t* fieldsizes; // Used to reconstruct Basesolution struct from master solution
  //  void* tag;
  //  JuliaRooted* rooted_julia_obj = nullptr;

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

 
    // Get optimization data pointer
    jl_value_t* get_data_ptr() const 
      // { return rooted_julia_obj->get();}
       { return opt_data;}

    jl_value_t* get_cont_sol() const { return cont_sol;}

    jl_value_t* get_base_sol() const { return base_sol;}

    // Set optimization data pointer
    void set_data_ptr(jl_value_t* _ptr)
    {
        if (opt_data != nullptr)
        {
            std::cerr << "Data pointer already set! Rank: " << rank << std::endl;
            exit(EXIT_FAILURE);
        }
        opt_data = _ptr;

    //rooted_julia_obj = new JuliaRooted(opt_data);       
    // tag = jl_gc_preserve(opt_data, NULL);

    }

    void set_cont_sol(jl_value_t* ptr) { cont_sol = ptr;}

    void set_base_sol(jl_value_t* ptr) { base_sol = ptr;}

    jl_value_t* read_data() 
    {

       std::string exajugo_path = std::getenv("PATH_TO_EXAJUGO");
       std::string example_path = exajugo_path+"/examples/"+instance+"/";
       std::string userCommand= "rm -Rf "+instance+"; mkdir -p "+instance;
       //int result = 
       system(userCommand.c_str());

       std::string raw_file = (example_path+"case.raw"); 
       std::string rop_file = (example_path+"case.rop"); 
       std::string con_file = (example_path+"case.con");
      
       jl_value_t* jl_raw_file = jl_cstr_to_string(raw_file.c_str()); 
       jl_value_t* jl_rop_file = jl_cstr_to_string(rop_file.c_str()); 
       jl_value_t* jl_con_file = jl_cstr_to_string(con_file.c_str());

       return jl_call3(jl_load_ACOPF, jl_raw_file, jl_rop_file, jl_con_file);
       }

public:
  
    int get_max_iter() { return max_iter; }

    JL_Interface(const std::string& , const int _max_it=100);

    // Destructor
    ~JL_Interface() {  release_buffer(); 
   // jl_gc_enable(1); 
   /// jl_gc_unpreserve(tag);
   // delete rooted_julia_obj;
     }

    // Initialize MPI
    void init_MPI();
    
    // Send Julia object via MPI
    void send_MPI_data(jl_value_t* _dt_ptr, int tag = 0, bool block = true);
     
    // Receive Julia object via MPI
    jl_value_t* receive_MPI_data(int tag = 0, bool block = true);

    // Send base solution
    void send_solution() {  send_MPI_data(get_base_sol(), 99, false);  }

    // Receive base solution
    void receive_solution()  { set_base_sol(receive_MPI_data(99));  }

    void getCost(double& rval) { rval =  jl_unbox_float64(jl_call1(jl_getCost, get_cont_sol())); }

    void solve_contingency_recourse(int i, double& rval) 
    {  
        receive_solution();
        solve_contingency_prob(i);  //cont_sol
        getCost(rval);

    }

    // Solve contingency problem
    void solve_contingency_prob(int i)
    {

       set_cont_sol(jl_call3(jl_solve_contingency_pridec, get_data_ptr(), jl_box_int64(i+1), get_base_sol()));

       std::string fsolcont="contingency_"+std::to_string(i+1) +".csv";

       save_jl_array(jl_save_cont_solution, fsolcont, get_cont_sol());

    }
    
    jl_value_t* jl_array(double *_ptr, int _size);


    // Get gradient
    void getGradient(double* x_vec)
    {
        jl_call2(jl_getGradient, get_cont_sol(), jl_array(x_vec, getDim()));
    }

    // Get objective value
    double getObjective()
    {
        return jl_unbox_float64(jl_call1(jl_getObjective, get_base_sol()));
    }

    // Get solution
    void getSolution(double* x_vec)
    {
       jl_value_t* jl_x =  jl_array(x_vec, getDim());
       jl_call2(jl_getSolution, get_base_sol(), (jl_value_t*)jl_x);
    }

    int64_t number_of_contingencies() const { return jl_unbox_int64(jl_call1(jl_number_of_contingencies, get_data_ptr())); }

    int64_t number_of_columns() const 
    {
         return jl_unbox_int64(jl_call1(jl_number_of_columns, get_data_ptr()));
    }

    // Get optimization problem dimension
    int64_t getDim() const
    {
       return jl_unbox_int64(jl_call1(jl_getDim, get_data_ptr()));
    }
   
   // Function to save a Julia array to a CSV file
    void save_jl_array(jl_function_t* jl_save,
                       const std::string& filename, jl_value_t* array_ptr) 
    {
       std::string fullpath = instance + "/" + filename;
       jl_value_t* jl_fname = jl_cstr_to_string(fullpath.c_str());
       jl_call3(jl_save, jl_fname, get_data_ptr(), array_ptr);
    }

    // Solve base optimization problem
    void solve_base_case_with_recourse(double *grad, double *hess) 
    {
       jl_value_t* jl_grad = jl_array(grad, getDim());
       jl_value_t* jl_hess = jl_array(hess, getDim());

// this is necessary to root pointers tp protect from Julia GC
       JL_GC_PUSH2(&jl_grad, &jl_hess);

//       jl_value_t* ptr_rderivatives = jl_call2(jl_get_recourse_derivatives, jl_grad, jl_hess);
       jl_value_t* ptr_rderivatives = 
          jl_call3(jl_get_recourse_derivatives, jl_grad, jl_hess, jl_box_int64(getDim()));

       set_base_sol(jl_call3(jl_solve_base_case_recourse, get_data_ptr(), get_base_sol(), ptr_rderivatives));

       save_jl_array(jl_save_array, "gradient.csv", (jl_value_t*)jl_grad);
       save_jl_array(jl_save_array, "hessian.csv", (jl_value_t*)jl_hess);

       JL_GC_POP();
     }

     void solve_base() 
     {
        set_base_sol(jl_call1(jl_solve_base_case, get_data_ptr())); 
     }


     void solve_base(double *grad, double *hess) 
     {
        if (get_base_sol() == nullptr)
           solve_base();
        else
           solve_base_case_with_recourse(grad, hess);

       save_jl_array(jl_save_solution, "solution.csv", get_base_sol());
     }

};

#endif

