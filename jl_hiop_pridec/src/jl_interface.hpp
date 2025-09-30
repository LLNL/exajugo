
#ifndef JL_INTERFACE_HPP

#define JL_INTERFACE_HPP

#pragma once 

#include <julia.h>

#include <string>
#include <iostream>
#include <vector>

#include <memory>

#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <sstream>

#include <filesystem>
//namespace fs = std::filesystem;
extern const char preferred_separator;

extern jl_function_t* jl_load_ACOPF;
extern jl_function_t* jl_load_ACOPF_instance;
extern jl_function_t* jl_copy_ACOPF;

extern jl_function_t* jl_number_of_contingencies;
extern jl_function_t* jl_number_of_columns;

extern jl_function_t* jl_solve_base_case ;

extern jl_function_t* jl_solve_base_case_recourse;
extern jl_function_t* jl_solve_base_case_recourse_sparse;
extern jl_function_t* jl_get_recourse_derivatives;
extern jl_function_t* jl_get_recourse_sparse;
extern jl_function_t* jl_get_sparse_matrix_index_wrap;
extern jl_function_t* jl_get_sparse_matrix_wrap;

extern jl_function_t* jl_getModel;
extern jl_function_t* jl_getDim;
extern jl_function_t* jl_getObjective;
extern jl_function_t* jl_getSolution;

extern jl_function_t* jl_solve_contingency_pridec;
extern jl_function_t* jl_getCost;
extern jl_function_t* jl_getGradient;
extern jl_function_t* jl_gmultiply_array;
    
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
extern jl_function_t* jl_save_sparse_matrix;

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
    std::string outputDir;


std::string buildOutputPath(const std::string& fileName, int cont_id) {
    std::ostringstream dir;
    dir << outputDir << "/" << cont_id;
    mkdir(dir.str().c_str(), 0777); // creates directory if it doesn't exist

    std::ostringstream filePath;
    filePath << dir.str() << "/" << fileName << ".csv";
    return filePath.str();
}


protected:

    JL_Pointer opt_data;   // Julia object for optimization data
    JL_Pointer base_sol;   // Julia object for base solution
    JL_Pointer cont_sol;   // Julia object for contingency solution
    JL_Pointer fieldsizes; // Used to reconstruct Basesolution struct from master solution

    double grad_multiplier;
    double hess_multiplier;
    int iter;

    JL_Pointer sparse_index;
    JL_Pointer sparse_hessian;

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

    const char* env_path = std::getenv("PATH_TO_INSTANCES");
    assert(env_path != nullptr); // This will catch unset env variable
    std::string exajugo_path(env_path);

       assert(!exajugo_path.empty());
       assert(!instance.empty());
       assert(preferred_separator != '\0');

       //std::string example_path = exajugo_path+instance+fs::path::preferred_separator;
       std::string example_path = exajugo_path+instance+preferred_separator;

       jl_value_t* jl_opt_instance = jl_cstr_to_string(instance.c_str());

       return jl_call1(jl_load_ACOPF_instance, jl_opt_instance);
       }

    jl_value_t* get_field_data()
    {
        return jl_call1(jl_define_array_lengths, opt_data.get());
    }

public:
  
    void set_grad_multiplier(double _val) { grad_multiplier =_val; }
    void set_hess_multiplier(double _val) 
    { 
        hess_multiplier =_val; 
        std::cout<<" hess_multiplier set to "<<hess_multiplier<<std::endl;
    }

    int get_max_iter() { return max_iter; }

    JL_Interface(const std::string&, const std::string&, const int _max_it=100);

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
    void send_solution() {  if (nproc>1) send_MPI_data(base_sol.get(), 99);  }

    // Receive base solution
    void receive_solution()  { if (nproc>1) base_sol.set(receive_MPI_data(99));  }

    void getCost(double& rval) { rval =  jl_unbox_float64(jl_call1(jl_getCost, cont_sol.get())); }

    void solve_contingency_recourse(int i, double& rval) 
    {  
        if (i<nproc-1)
           receive_solution(); 

        solve_contingency_prob(i);  //cont_sol
        getCost(rval);

    }

    // Solve contingency problem
    void solve_contingency_prob(int i)
    {
       int cont_id = i+1;

       cont_sol.set(jl_call3(jl_solve_contingency_pridec, opt_data.get(), jl_box_int64(cont_id), base_sol.get()));
       save_jl_array(jl_save_cont_solution, "solution_"+std::to_string(cont_id), cont_sol.get(), cont_id);

    }
    
    jl_value_t* jl_array(double *_ptr, int _size);
    jl_value_t* jl_array(int64_t *_ptr, int _size);
    jl_value_t* jl_array_mult(double *_ptr, int _size, double _mult);

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
                       const std::string& filename, jl_value_t* array_ptr, int cont_id=0) 
    {
       std::string fullpath = buildOutputPath(filename, cont_id); 

       jl_value_t* jl_fname = jl_cstr_to_string(fullpath.c_str());
       jl_call3(jl_save, jl_fname, opt_data.get(), array_ptr);
    }

   // Function to save a Julia array to a CSV file
    void save_jl_array(const std::string& filename, jl_value_t* array_ptr, int _iter) 
    {
       std::string fullpath = buildOutputPath(filename, 0); 

       jl_value_t* jl_fname = jl_cstr_to_string(fullpath.c_str());
       jl_call3(jl_save_sparse_matrix, jl_fname, array_ptr, jl_box_int64(_iter));
    }

    void sparse_matrix_wrap(int64_t *_rows, int64_t *_cols, double *_vals, int64_t _nelements)
    {
       jl_value_t* jl_rows = jl_array(_rows, _nelements); 
       jl_value_t* jl_cols = jl_array(_cols, _nelements); 
       jl_value_t* jl_vals = jl_array(_vals, _nelements);

      // JL_GC_PUSH1(&jl_vals);

       sparse_index.set(jl_call3(jl_get_sparse_matrix_index_wrap, jl_rows, jl_cols, jl_box_int64(_nelements)));

       sparse_hessian.set(jl_call3(jl_get_sparse_matrix_wrap, sparse_index.get(), jl_vals, jl_box_int64(_nelements)));
      // JL_GC_POP();

    }

    void denseToSparse(const double* hess, size_t nrows, int64_t*& _rows, int64_t*& _cols, double*& _vals, int64_t& nnz) 
    {
        nnz = 0;
        // First pass: count nonzeros
        for (size_t i = 0; i < nrows; ++i)
            ++nnz;

        // Allocate arrays
        _rows = (int64_t*) malloc(nnz * sizeof(int64_t));
        _cols = (int64_t*) malloc(nnz * sizeof(int64_t));
        _vals = (double*)  malloc(nnz * sizeof(double));
        // Second pass: fill arrays
        size_t idx = 0;
        for (size_t i = 0; i < nrows; ++i) 
        {
                double val = hess[i];
                if (val != 0.0) 
                {
                    _rows[idx] = static_cast<int64_t>(i);
                    _cols[idx] = static_cast<int64_t>(i);
                    _vals[idx] = val;
                    ++idx;
                }
            }
    }

    void free_arrays(int64_t*& _rows, int64_t*& _cols, double*& _vals) 
    {
         free(_rows); _rows=nullptr;
         free(_cols); _cols=nullptr;
         free(_vals); _vals=nullptr;
    }

    // Solve base optimization problem
    void test_solve_base_case_with_recourse(double *grad, double *hess) 
    {
       jl_value_t* jl_grad = jl_array_mult(grad, getDim(), grad_multiplier);
     //  jl_value_t* jl_hess = jl_array_mult(hess, getDim(), hess_multiplier);

       int64_t* _rows;
       int64_t* _cols;
       double* _vals;
       int64_t nnz=getDim();
       iter+=1;

      // convert hessian to sparse for testing
       denseToSparse(hess, getDim(), _rows, _cols, _vals, nnz);

       sparse_matrix_wrap(_rows, _cols, _vals, nnz);

// this is necessary to root pointers tp protect from Julia GC
       JL_GC_PUSH1(&jl_grad);

       jl_value_t* ptr_rderivatives = 
          jl_call3(jl_get_recourse_sparse, jl_grad, sparse_hessian.get(), jl_box_int64(getDim()));

       base_sol.set(jl_call3(jl_solve_base_case_recourse_sparse, opt_data.get(), base_sol.get(), ptr_rderivatives));

       save_jl_array(jl_save_array, "gradient", (jl_value_t*)jl_grad);
       save_jl_array("hessian", sparse_hessian.get(), iter);

       JL_GC_POP();

       free_arrays(_rows, _cols, _vals);

     }

    // Solve base optimization problem
    void solve_base_case_with_recourse(double *grad, double *hess) 
    {
       jl_value_t* jl_grad = jl_array_mult(grad, getDim(), grad_multiplier);
       jl_value_t* jl_hess = jl_array_mult(hess, getDim(), hess_multiplier);

       // this is necessary to protect pointers from Julia GC
       JL_GC_PUSH2(&jl_grad, &jl_hess);

       jl_value_t* ptr_rderivatives = 
          jl_call3(jl_get_recourse_derivatives, jl_grad, jl_hess, jl_box_int64(getDim()));

       base_sol.set(jl_call3(jl_solve_base_case_recourse, opt_data.get(), base_sol.get(), ptr_rderivatives));

       save_jl_array(jl_save_array, "gradient", (jl_value_t*)jl_grad);
       save_jl_array(jl_save_array, "hessian", (jl_value_t*)jl_hess);

       JL_GC_POP();
     }

     void solve_base() 
     {
        base_sol.set(jl_call1(jl_solve_base_case, opt_data.get())); 

     }
    
     bool success() { return base_sol.get() != nullptr; }

     void solve_base(double *grad, double *hess) 
     {
        if (base_sol.get() == nullptr)
           solve_base();
        else
           solve_base_case_with_recourse(grad, hess);
           //test_solve_base_case_with_recourse(grad, hess);
       

       std::cout<<" --- BASE CASE SOLVED! --- \n";
       try {

         jl_value_t* solptr = base_sol.get();

         if (solptr == nullptr)
         {
             std::cout<<" --- Solution could not be retrieved! ---\n"; 
             return;
          }
          std::cout<<" --- sol pointer: "<< solptr <<" ---\n"; 
          save_jl_array(jl_save_solution, "solution", solptr);
        } 
        catch (...) { // Catch-all handler
          std::cerr << "Error saving solution!" << std::endl;
         }
       std::cout<<" --- SOLUTION saved! ---\n";

     }

};

#endif

