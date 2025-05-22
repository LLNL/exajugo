
#ifndef JL_INTERFACE_HPP

#define JL_INTERFACE_HPP

#include <string>
#include <iostream>
#include <vector>

#include <julia.h>

#include <mutex>
extern std::mutex julia_mutex;

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

extern jl_function_t* jl_serialize_object;
extern jl_function_t* jl_deserialize_object;


void include_jl_functions();

// Define logging macros for better debugging
//#define LOG_DEBUG(msg) if (debug_enabled) { std::cout << "[DEBUG] " << msg << std::endl; }
//#define LOG_INFO(msg) if (info_enabled) { std::cout << "[INFO] " << msg << std::endl; }

// JL_Interface class definition
class JL_Interface
{
private:
    int source_process = 0; // Default source process for MPI communication
    int rank;               // MPI rank
    int nproc;              // Number of MPI processes
    int size_buffer;        // Size of the buffer
    void* data_ptr;         // Pointer to data
    uint8_t* send_buffer;   // Buffer for MPI communication
    std::string instance;   // Instance name
//    bool debug_enabled = true; // Debug flag
 //   bool info_enabled = true;  // Info flag

protected:
    jl_value_t* opt_data;   // Julia object for optimization data
    jl_value_t* base_sol;   // Julia object for base solution
    jl_value_t* cont_sol;   // Julia object for contingency solution

    // Helper method to allocate buffer
    uint8_t* alloc_buffer(int size)
    {
       release_buffer();
        uint8_t* buffer = (uint8_t*)malloc(size);
        if (!buffer)
        {
            std::cerr << "Failed to allocate buffer!" << std::endl;
            exit(EXIT_FAILURE);
        }
        return buffer;
    }

       // Release the buffer after sending
 
    // Helper method to release buffer
    void release_buffer()
    {
        if (send_buffer)
        {
            free(send_buffer);
            send_buffer = nullptr;
        }
    }

    // Set optimization data pointer
    void set_data_ptr(jl_value_t* _ptr)
    {
        if (opt_data != nullptr)
        {
            std::cerr << "Data pointer already set! Rank: " << rank << std::endl;
            exit(EXIT_FAILURE);
        }
        opt_data = _ptr;

        // Assign opt_data to the array
      //  jl_array_ptr_set(gc_protected_array, 0, opt_data);
         //  JL_GC_PUSH1(opt_data);

    }

    // Get optimization data pointer
    jl_value_t* get_data_ptr() const { return opt_data; }

    // Get base solution pointer
    jl_value_t* get_base_sol() const { return base_sol; }

    // Set base solution pointer
    void set_base_sol(jl_value_t* ptr)
    {
        if (ptr == nullptr)
        {
            std::cerr << "Error: Null pointer received for base solution!" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (base_sol != nullptr)
        {
          //  JL_GC_POP();
            //LOG_DEBUG("Object popped!");
        }

        base_sol = ptr;

        if (base_sol != nullptr)
        {
 //           JL_GC_PUSH1(base_sol);
           // LOG_DEBUG("Object pushed! Rank: " << rank << " :: " << base_sol);
        }
    }

    // Serialize Julia object into a buffer
    void* serialize_data(int& size, jl_value_t* _data_ptr)
    {
        jl_value_t* data_bytes = jl_call1(jl_serialize_obj, _data_ptr);
        size = jl_unbox_int64(jl_call1(jl_get_data_size, data_bytes));

        uint8_t* buffer = alloc_buffer(size);
        jl_call2(jl_get_data_bytes, data_bytes, jl_array(buffer, size));

        return buffer;
    }

    // Deserialize buffer into Julia object
    jl_value_t* deserialize_data(uint8_t* buffer, int size)
    {
        return jl_call1(jl_deserialize_obj, jl_array(buffer, size));
    }

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
    // Constructor
    JL_Interface() : size_buffer(0), data_ptr(nullptr), send_buffer(nullptr),
                                   opt_data(nullptr), base_sol(nullptr), cont_sol(nullptr)
    {
        instance = "9bus"; // Default instance name
        include_jl_functions(); // Load Julia functions

        init_MPI(); // Initialize MPI

        if (rank == source_process)
        {
            // Read data and set the optimization data pointer
            set_data_ptr(read_data());
            //LOG_INFO("Data read and set: " << get_data_ptr());

            // Send the instance to other processes
            send_instance();

        }
        else
        {
            // Receive the instance from the source process
            receive_instance();
        }
    }

    // Destructor
    ~JL_Interface() { release_buffer();
   // jl_gc_preserve_pop(1); 
}

    // Initialize MPI
    void init_MPI();
    // Send Julia object via MPI
    void send_MPI_data(jl_value_t* _dt_ptr, int tag = 0, bool block = true);
    // Receive Julia object via MPI
    jl_value_t* receive_MPI_data(int tag = 0, bool block = true);

    // Send optimization instance
    void send_instance() { send_MPI_data(get_data_ptr(), 77); }

    // Receive optimization instance
    void receive_instance()
    {
        set_data_ptr(receive_MPI_data(77));

    }

    // Send base solution
    void send_solution() { send_MPI_data(get_base_sol(), 99, false); }

    // Receive base solution
    void receive_solution()
    {

        set_base_sol(receive_MPI_data(99, false));

        //LOG_INFO("Base solution received.");
    }
    void getCost(double& rval) { rval =  jl_unbox_float64(jl_call1(jl_getCost, cont_sol)); }


   void solve_contingency_recourse(int i, double& rval) 
   {  
      receive_solution();
      solve_contingency_prob(i, base_sol);  //cont_sol
      getCost(rval);
    }

    // Solve contingency problem
    void solve_contingency_prob(int i, jl_value_t* bsol)
    {
        //LOG_INFO("Solving contingency problem: Index: " << i + 1 << " Rank: " << rank);

        int cindex = i + 1;
        jl_value_t* jl_cindex = jl_box_int64(cindex);

        cont_sol = jl_call3(jl_solve_contingency_pridec, get_data_ptr(), jl_cindex, bsol);


     std::string fname=instance+"/contingency_"+std::to_string(i+1) +".csv";
     jl_value_t*  jl_fname = jl_cstr_to_string(fname.c_str());

    // save solution     
     jl_call3(jl_save_cont_solution, jl_fname, get_data_ptr(), cont_sol);


        //LOG_INFO("Contingency problem solved: Index: " << i + 1);
    }
    jl_value_t* jl_array(uint8_t *_ptr, int _size);
  jl_value_t* jl_float_array(double* x_vec);


    // Get gradient
    void getGradient(double* x_vec)
    {
        jl_call2(jl_getGradient, cont_sol, jl_float_array(x_vec));
    }

    // Get objective value
    double getObjective()
    {
        return jl_unbox_float64(jl_call1(jl_getObjective, base_sol));
    }

    // Get solution
    void getSolution(double* x_vec)
    {
        jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);
        jl_array_t* jl_x = jl_ptr_to_array_1d(array_type, x_vec, getDim(), 0);

        jl_call2(jl_getSolution, get_base_sol(), (jl_value_t*)jl_x);
    }

    int64_t number_of_contingencies() const { return jl_unbox_int64(jl_call1(jl_number_of_contingencies, get_data_ptr())); }
    int64_t number_of_columns() const { return jl_unbox_int64(jl_call1(jl_number_of_columns, get_data_ptr())); }


    // Get optimization problem dimension
    int64_t getDim() const
    {
        return jl_unbox_int64(jl_call1(jl_getDim, get_data_ptr()));
    }
        // Solve base optimization problem
void solve_base_case_with_recourse(double *grad, double *hess) 
    {
       jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);
       jl_array_t* jl_grad = jl_ptr_to_array_1d(array_type, grad, getDim(), 0);

       // this line is not needed: can use the same type above 'array_type'
       jl_value_t* array_typeh = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);
       jl_array_t* jl_hess = jl_ptr_to_array_1d(array_typeh, hess, getDim(), 0);

       jl_value_t* ptr_rderivatives = jl_call2(jl_get_recourse_derivatives, (jl_value_t *)jl_grad, (jl_value_t *)jl_hess);

       set_base_sol(jl_call3(jl_solve_base_case_recourse, get_data_ptr(), get_base_sol(), ptr_rderivatives));
 
       // save gradient
       std::string fname=instance+"/gradient.csv";
       jl_value_t*  jl_fname = jl_cstr_to_string(fname.c_str());
       jl_call3(jl_save_array, jl_fname, get_data_ptr(), (jl_value_t *)jl_grad);

       // save hessian
       std::string fnameh=instance+"/hessian.csv";
       jl_value_t*  jl_fnameh = jl_cstr_to_string(fnameh.c_str());
       jl_call3(jl_save_array, jl_fnameh, get_data_ptr(), (jl_value_t *)jl_hess);

   }

    void solve_base() 
    {
         set_base_sol(jl_call1(jl_solve_base_case, get_data_ptr())); 
    }


      jl_value_t* solve_base(double *grad, double *hess) 
  {
     if (get_base_sol() == nullptr)
       solve_base();
     else
        solve_base_case_with_recourse(grad, hess);

     std::string fname=instance+"/solution.csv";
     jl_value_t*  jl_fname = jl_cstr_to_string(fname.c_str());

    // save solution     
     jl_call3(jl_save_solution, jl_fname, get_data_ptr(), base_sol);

     return base_sol;
   }

};

#endif

