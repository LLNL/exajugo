

#include <vector>
#include <cstring>  //for memcpy

#include "jl_interface.hpp"

#define HIOP_USE_MPI

#ifdef HIOP_USE_MPI
#include "mpi.h"

#else
#define MPI_Comm int
#define MPI_COMM_WORLD 0
#endif


jl_function_t* jl_load_ACOPF;
jl_function_t* jl_copy_ACOPF;

jl_function_t* jl_number_of_contingencies;
jl_function_t* jl_number_of_columns;

jl_function_t* jl_solve_base_case;

jl_function_t* jl_solve_base_case_recourse;
jl_function_t* jl_get_recourse_derivatives;

jl_function_t* jl_getModel;
jl_function_t* jl_getDim;
jl_function_t* jl_getObjective;
jl_function_t* jl_getSolution;

jl_function_t* jl_solve_contingency_pridec;
jl_function_t* jl_getCost;
jl_function_t* jl_getGradient;
    
jl_function_t*  deepcopy_func;
jl_function_t*  jl_serialize_obj;
jl_function_t*  jl_deserialize_obj;

jl_function_t*  jl_get_data_size;
jl_function_t*  jl_get_data_bytes;

jl_function_t* jl_debug_base_case;
jl_function_t* jl_debug_array;

jl_function_t* jl_save_solution;
jl_function_t* jl_save_cont_solution;
jl_function_t* jl_save_array;

jl_function_t* jl_struct_to_array_generic;
jl_function_t* jl_array_to_struct;

jl_function_t* jl_full_solution_dim;
jl_function_t* jl_define_array_lengths;

jl_function_t* jl_get_data_ptr;

/*
jl_value_t* opt_data=nullptr;   // Julia object for optimization data
jl_value_t* base_sol=nullptr;   // Julia object for base solution
jl_value_t* cont_sol=nullptr;   // Julia object for contingency solution
jl_value_t* fieldsizes=nullptr;

*/

void include_jl_functions()
{
    const char* julia_file_path = "hiop.jl";

    std::string command = "include(\"" + std::string(julia_file_path) + "\")"; 

    jl_eval_string(command.c_str());
    
    jl_load_ACOPF = jl_get_function(jl_main_module, "load_ACOPF");
    jl_copy_ACOPF = jl_get_function(jl_main_module, "copy_ACOPF");

    jl_number_of_contingencies = jl_get_function(jl_main_module, "number_of_contingencies");
    jl_number_of_columns = jl_get_function(jl_main_module, "number_of_columns");
    jl_solve_base_case = jl_get_function(jl_main_module, "solve_base_case");

    jl_solve_base_case_recourse = jl_get_function(jl_main_module, "solve_base_case_recourse");
    jl_get_recourse_derivatives = jl_get_function(jl_main_module, "get_recourse_derivatives");

    jl_getModel = jl_get_function(jl_main_module, "getModel");
    jl_getDim = jl_get_function(jl_main_module, "getDim");
    jl_getObjective = jl_get_function(jl_main_module, "getObjective");
    jl_getSolution = jl_get_function(jl_main_module, "getSolution");

    jl_solve_contingency_pridec = jl_get_function(jl_main_module, "solve_contingency_pridec");
    jl_getCost = jl_get_function(jl_main_module, "getCost");
    jl_getGradient = jl_get_function(jl_main_module, "getGradient");
    
    deepcopy_func = jl_get_function(jl_base_module, "deepcopy");

    jl_serialize_obj = jl_get_function(jl_main_module, "serialize_obj");
    jl_deserialize_obj = jl_get_function(jl_main_module, "deserialize_obj");

    jl_get_data_size = jl_get_function(jl_main_module, "get_data_size");
    jl_get_data_bytes = jl_get_function(jl_main_module, "get_data_bytes");
    jl_debug_base_case = jl_get_function(jl_main_module, "debug_base_case");
    jl_debug_array = jl_get_function(jl_main_module, "debug_array");

    jl_save_solution = jl_get_function(jl_main_module, "save_solution");
    jl_save_cont_solution = jl_get_function(jl_main_module, "save_cont_solution");
    jl_save_array = jl_get_function(jl_main_module, "save_array");

    jl_struct_to_array_generic = jl_get_function(jl_main_module, "struct_to_array_generic!");
    jl_array_to_struct = jl_get_function(jl_main_module, "array_to_struct");

    jl_full_solution_dim = jl_get_function(jl_main_module, "full_solution_dim");
    jl_define_array_lengths = jl_get_function(jl_main_module, "define_array_lengths");

    jl_get_data_ptr = jl_get_function(jl_main_module, "get_data_ptr");
}


jl_value_t* JL_Interface::jl_array(double *_ptr, int _size)
{
    jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);
    return (jl_value_t*) jl_ptr_to_array_1d(array_type, _ptr, _size, 0);
} 


// Constructor
JL_Interface::JL_Interface(const std::string& _inst, const int _max_it) 
   : max_iter(_max_it), size_buffer(0), data_buffer(nullptr), instance(_inst), 
   opt_data(nullptr), base_sol(nullptr),  cont_sol(nullptr), fieldsizes(nullptr)
{
   // jl_gc_enable(0);
    include_jl_functions(); // Load Julia functions

    init_MPI(); // Initialize MPI
    set_data_ptr(read_data());
 
    fieldsizes=jl_call1(jl_define_array_lengths, get_data_ptr()); 
    size_buffer = jl_unbox_int64(jl_call1(jl_full_solution_dim, fieldsizes));
}

// Send Julia object via MPI
void JL_Interface::send_MPI_data(jl_value_t* _dt_ptr, int tag, bool block)
{
    data_buffer = alloc_buffer(size_buffer);
    jl_value_t* jl_data_buffer_s = jl_array(data_buffer, size_buffer);

    jl_call2(jl_struct_to_array_generic, _dt_ptr, jl_data_buffer_s);

    // Send the data using MPI
    for (int i = 1; i < nproc; i++)
    {
        if (block)
        {
            MPI_Send(data_buffer, size_buffer, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Request request;
            MPI_Isend(data_buffer, size_buffer, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &request);
            // MPI_Wait(&request, MPI_STATUS_IGNORE);
            // int flag;
            // MPI_Test(&request, &flag, MPI_STATUS_IGNORE);
        }
    }
}

void JL_Interface::init_MPI()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
}

// Receive master solution via MPI
jl_value_t* JL_Interface::receive_MPI_data(int tag, bool block)
{
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);

    int data_size;
    MPI_Get_count(&status, MPI_DOUBLE, &data_size);

    data_buffer = alloc_buffer(data_size);
    MPI_Recv(data_buffer, data_size, MPI_DOUBLE, status.MPI_SOURCE, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    jl_value_t* jl_data_buffer= jl_array(data_buffer, data_size);
    jl_value_t* received_data = jl_call3(jl_array_to_struct, get_data_ptr(), jl_data_buffer, fieldsizes);

    return received_data;
}



