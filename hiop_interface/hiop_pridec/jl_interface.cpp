

#include <vector>
#include <cstring>  //for memcpy

#include "jl_interface.hpp"

//#include <julia.h>

std::mutex julia_mutex;

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
jl_function_t* jl_save_solution;
jl_function_t* jl_save_cont_solution;
jl_function_t* jl_save_array;

jl_function_t* jl_serialize_object;
jl_function_t* jl_deserialize_object;


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

    jl_save_solution = jl_get_function(jl_main_module, "save_solution");
    jl_save_cont_solution = jl_get_function(jl_main_module, "save_cont_solution");
    jl_save_array = jl_get_function(jl_main_module, "save_array");

//    jl_serialize_object = jl_get_function(jl_main_module, "serialize_object");
 //   jl_deserialize_object = jl_get_function(jl_main_module, "deserialize_object");

}

jl_value_t* JL_Interface::jl_array(uint8_t *_ptr, int _size)
{
    jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_uint8_type, 1);
   return (jl_value_t*) jl_ptr_to_array_1d(array_type, _ptr, _size, 0);
 } 

  jl_value_t* JL_Interface::jl_float_array(double* x_vec)
   {
       jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);
       return (jl_value_t*) jl_ptr_to_array_1d(array_type, x_vec, getDim(), 0);
   }


    // Send Julia object via MPI
    void JL_Interface::send_MPI_data(jl_value_t* _dt_ptr, int tag, bool block)
    {
        //LOG_INFO("Sending MPI data: Rank: " << rank << " Tag: " << tag);

        // Serialize the Julia object
        jl_value_t* serialized_data = jl_call1(jl_serialize_obj, _dt_ptr);
        if (jl_exception_occurred())
        {
            std::cerr << "Julia exception occurred during serialization!" << std::endl;
            exit(EXIT_FAILURE);
        }

        // Get serialized data as a byte array
        jl_array_t* data_array = reinterpret_cast<jl_array_t*>(serialized_data);
        size_t data_size = jl_array_len(data_array);
        void* data_ptr = jl_array_data(data_array);

        // Allocate and copy serialized data into the send buffer
        send_buffer = alloc_buffer(data_size);
        memcpy(send_buffer, data_ptr, data_size);

//        LOG_DEBUG("Serialized data copied to send buffer: Rank: " << rank << " Tag: " << tag);

        // Send the data using MPI
        for (int i = 1; i < nproc; i++)
        {
            if (block)
            {
                MPI_Send(send_buffer, data_size, MPI_BYTE, i, tag, MPI_COMM_WORLD);
            }
            else
            {
                MPI_Request request;
                MPI_Isend(send_buffer, data_size, MPI_BYTE, i, tag, MPI_COMM_WORLD, &request);
               // MPI_Wait(&request, MPI_STATUS_IGNORE);
                int flag;
                MPI_Test(&request, &flag, MPI_STATUS_IGNORE);
            }
        }

    }

    void JL_Interface::init_MPI()
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    }

    // Receive Julia object via MPI
    jl_value_t* JL_Interface::receive_MPI_data(int tag, bool block)
    {
//        LOG_INFO("Receiving MPI data: Rank: " << rank << " Tag: " << tag);

        MPI_Status status;
        MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);

        int data_size;
        MPI_Get_count(&status, MPI_BYTE, &data_size);

        uint8_t* recv_buffer = alloc_buffer(data_size);
        MPI_Recv(recv_buffer, data_size, MPI_BYTE, status.MPI_SOURCE, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        jl_value_t* received_data = deserialize_data(recv_buffer, data_size);

        free(recv_buffer);
        return received_data;
    }




