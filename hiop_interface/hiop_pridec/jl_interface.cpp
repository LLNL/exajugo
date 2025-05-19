

#include "jl_interface.hpp"

//#include <julia.h>


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
jl_function_t* jl_save_array;



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
    jl_save_array = jl_get_function(jl_main_module, "save_array");

}

jl_value_t* JL_Interface::jl_array(uint8_t *_ptr, int _size)
{
    jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_uint8_type, 1);
   return (jl_value_t*) jl_ptr_to_array_1d(array_type, _ptr, _size, 0);
 } 


void JL_Interface::release_buffer()
{
      if (send_buffer != nullptr)
      {
         delete []send_buffer; 
         send_buffer=nullptr;
        }
 
}

uint8_t *JL_Interface::alloc_buffer(int size)
{
    if (size != size_buffer) 
    	release_buffer();

    if (send_buffer == nullptr)
    {
        send_buffer=new uint8_t[size];
        size_buffer=size;
     }

   return send_buffer;

}

   JL_Interface::JL_Interface(): cont_sol(nullptr), base_sol(nullptr), send_buffer(nullptr), size_buffer(0)
    { 

      init_MPI();

      if (rank==source_process)
      {
         opt_data = read_data(); 
         send_instance();
      }
      else
        receive_instance();

   }


void JL_Interface::send_MPI_data(jl_value_t* _dt_ptr, int tag, bool block)
{
 
   int size;

   send_buffer= (uint8_t *)serialize_data(size, _dt_ptr);  //jl_call2(jl_serialize_obj, pobj, refsz);
  
  if (tag>0)
    for (int i=1;i<nproc;i++)
    {

       if (block)
         MPI_Send(send_buffer, size, MPI_CHAR, i, tag,MPI_COMM_WORLD);
       else
       {
         MPI_Request request;

         MPI_Isend(send_buffer, size, MPI_CHAR, i, tag,MPI_COMM_WORLD,&request);       
       }

    }

 else
    MPI_Bcast(send_buffer, size, MPI_CHAR, source_process, MPI_COMM_WORLD);


  }

void JL_Interface::init_MPI() 
{    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	source_process=0; 
}


jl_value_t* JL_Interface::receive_MPI_data(int tag, bool block)
{

      MPI_Status status;
      MPI_Probe(source_process, tag, MPI_COMM_WORLD, &status);

      int count;
      MPI_Get_count(&status, MPI_CHAR, &count);

      uint8_t *retptr = alloc_buffer(count); //(uint8_t *)malloc(count*sizeof(uint8_t));

      MPI_Recv(retptr, count, MPI_CHAR, source_process, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      jl_value_t* pdtdata= deserialize_data(retptr, count);
   
      return pdtdata;
 }









