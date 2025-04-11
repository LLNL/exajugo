
#include <iostream>
#include <julia.h>

JULIA_DEFINE_FAST_TLS // Required for thread-local storage in Julia

int main() {

    // Step 1: Initialize Julia
    jl_init();

    const char* julia_file_path = "hiop.jl";

    std::string command = "include(\"" + std::string(julia_file_path) + "\")";
    jl_eval_string(command.c_str());

    jl_function_t* load_ACOPF = jl_get_function(jl_main_module, "load_ACOPF");
    jl_function_t* load_ACOPF_dir = jl_get_function(jl_main_module, "load_ACOPF_dir");

    jl_function_t* number_of_contingencies = jl_get_function(jl_main_module, "number_of_contingencies");
    jl_function_t* solve_base_case = jl_get_function(jl_main_module, "solve_base_case");

    std::string exajugo_path = std::getenv("PATH_TO_EXAJUGO");
    std::string inst_dir = (exajugo_path+"/examples/500bus/");
    
    jl_value_t* jl_inst_dir = jl_cstr_to_string(inst_dir.c_str());

    jl_value_t *opt_data = jl_call1(load_ACOPF_dir, jl_inst_dir);

    jl_value_t *ncont = jl_call1(number_of_contingencies, opt_data);

    int64_t jl_ncont = jl_unbox_int64(ncont);

    std::cout<< " # of contingencies: "<<jl_ncont<<std::endl;

    jl_value_t *val = jl_call1(solve_base_case, opt_data);
    
    jl_atexit_hook(0);
    return 0;
}

