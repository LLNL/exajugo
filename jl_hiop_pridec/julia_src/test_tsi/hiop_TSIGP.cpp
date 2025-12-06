
#include <iostream>
#include <julia.h>

JULIA_DEFINE_FAST_TLS // Required for thread-local storage in Julia

int main() {

    // Step 1: Initialize Julia
    jl_init();

    const char* julia_file_path = "hiop.jl";

    std::string command = "include(\"" + std::string(julia_file_path) + "\")";
    jl_eval_string(command.c_str());

    const char* tsi_julia_file_path = "test_TSACOPF.jl";

    std::string tsi_command = "include(\"" + std::string(tsi_julia_file_path) + "\")";
    jl_eval_string(tsi_command.c_str());

    jl_function_t* test_TSI = jl_get_function(jl_main_module, "test_TSI");
    jl_function_t* load_tsmodel = jl_get_function(jl_main_module, "load_tsmodel");

    std::string  path_to_data = "./data";

// GP model
    std::string model_file = (path_to_data+"/model_state_DSPP_500_0.990_0.884_28.9_9.0_4.2_20000_318.0_300_6_0.021_0.00075.pth");
    std::string data_record = (path_to_data+"/data_record.mat");
 
    jl_value_t* jl_model_file = jl_cstr_to_string(model_file.c_str());
    jl_value_t* jl_data_record = jl_cstr_to_string(data_record.c_str());

    jl_value_t *tsi_model = jl_call2(load_tsmodel, jl_model_file, jl_data_record);

    std::string pf_limit_file = (path_to_data+"/pf_new.mat");
    std::string case_path = (path_to_data+"/ACTIVSg500");

    jl_value_t* jl_pf_limit_file = jl_cstr_to_string(pf_limit_file.c_str());
    jl_value_t* jl_case_path = jl_cstr_to_string(case_path.c_str());

// TSI test
    jl_value_t *tsi = jl_call3(test_TSI, tsi_model, jl_case_path, jl_pf_limit_file);
    
    jl_atexit_hook(0);
    return 0;
}

