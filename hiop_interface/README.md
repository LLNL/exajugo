
Contents:

 - load.sh: sets environment variables and loads julia module;
 
 - compiler.sh arg1: compiles a C++ file that uses a julia module:
    arg1: C++ input file
    IMPORTANT: the executable "run" is generated. 
    
 - hiop_9bus.cpp: C++ code that tests the C++-Julia interface for the 9-bus case;
 
 - hiop_500bus.cpp: C++ code that tests the C++-Julia interface for the 500-bus case;

 - hiop_gp.cpp: C++ code that tests the C++-Julia interface for the GP model;

 - hiop.jl: julia interface that contains the julia functions that will be called from C++.

    Functions:

       - load_ACOPF_dir: loads the case given a directory and returns a reference to SCACOPFdata; 
       - load_ACOPF: loads the case given a 3 files (in sequence: raw, rop and con) and returns a reference to SCACOPFdata;
       - solve_base_case: solves base case given a reference (returned by load_ACOPF or load_ACOPF_dir);
       - number_of_contingencies: returns the number of contingencies;
       - load_tsmodel: loads the GP model and returns a reference to a type TSI_GPmodel (holds the GP model as well as the data);
       - test_model: tests dereferencing the model (gets a reference to a type TSI_GPmodel);

 - test_cases.jl: test cases 

       - test_9: returns data for 9-bus system for testing;
       - test_500: returns data for 500-bus system for testing.

 - exajugo_call.jl: intgerface funtions for TSACOPF to test the TSI interface

 - test_TSACOPF.jl: TSI test function that uses the GP model interface

 - test_TSI500.jl: TSI test interface for the 500-bus case.

 
Before running, set the following environment variable:

  - PATH_TO_EXAJUGO: path to exajugo installation; 
  - PATH_TO_TSSLOPE: path to tsSLOPE installation.

To test the 9bus case from julia (replace 9bus with 500bus for testing the 500-bus case):

  2) julia

  In the julia REPL:

    julia> include("test_cases.jl") 
    julia> ptr = test_bus9()  
    julia> solve_base_case(ptr)
 

Steps to run the 9bus case from C++ (replace 9bus with 500bus for testing the 500-bus case):

  1) set the environment variable PATH_TO_EXAJUGO
  2) compile:
      ./compile.sh hiop_9bus.cpp

  3) run:
       ./run  
