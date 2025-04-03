
Contents:

 - load.sh: sets environment variables and loads julia module;
 
 - compiler.sh arg1: compiles a C++ file that uses a julia module:
    arg1: C++ input file
    IMPORTANT: the executable "run" is generated. 
    
 - hiop.cpp: C++ code that tests the C++-Julia interface;

 - hiop.jl: julia interface that contains the julia functions that will be called from C++.
    IMPORTANT: the main function that loads the problems (load_case) returns a pointer.

Before running, set the following environment variable:
  - PATH_TO_EXAJUGO: path to exajugo installation 

To test from julia:

  2) julia

  In the julia REPL:

    julia> include("hiop.jl") 
    julia> ptr = test_bus9()  
    julia> solve_base_case(ptr)
 

Steps to run:

  1) set the environment variable PATH_TO_EXAJUGO
  2) compile:
      ./compile.sh hiop.cpp

  3) run:
       ./run  
