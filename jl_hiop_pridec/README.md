
## Contents

C++ pridec interface:

 - jl_NlpPriDecDriver.cpp: entry point for the program;
 - jl_NlpPriDec.hpp: pridec interface class header;
 - jl_NlpPriDec.cpp: pridec interface class implementation.

C++ classes that interface with julia code in gthe file "hiop.jl":

 - jl_interface.hpp: C++-julia interface class header;
 - jl_interface.cpp: C++-julia interface class implementation.

Julia code:

 - hiop.jl: contains the functions that interface with exajugo;
 - generate_sbatch.jl: generates sbatch file for slurm;
 - test_cases.jl: functions that load test instances;
 - test.jl: function that loads bus9 case and solves the base case.

Scripts:

 - load.sh: load the necessary environment variables to compile/run the code;
 - compile.sh: compiles the run;
 - generate_batch.sh: calls julia with generate_sbatch.jl.

Script templates:

 - template.sbatch: this script is used by generate_sbatch.jl. 


## Environment configuration


1) Set the following environment variables:
 
    - PATH_TO_EXAJUGO: path to EXAJUGO;
    - PATH_TO_INSTANCES: path to SACOPF instances;
    - PATH_TO_HSLLIB: path to HSLLIB;
    - HIOP_INSTALL_DIR: path to HIOP installation;

2) Execute on bash:

    source load.sh

## Compilation

  run: ./compile.sh


## Execution

   Before executing, the code must be compiled (see **Compilation**)

   1) Configure the environment (see **Environment configuration**)

   2) Generate the batch:

      ./generate_batch.sh  instance time [contingency]

      - instance: bus system
      - time execution time 
      - contingency: contingency file name without the extension, if contingency file name is NOT 'case.con'

   3) Run the command displaied in last line of the output:

      sbatch sub_<case>.sbatch

## Outputs

   All output is written to a subdirectory of the directory "output".

   The root output directory for a given case: case_$(date +%Y%m%d_%H%M%S)

   All output files are saved in different directories according to the ranks of the processes: rank_0, rank_1, ...

   Outputs for rank 0 for a given case (bnus system):

      - solution.csv: solution for each iteration;
      - solution_iterations.csv: objective value and execution time for each iteration;
      - log_case_JOBNUMBER_0.out: screen output;
      - log_case_JOBNUMBER_0.err: slurm errors;
      - hessian.csv: hessian value for each iteration;
      - gradient.csv: gradient value for each iteration.
   
   Outputs for rank i > 0 for a given case (bnus system):

      - contingency_i.csv: contingency subproblem solution for each iteration;
      - contingency_i_iterations.csv: contingency subproblem objective value and execution time for each iteration;
      - log_case_JOBNUMBER_i.out: screen output;
      - log_case_JOBNUMBER_i.err: slurm errors.
