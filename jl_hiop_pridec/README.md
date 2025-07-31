
## Contents

Directories:
 - src: C++ source files
 - julia_src: julia source files
 - test_tsi: test tsi constraints

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

  The directory "sbatch_templates" contains the templates to run the code on different machines: 
 
      - ruby.sbatch: contains configuration to run on ruby; 
      - dane.sbatch: contains configuration to run on dane; 
      - hiop.sh: this will be added to any template if it contains the string $HIOP_SH 
          (i.e., the contains of the file hiop.sh will be inserted in the place of the string - except for the first line)

  IMPORTANT: if you wish to override the templates above, write you template on the file default.sbatch and put it in the directory "sbatch_templates" 


## Environment configuration

1) Set the following environment variables:
 
    - PATH_TO_EXAJUGO: path to EXAJUGO;
    - PATH_TO_INSTANCES: path to SACOPF instances;
    - HIOP_INSTALL_DIR: path to HIOP installation;

    The environment variable NORMALIZE_X is used to normalize contingency gradients (default is true):

    export NORMALIZE_X=false  # do not normalize gradients

    if you are using TSI constraints, you have to install tsslope. Clone the following repo:

       https://github.com/SLOPE-grid/tsSLOPE.git

    Set the environment variable PATH_TO_TSSLOPE to the directory containing the cloned repo.

2) Execute on bash:

    source load.sh

## Compilation

  run: ./compile.sh

  creates "build" directory: stores object and executable files


## Execution

   Before executing, the code must be compiled (see **Compilation**)

   1) Configure the environment (see **Environment configuration**)

   2) Generate the batch:

      ./generate_batch.sh  instance time [contingency]

      - instance: bus system
      - time execution time 
      - contingency: contingency file name without the extension, if contingency file name is NOT 'case.con'

      The generated script will be written to "output/scripts/"

   3) After the execution, the generated script will be shown in the screen. If you wish to submit it, press ENTER. Otherwise, press any key. In case you pressed any key OTHER THAN ENTER, you will see this:

    --- Generated batch file: output/scripts/sub_<case>.sbatch ---

    --- RUN: sbatch output/scripts/sub_<case>.sbatch ---

   4) If you wish to submit Run the command displaied in last line of the output:

      sbatch output/scripts/sub_<case>.sbatch

## Checking the current run

   run: ./check.sh

    - this will show you the number of iterations (problems solved) completed by each rank of the most recent case run (or currently running).

    The following columns will be shown:

    rank  # of iterations  date created         last update 

   IMPORTANT: if you the case name as argument, it will filter based on the case name (in addition to the date).

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
