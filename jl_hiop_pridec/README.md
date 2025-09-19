
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
 - run_interactive.jl: generates script to run interactively in a debug node for slurm;
 - test_cases.jl: functions that load test instances;
 - test.jl: function that loads bus9 case and solves the base case.

Scripts:

 - load.sh: load the necessary environment variables to compile/run the code;
 - compile.sh: compiles the run;
 - generate_batch.sh: calls julia with generate_sbatch.jl.
 - run_interactive.sh: calls julia with run_interactive.jl.

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

    The environment variable NORMALIZE_X is used to normalize contingency gradients (default is false):

    export NORMALIZE_X=true  # do not normalize gradients

    if you are using TSI constraints, you have to install tsslope. Clone the following repo:

       https://github.com/SLOPE-grid/tsSLOPE.git

    Set the environment variable PATH_TO_TSSLOPE to the directory containing the cloned repo.

2) Execute on bash:

    source load.sh

## Compilation

  run: ./compile.sh

  creates "build" directory: stores object and executable files


## Execution (job submission)

   Before executing, the code must be compiled (see **Compilation**)

   1) Configure the environment (see **Environment configuration**)

   2) Generate the batch:

      ./generate_batch.sh  instance time [contingency]

      - instance: bus system
      - time execution time 
      - contingency: contingency file name without the extension, if contingency file name is NOT 'case.con'

      The generated script will be written to "output/scripts/"

   3) During execution, you will be able to enter the number of tasks (>=2), as then following will be printed:

         "Enter # of tasks or press ENTER to set # of tasks= # of contingencies+1):"

         If you press ENTER, then the script will automatically set the # of processors = # of contingencies+1

      IMPORTANT: if you set the environment variable NTASKS, its value will be used as the # of processors.


   4) After the execution, the generated script will be shown in the screen. If you wish to submit it, press ENTER. Otherwise, press any key. In case you pressed any key OTHER THAN ENTER, you will see this:

    --- Generated batch file: output/scripts/sub_<case>.sbatch ---

    --- RUN: sbatch output/scripts/sub_<case>.sbatch ---

   5) If you wish to run the command displaied in last line of the output:

      sbatch output/scripts/sub_<case>.sbatch


## Execution (interactive pdebug node allocation)

   Before executing, the code must be compiled (see **Compilation**)

   1) Configure the environment (see **Environment configuration**)

   2) Generate the batch:

      ./run_interactive.sh instance [contingency]

      - instance: bus system
      - contingency: contingency file name without the extension, if contingency file name is NOT 'case.con'

      The generated script will be written to "output/scripts/"

   3) After the execution, the generated script will be shown in the screen. If you wish to run it, press ENTER. Otherwise, press any key. In case you pressed any key OTHER THAN ENTER, you will see this:

    --- Running command: ./output/scripts/run_<case>.sh  ---
    Output directory: output/$OUTPUT_DIR

   4) If you wish to run the command displaied in last line of the output:

      ./output/scripts/run_<case>.sh

   IMPORTANT: before running the code interactively, you must allocate a node with the necessary number of processors using the following command:

         salloc -N1 -nNTASKS -ppdebug

   where NTASKS = number of contingencies + 1


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
   
   Outputs are written to two different directories:

      - rank_logs: logs of all ranks are written to this directory, written to different directories according to rank;

      - problems: solution and iteration data for the master problem as well as for each contingency problem are saved to this directory in different directories, according to the contingency number (0 for master).

   Logs: the logs saved to rank_logs directory are the following:
      - log_case_JOBNUMBER_0.out: screen output;
      - log_case_JOBNUMBER_0.err: slurm errors;

   Iteration data:
 
     Master problem:
 
         - solution.csv: solution for each iteration;
         - iterations.csv: objective value and execution time for each iteration;
         - hessian.csv: hessian value for each iteration;
         - gradient.csv: gradient value for each iteration.

     Contingency problems:
   
         - solution_i.csv: contingency subproblem solution for each iteration;
         - iterations_i.csv: contingency subproblem objective value and execution time for each iteration;
   
## Solver options

    You can change the ipopt solver options for the base case and contingency cases separately by setting the following environment variables:

      - BASE_CASE_OPTIONS: name of the file that contains the ipopt options for the base case problem;
      - CONTINGENCY_CASE_OPTIONS: name of the file that contains the ipopt options for the contingency problems.

   The files format for both files above are a text file with two columns: option value (without header)
      Example:
           sb yes
           tol 1e-6


       





