Before compiling:

source config.sh

Entry point: jl_NlpPriDecDriver.cpp 


Class JL_Interface was created to interface Julia objects with pridec


Execution:

srun -n 4  ./jl_NlpPriDec.exe

OR:

srun -n 4 --verbose ./jl_NlpPriDec.exe
