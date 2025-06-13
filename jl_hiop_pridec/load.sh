#!/bin/bash

export MODULEPATH=/usr/workspace/llopt/modulefiles/:$MODULEPATH
export LD_LIBRARY_PATH=/usr/workspace/llopt/julia/julia-1.10.2/lib:/usr/workspace/llopt/julia/julia-1.10.2/lib/julia/:$LD_LIBRARY_PATH

module load julia

export JULIA_DIR=/usr/workspace/llopt/julia/julia-1.10.2/



