#!/bin/bash

export JULIA_VERSION=julia-1.11.5
export MODULEPATH=/usr/workspace/llopt/modulefiles/:$MODULEPATH
export LD_LIBRARY_PATH=/usr/workspace/llopt/julia/$JULIA_VERSION/lib:/usr/workspace/llopt/julia/$JULIA_VERSION/lib/julia/:$LD_LIBRARY_PATH

module load julia

export JULIA_DIR=/usr/workspace/llopt/julia/$JULIA_VERSION/



