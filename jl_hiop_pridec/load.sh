#!/bin/bash

export PATH_TO_HSLLIB="/usr/WS2/hiop/software/spack/opt/spack/linux-rhel8-icelake/gcc-10.3.1/coinhsl-2015.06.23-nszs3vct4qvtgythawqu2kdue5zx5pjj/lib/"

export JULIA_PROJECT=$PATH_TO_EXAJUGO 

export HIOP_INCLUDE_DIR="$HIOP_INSTALL_DIR/include"
export HIOP_LIBRARY_DIR="$HIOP_INSTALL_DIR/lib"

export JULIA_VERSION_NUMBER=1.11.5
export JULIA_VERSION=julia-$JULIA_VERSION_NUMBER
export MODULEPATH=/usr/workspace/llopt/modulefiles/:$MODULEPATH
export LD_LIBRARY_PATH=/usr/workspace/llopt/julia/$JULIA_VERSION/lib:/usr/workspace/llopt/julia/$JULIA_VERSION/lib/julia/:$LD_LIBRARY_PATH

export JULIA_DIR=/usr/workspace/llopt/julia/$JULIA_VERSION/

module load mkl
module load julia/v$JULIA_VERSION_NUMBER

export LD_LIBRARY_PATH="$PATH_TO_HSLLIB:${LD_LIBRARY_PATH}"



