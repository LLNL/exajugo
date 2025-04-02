#!/bin/bash

source ./load.sh

g++ -o run -fPIC -I$JULIA_DIR/include/julia -L$JULIA_DIR/lib -Wl,-rpath,$JULIA_DIR/lib $1 -ljulia -static-libstdc++ --std=c++11

