cmake_minimum_required(VERSION 3.18)
project(jl_NlpPriDec C CXX Fortran)

find_package(MPI REQUIRED COMPONENTS C CXX Fortran)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

include_directories($ENV{HIOP_INCLUDE_DIR})
link_directories($ENV{HIOP_LIBRARY_DIR})

add_executable(jl_NlpPriDec.exe jl_NlpPriDec.cpp jl_interface.cpp jl_NlpPriDecDriver.cpp)

target_include_directories(jl_NlpPriDec.exe PRIVATE
    ${MPI_C_INCLUDE_DIRS}
    $ENV{JULIA_DIR}include/julia
)

target_link_libraries(jl_NlpPriDec.exe
    hiop
    ${MPI_C_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    -L$ENV{JULIA_DIR}lib
    -Wl,-rpath,$ENV{JULIA_DIR}lib
    -ljulia
    -static-libstdc++
    -std=gnu++17
    -lstdc++fs
)