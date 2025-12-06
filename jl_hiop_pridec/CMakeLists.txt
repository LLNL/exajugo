cmake_minimum_required(VERSION 3.18)
project(jl_NlpPriDec C CXX Fortran)

set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set ABI globally
add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=1)

find_package(MPI REQUIRED COMPONENTS C CXX Fortran)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

include_directories($ENV{HIOP_INCLUDE_DIR})
link_directories($ENV{HIOP_LIBRARY_DIR})

add_executable(jl_NlpPriDec.exe
    src/jl_NlpPriDec.cpp
    src/jl_interface.cpp
    src/jl_NlpPriDecDriver.cpp
)

target_include_directories(jl_NlpPriDec.exe PRIVATE
    ${MPI_C_INCLUDE_DIRS}
    $ENV{JULIA_DIR}include/julia
)

target_link_directories(jl_NlpPriDec.exe PRIVATE
    $ENV{JULIA_DIR}lib
)

target_link_libraries(jl_NlpPriDec.exe
    hiop
    ${MPI_CXX_LIBRARIES} # Use CXX libraries instead of C
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    julia
    m
    dl
    pthread
)

set_target_properties(jl_NlpPriDec.exe PROPERTIES
    BUILD_RPATH "$ENV{JULIA_DIR}lib"
    INSTALL_RPATH "$ENV{JULIA_DIR}lib"
    LINKER_LANGUAGE CXX
)
