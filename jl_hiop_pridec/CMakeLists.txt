

add_compile_options("-I" $ENV{JULIA_DIR}include/julia "-ljulia" "-static-libstdc++")
add_link_options("-L" $ENV{JULIA_DIR}lib "-Wl,-rpath,$ENV{JULIA_DIR}lib" "-ljulia" "-static-libstdc++")

#add_link_options("-I" $ENV{JULIA_DIR}include/julia "-L" $ENV{JULIA_DIR}lib "-ljulia" "-static-libstdc++")
#target_link_libraries(NlpPriDecEx1.exe "-Wl,--disable-new-dtags")

#add_compile_options("-fPIC" "-I" $ENV{JULIA_DIR}include/julia "-ljulia" "-static-libstdc++")
#add_link_options("-fPIC" "-L" $ENV{JULIA_DIR}lib "-Wl,-rpath,$ENV{JULIA_DIR}lib" "-ljulia" "-static-libstdc++")
#  add_executable(jl_NlpPriDec.exe jl_NlpPriDec.cpp jl_interface.cpp jl_NlpPriDecDriver.cpp)


if(HIOP_USE_MPI)
  add_executable(jl_NlpPriDec.exe jl_NlpPriDec.cpp jl_interface.cpp jl_NlpPriDecDriver.cpp)
  target_link_libraries(jl_NlpPriDec.exe HiOp::HiOp)
endif()


##########################################################
# CMake Tests
##########################################################
if(HIOP_USE_MPI)
  add_test(NAME NlpPriDec1_1 COMMAND ${RUNCMD} "$<TARGET_FILE:jl_NlpPriDec.exe>" "-selfcheck")
  add_test(NAME NlpPriDec1_mpi COMMAND ${MPICMD} -n 2 "$<TARGET_FILE:jl_NlpPriDec.exe>" "-selfcheck")
endif(HIOP_USE_MPI)


