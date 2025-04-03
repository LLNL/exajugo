
using Pkg
Pkg.activate((ENV["PATH_TO_EXAJUGO"]))
push!(LOAD_PATH, string(ENV["PATH_TO_EXAJUGO"], "/modules"))

push!(LOAD_PATH, string(ENV["PATH_TO_EXAJUGO"], "/modules/SCACOPFSubproblems.jl"))

using DataFrames
using Ipopt, JuMP
using SCACOPFSubproblems

include(string(ENV["PATH_TO_EXAJUGO"], "/modules/SCACOPFSubproblems/starting_point.jl"))
include(string(ENV["PATH_TO_EXAJUGO"], "/modules/corejugo/solution_evaluator.jl"))


function test_9bus()

    prob="/p/lustre1/santiago/COSMIN/GITHUB/exajugo/examples/9bus/";
    rawfile=joinpath(prob,"case.raw"); confile=joinpath(prob,"case.con"); ropfile=joinpath(prob,"case.rop");

    return load_ACOPF(rawfile, ropfile, confile)

end


function load_ACOPF(rawfile, ropfile, confile)

   opt_data=SCACOPFdata(raw_filename=rawfile, rop_filename=ropfile, con_filename=confile)
   return Ref(opt_data)
 
end


function solve_base_case(ptr)

    psd = ptr[]
    opt = optimizer_with_attributes(Ipopt.Optimizer, "sb" => "no")

    sol= solve_basecase(psd, opt)
    return 1

end


function number_of_contingencies(ptr)

   obj = ptr[]
   return length(obj.cont_labels)

end


function solve_contingency(ptr, i, ptr_grad, ptr_x)


end


