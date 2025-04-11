
using Pkg
Pkg.activate((ENV["PATH_TO_EXAJUGO"]))
push!(LOAD_PATH, string(ENV["PATH_TO_EXAJUGO"], "/modules"))

push!(LOAD_PATH, string(ENV["PATH_TO_EXAJUGO"], "/modules/SCACOPFSubproblems.jl"))

using DataFrames
using Ipopt, JuMP
using SCACOPFSubproblems

include(string(ENV["PATH_TO_EXAJUGO"], "/modules/SCACOPFSubproblems/starting_point.jl"))
include(string(ENV["PATH_TO_EXAJUGO"], "/modules/corejugo/solution_evaluator.jl"))

include(string(ENV["PATH_TO_TSSLOPE"], "/init.jl"))

# GM model structure

struct TSI_GPmodel
  GPmodel::Dict{Any,Any}
  data::Matrix{Float32}
  TSI::Vector{Float32}
end


function load_tsmodel(model_file, data_record)

  GPmodel, data, TSI = tsslope_lib.load_model(model_file, data_record);

  return Ref(TSI_GPmodel(GPmodel, data, TSI))

end

function test_model(refmodel)

   TSIGPmodel = refmodel[]
   println(" Object derefecened successfully! ")
   return 0;

end

function load_ACOPF_dir(prob)

   opt_data=SCACOPFdata(prob)
   return Ref(opt_data)
 
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


