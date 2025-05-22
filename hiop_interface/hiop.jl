

#--
#using Pkg; Pkg.add(["MPI", "Serialization", "DataFrames"])


push!(LOAD_PATH, string(ENV["PATH_TO_EXAJUGO"], "/modules"))
push!(LOAD_PATH, string(ENV["PATH_TO_EXAJUGO"], "/modules/SCACOPFSubproblems.jl"))


using CSV
using DataFrames
using Ipopt, JuMP
using SCACOPFSubproblems
using Serialization

include(string(ENV["PATH_TO_EXAJUGO"], "/modules/SCACOPFSubproblems/starting_point.jl"))
include(string(ENV["PATH_TO_EXAJUGO"], "/modules/corejugo/solution_evaluator.jl"))


#--- python tslope desabled
#include(string(ENV["PATH_TO_TSSLOPE"], "/init.jl"))


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

#-- used for debug
function debug_base_case(ptr)

   println(" Julia structure: ")
   println(ptr[])
   println(" -----")
end

function build_model(ptr)

   return Ref(basecase_model_pridec(ptr[]))

end

function get_optimimizer()
 
    return optimizer_with_attributes(Ipopt.Optimizer, "sb" => "no", "print_level" => 0)

end


function getModel(ptr)

   return Ref(basecase_model_pridec(ptr[], get_optimimizer()))

end

struct PriDecBasecaseSolution

   solution::BasecaseSolution
   model::Model
   PriDecBasecaseSolution(_sol::BasecaseSolution, _m::Model) = new(_sol, _m)
end

struct PriDecContingencySolution

   solution::BasecaseSolution
   model::Model
   PriDecBasecaseSolution(_sol::BasecaseSolution, _m::Model) = new(_sol, _m)
end


function getCost(ptr)
   return ptr[].cont_cost
end

function getGradientDim(ptr)
   return length(ptr[].cont_grad)
end

function getGradient(ptr, x)
   for (i,v) in enumerate(ptr[].cont_grad)
       x[i] = v
   end
end


function getSolution(ptr, x)
   #println("getSolution(ptr, x)getSolution(ptr, x)")
   #println(ptr[])
   #println("ptr[].p_g")
   #println(ptr[].p_g)
   #println(" -- END--9")
   for (i,v) in enumerate(ptr[].p_g)
       x[i] = v
   end
end


function getObjective(ptr)
   return ptr[].base_cost
end

function getDim(ptr)
   return nrow(ptr[].G)
end


function copy_ACOPF(ptr)
   return Ref(deepcopy(ptr[]))
end

struct ObjectBytes

   size::Int64
   bytes::Vector{UInt8}
   ObjectBytes(_s::Int64, _d::Vector{UInt8}) = new(_s, _d)
end

function get_data_size(ptr)
    return ptr[].size

end

function get_data_bytes(ptr, xdata)
    for (i,v) in enumerate(ptr[].bytes)
       xdata[i] = v
   end
end

function serialize_obj(ptr)

    #println("ptrint")
    io = IOBuffer()
    serialize(io, ptr[])
    bytes = take!(io)

    return Ref(ObjectBytes(length(bytes),bytes))
end

function deserialize_obj(bytes) 

    io = IOBuffer(bytes)
    received_data = deserialize(io)

    return Ref(received_data)
end


function process_data(data; callback::T=nothing,  tmpcallback::T=nothing) where {T <: Union{Nothing, Function}}
    if callback == nothing
        return data
    else
        return callback(data)
    end
end

struct RecourseDerivatives2
    gradient::Ptr{Cdouble}
    hessian::Ptr{Cdouble}
    RecourseDerivatives(_grad, _hess) = new(_grad, _hess)
end


struct RecourseDerivatives
    gradient::Ref{Vector{Float64}}
    hessian::Ref{Vector{Float64}}
    RecourseDerivatives(_grad, _hess) = new(_grad, _hess)
end

function get_recourse_derivatives(grad, hess)
    return Ref(RecourseDerivatives(grad, hess))
end

using LinearAlgebra



function save_array(file_path, ptr, grad)

    save_opt_data(file_path, ("norm"=>norm(grad)), [round(x, digits=5) for x in grad])

end

function save_cont_solution(file_path, ptr, prev_sol)

    save_opt_data(file_path, ("objective"=>prev_sol[].cont_cost), [round(x, digits=5) for x in prev_sol[].p_g])

end

function save_solution(file_path, ptr, prev_sol)

    save_opt_data(file_path, ("objective"=>prev_sol[].base_cost), [round(x, digits=5) for x in prev_sol[].p_g])

end

function save_opt_data(file_path, kval::Pair{String, Float64}, sol)

  fd_name = kval.first
  value = kval.second
  if isfile(file_path)
        # If the file exists, read it
        existing_data = CSV.read(file_path, DataFrame)
        
        # Determine the number of rows in the existing file
        num_rows = size(existing_data, 1)
        
        columns = Dict(
        :iteration => num_rows,
        Symbol(fd_name) => value,
        [Symbol("x_$i") => sol[i] for i in 1:length(sol)]...)

        # Create a new row with iteration set to num_rows and the given objective
        new_data = DataFrame(columns)
        
        # Append the new row to the existing data
        updated_data = vcat(existing_data, new_data)
        
        # Write the updated data back to the file
        CSV.write(file_path, updated_data)
    else
        # If the file does not exist, create it with iteration set to 0 and the given objective
        columns = Dict(
        :iteration => 0,
        Symbol(fd_name) => value,
        [Symbol("x_$i") => sol[i] for i in 1:length(sol)]...)

        new_data = DataFrame(columns)
        CSV.write(file_path, new_data)
    end

end

function solve_base_case_recourse(ptr, prev_sol, ptr_rderivaties)

   gnorm = norm(ptr_rderivaties[].gradient[])
   G = ptr_rderivaties[].gradient[]/gnorm
   hnorm = norm(ptr_rderivaties[].hessian[])
   H = ptr_rderivaties[].hessian[]/hnorm

   recourse_fx = (args...) ->  begin x = collect(args); (1/2)*(x.^2)'H + (G - H.*x)'x end
   recourse_gx = (argG, args...) ->   begin  x = collect(args);  argG .= G.*x;  end
   recourse_Hx = (argH, args...) ->   begin  x = collect(args); argH[diagind(argH)].=H;  end

    return Ref(solve_basecase(ptr[], get_optimimizer(), 
              recourse_f=recourse_fx, recourse_g=recourse_gx, recourse_H=recourse_Hx,
              previous_solution=prev_sol[])[1])

end

function solve_base_case(ptr)

   return Ref(solve_basecase(ptr[], get_optimimizer())[1])

end


function number_of_contingencies(ptr)

   return length(ptr[].cont_labels)

end

function number_of_columns(ptr)

   return nrow(ptr[].G)

end


function solve_contingency_pridec(ptr, i::Int64, ptr_basesol)
 
   ptr_basesol[].psd_hash = hash(ptr[])
   return Ref(solve_contingency(ptr[], i, ptr_basesol[], get_optimimizer(),quadratic_relaxation_k=0.5))


end


