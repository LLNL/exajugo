

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

OPF_DATA=nothing;

function get_data_ptr()

   return OPF_DATA

end

function load_ACOPF_dir(prob)

   opt_data=SCACOPFdata(prob)
   global OPF_DATA
   OPF_DATA =Ref(opt_data)

   return OPF_DATA
end

function load_ACOPF(rawfile, ropfile, confile)

   opt_data=SCACOPFdata(raw_filename=rawfile, rop_filename=ropfile, con_filename=confile)
   global OPF_DATA
   OPF_DATA =Ref(opt_data)

   return OPF_DATA

#   return Ref(opt_data)
 
end

#-- used for debug
function debug_base_case(ptr)

   println(" Julia structure: ")
   println(ptr[])
   println(" -----")
end

function debug_array(x)

   println(" --- debug_array ---")
   println(x)
   println("----")
   for (i,v) in enumerate(x)
       println(x[i])
   end
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

   #gnorm = norm(ptr_rderivaties[].gradient[])
   #G = ptr_rderivaties[].gradient[]/gnorm


function getGradient(ptr, x)
   for (i,v) in enumerate(ptr[].cont_grad)
       x[i] = v
   end
   x .= normalize(x)
end

#function array_to_struct(prob::Ref{SCACOPFdata}, arr, array_lengths::Ref{Dict{Symbol, Int}})

MASTER_SOL = nothing

function array_to_struct(prob::Ref{SCACOPFdata}, arr, array_lengths::Ref{Dict{Symbol, Int}})

    global MASTER_SOL
    MASTER_SOL= array_to_struct_generic(prob, arr, BasecaseSolution, array_lengths)
    return MASTER_SOL

end

#::Vector{Float64}
#function array_to_struct_generic(prob::Ref{SCACOPFdata}, arr, T::Type, array_lengths::Ref{Dict{Symbol, Int}})

function array_to_struct_generic(prob, arr, T::Type, array_lengths)
    field_types = fieldtypes(T)
    field_names = fieldnames(T)

    
    # Prepare to extract fields from the array
    fields = []
    idx = 1  # Start index for the array
    
    # Skip the first field
    for (name, field_type) in zip(field_names[2:end], field_types[2:end])
        if field_type <: AbstractVector  # If the field is an array

            len = array_lengths[][name]  # Get the length of this array field
            push!(fields, arr[idx:idx+len-1])  # Extract the array slice
            idx += len

        else  # If the field is a scalar
            push!(fields, arr[idx])  # Extract the scalar value
            idx += 1

        end
    end
    
    # Ensure all array elements are consumed
    @assert idx - 1 == length(arr) "Array size does not match structure requirements"
    
    # Construct the structure
    return Ref(T(prob[], fields...))
end

function struct_to_array_generic!(s::Ref{T}, arr::Vector{Float64}) where T
    idx = 1  # Start writing at the first index of the array
    
    # Skip the first field
    for name in fieldnames(T)[2:end]
        value = getfield(s[], name)  # Dereference the Ref to access the value
        if value isa AbstractVector  # If the field is an array
            for v in value
                arr[idx] = v  # Write each element of the array into arr
                idx += 1  # Move to the next index
            end
        else  # If the field is a scalar
            arr[idx] = value  # Write the scalar value into arr
            idx += 1  # Move to the next index
        end
    end

    # Ensure we don't exceed the allocated size
    @assert idx - 1 <= length(arr) "Array size exceeded during struct conversion!"
end

#:psd_hash, :v_n, :theta_n, :b_s, :p_g, :q_g, :base_cost, :recourse_cost)

FIELD_SIZES_DICT=nothing

function define_array_lengths(prob::Ref{SCACOPFdata})
    global FIELD_SIZES_DICT
    FIELD_SIZES_DICT= Ref(Dict(
        :v_n => size(prob[].N, 1),  # Same size as the number of rows in `prob.N`
        :theta_n => size(prob[].N, 1),  # Same size as the number of rows in `prob.N`
        :b_s => size(prob[].SSh, 1),  # Same size as the number of rows in `prob.SSh`
        :p_g => size(prob[].G, 1),  # Same size as the number of rows in `prob.N`
        :q_g => size(prob[].G, 1),  # Same size as the number of rows in `prob.N`
        :base_cost => 1,  # Assuming scalar size for `base_cost`
        :recourse_cost => 1  # Assuming scalar size for `recourse_cost`
    ))
    return FIELD_SIZES_DICT
end


function full_solution_dim(fieldsizes)

    return sum(values(fieldsizes[]))

end

function getSolution(ptr, x)

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

function serialize_obj(obj)
    buffer = IOBuffer()
    obj_copy = deepcopy(obj)  # Create a deep copy of the object
    serialize(buffer, obj_copy)  # Serialize the copy
    return take!(buffer)  # Return serialized data and the copy
end


function deserialize_obj(data::Vector{UInt8})
    buffer = IOBuffer(data)
    obj = deserialize(buffer)  # Deserialize the object
    return deepcopy(obj)  # Return a deep copy of the deserialized object
end


function serialize_obj_OLD(ptr)

    io = IOBuffer()
    serialize(io, ptr[])
    bytes = take!(io)

    return Ref(ObjectBytes(length(bytes),bytes))
end

function deserialize_obj_OLD(bytes) 

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

#-- 

struct RecourseDerivativesRef
    gradient::Ref{Vector{Float64}}
    hessian::Ref{Vector{Float64}}
    RecourseDerivativesRef(_grad, _hess) = new(Ref(_grad), Ref(_hess))
end

function get_recourse_derivatives_ref(grad, hess)
    return Ref(RecourseDerivativesRef(grad, hess))
end


struct RecourseDerivatives

    gradient::Vector{Float64}
    hessian::Vector{Float64}
   # RecourseDerivativesAlloc(_grad, _hess) = new(_grad, _hess)
    function RecourseDerivatives(_grad, _hess, _len) 

        grad_copy = Vector{Float64}(undef, _len)
        hess_copy = Vector{Float64}(undef, _len)

        grad_copy .= _grad[1:_len]
        hess_copy .= _grad[1:_len]

        new(grad_copy, hess_copy)
  
    end
end

function get_recourse_derivatives(grad, hess, _len)
    return Ref(RecourseDerivatives(grad, hess, _len))
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

#x .= normalize(x)

SOLUTION_WITH_RECOURSE=nothing

function solve_base_case_recourse(ptr, prev_sol, ptr_rderivaties)

#   G = ptr_rderivaties[].gradient[]
#   H = ptr_rderivaties[].hessian[]
   G = ptr_rderivaties[].gradient
   H = ptr_rderivaties[].hessian

   recourse_fx = (args...) ->  begin x = collect(args); (1/2)*(x.^2)'H + (G - H.*x)'x end
   recourse_gx = (argG, args...) ->   begin  x = collect(args);  argG .= G.*x;  end
   recourse_Hx = (argH, args...) ->   begin  x = collect(args); argH[diagind(argH)].=H;  end

   global SOLUTION_WITH_RECOURSE
    SOLUTION_WITH_RECOURSE= Ref(solve_basecase(ptr[], get_optimimizer(), 
              recourse_f=recourse_fx, recourse_g=recourse_gx, recourse_H=recourse_Hx,
              previous_solution=prev_sol[])[1])
   # allocated_bytes = Base.gc_bytes() 
   #println("SOLUTION_WITH_RECOURSE Memory allocated: ", allocated_bytes, " bytes")

    return SOLUTION_WITH_RECOURSE
end

#    SOLUTION_WITH_RECOURSE= Ref(solve_basecase(ptr[], get_optimimizer(), 
 #             recourse_f=recourse_fx, recourse_g=recourse_gx, recourse_H=recourse_Hx,
  #            previous_solution=prev_sol[])[1])

function solve_base_case(ptr)

   global SOLUTION_WITH_RECOURSE
   println("SOLUTION_WITH_RECOURSE\n\n")
   sol=solve_basecase(ptr[], get_optimimizer())[1]
   SOLUTION_WITH_RECOURSE= Ref(sol)
#   SOLUTION_WITH_RECOURSE= Ref(solve_basecase(ptr[], get_optimimizer())[1])

   #allocated_bytes = Base.gc_bytes()
   #println("BASE Memory allocated: ", allocated_bytes, " bytes")
   println("===>SOLUTION_WITH_RECOURSE\n\n")

   return SOLUTION_WITH_RECOURSE

end


function number_of_contingencies(ptr)

   return length(ptr[].cont_labels)

end

function number_of_columns(ptr)

   return nrow(ptr[].G)

end

CONT_SOL=nothing

function solve_contingency_pridec(ptr, i::Int64, ptr_basesol)
 
   ptr_basesol[].psd_hash = hash(ptr[])

   global CONT_SOL
   CONT_SOL = Ref(solve_contingency(ptr[], i, ptr_basesol[], get_optimimizer()))
   return CONT_SOL

end


