
for pkg in ["CSV", "Revise", "MKL"]
    try
        @eval using $(Symbol(pkg))
    catch
        import Pkg
        Pkg.add(pkg)
        @eval using $(Symbol(pkg))
    end
end

using MKL
using Dates

start_time = time()

#"linear_solver" => "ma27",

function get_optimimizer()

    return optimizer_with_attributes(Ipopt.Optimizer, "sb" => "yes")

end

# Helper function to read and parse environment variables
function get_env(varname::String, default)
    val = get(ENV, varname, nothing)
    if isnothing(val)
        return default
    end
    # Try to convert to the type of the default value
    try
        return convert(typeof(default), parse(typeof(default), val))
    catch
        return val  # For string values, just return as is
    end
end

base_case_options = Dict(
    "sb" => "yes",
    "tol" => 1e-6,
    "mu_superlinear_decrease_power" => 1.25,
    "mu_linear_decrease_factor" => 0.4,
    "max_iter" => 500,
    "print_user_options" => "yes"
)

function read_ipopt_options(filename)
    options = Dict{String, Any}()
    for line in eachline(filename)
        if isempty(strip(line)) || startswith(strip(line), "#")
            continue  # skip empty lines or comments
        end
        k, v = split(strip(line))
        k = String(k)
        v = String(v)
        # Try to parse numbers, otherwise keep as string
        try
            v_parsed = parse(Int, v)
            options[k] = v_parsed
        catch
            try
                v_parsed = parse(Float64, v)
                options[k] = v_parsed
            catch
                options[k] = v
            end
        end
    end
    return options
end

if haskey(ENV, "BASE_CASE_OPTIONS")
   filename = ENV["BASE_CASE_OPTIONS"]
   if isfile(filename)
      println(" --- File $filename with base case options loaded! ---")
      merge!(base_case_options, read_ipopt_options(filename))  # Update existing dict      
   else
      println(" --- File $filename with base case options not found!! ---")
   end
end

function get_optimimizer_base_case()

    return optimizer_with_attributes(Ipopt.Optimizer, base_case_options...)

end

function get_optimimizer_base_case_recourse()

    return get_optimimizer_base_case()

end


cont_case_options = Dict(
                            "sb" => "yes",
                            "tol" =>  1e-6,
                            "mu_superlinear_decrease_power" =>  1.25,
                            "mu_linear_decrease_factor" =>  0.4,
                            "bound_relax_factor" =>  1e-6,
                            "max_iter" =>  500,
                            "fixed_variable_treatment" =>  "relax_bounds",
                            "jacobian_regularization_value" => 1e-10,
                            "inf_pr_output" =>  "internal",
                            "acceptable_dual_inf_tol" =>  0.01,
                            "acceptable_constr_viol_tol" => 1e-6,
                            "acceptable_compl_inf_tol" => 0.01,
                            "acceptable_iter" => 1,
                            "print_user_options"  =>  "yes"
)

if haskey(ENV, "CONTINGENCY_CASE_OPTIONS")
   filename = ENV["CONTINGENCY_CASE_OPTIONS"]
   if isfile(filename)
      println(" --- File $filename with contingency case options loaded! ---")
      merge!(cont_case_options, read_ipopt_options(filename))  # Update existing dict      
   else
      println(" --- File $filename with contingency case options not found!! ---")
   end
end

function get_optimimizer_contingency()

    return optimizer_with_attributes(Ipopt.Optimizer, cont_case_options...)

end


function pointer_manager()

    # IdDict to hold references
    cpp_jl_pointers = IdDict{Int64, Any}()

    # Protect an object, deleting any previous object with the same id
    function hold_pointer(obj, id)

        if haskey(cpp_jl_pointers, id)
            release_pointer(id)
        end
        cpp_jl_pointers[id] = obj

    end

    # Unprotect (release) an object by id
    function release_pointer(id)
        delete!(cpp_jl_pointers, id)
    end

    return hold_pointer
end

#hold_pointer = pointer_manager()

#--
#using Pkg; Pkg.add(["MPI", "Serialization", "DataFrames"])

if !haskey(ENV, "PATH_TO_EXAJUGO")
    println("'PATH_TO_EXAJUGO' not set!")
    exit(1)
end

push!(LOAD_PATH, string(ENV["PATH_TO_EXAJUGO"], "/modules"))
push!(LOAD_PATH, string(ENV["PATH_TO_EXAJUGO"], "/modules/SCACOPFSubproblems.jl"))


using CSV
using DataFrames
using Ipopt, JuMP
using SCACOPFSubproblems
using Serialization

include(string(ENV["PATH_TO_EXAJUGO"], "/modules/SCACOPFSubproblems/starting_point.jl"))
include(string(ENV["PATH_TO_EXAJUGO"], "/modules/corejugo/solution_evaluator.jl"))

println(" loaded package: ", ENV["PATH_TO_EXAJUGO"])

if haskey(ENV, "PATH_TO_TSSLOPE")

    include(string(ENV["PATH_TO_TSSLOPE"], "/init.jl"))
    println("PATH_TO_TSSLOPE is set to: ", ENV["PATH_TO_TSSLOPE"])

else

    println("PATH_TO_TSSLOPE is not set.")

end


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
   OPF_DATA =Ref(opt_data)

   return OPF_DATA
end

function load_ACOPF(rawfile, ropfile, confile)

   opt_data=SCACOPFdata(raw_filename=rawfile, rop_filename=ropfile, con_filename=confile)
   OPF_DATA =Ref(opt_data)
   return OPF_DATA
 
end

function load_ACOPF_instance(case)

    return load_ACOPF(get_instance_files(case)...)

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

# Default value
normalize_x = false
grad_multiplier = 1.0
hess_multiplier = 1.0


# Check if environment variable exists and update if so
if haskey(ENV, "NORMALIZE_X")
    # Convert string to Bool: accept "true", "1", "yes" (case-insensitive) as true
    val = lowercase(ENV["NORMALIZE_X"])
    normalize_x = val in ["true", "1", "yes"]
end

if haskey(ENV, "GRAD_MULTIPLIER")
    val = lowercase(ENV["GRAD_MULTIPLIER"])
    parsed_val = tryparse(Float64, val)
    if parsed_val !== nothing
        grad_multiplier = parsed_val
    end
end

if haskey(ENV, "HESS_MULTIPLIER")
    val = lowercase(ENV["HESS_MULTIPLIER"])
    parsed_val = tryparse(Float64, val)
    if parsed_val !== nothing
        hess_multiplier = parsed_val
    end
end

function getGradient(ptr, x)
   for (i,v) in enumerate(ptr[].cont_grad)
       x[i] = v
   end
   normalize_x ? (x .= normalize(x)) : nothing
end

function multiply_array(x, factor)

    x .*= factor

end


function array_to_struct(prob::Ref{SCACOPFdata}, arr, array_lengths::Ref{Dict{Symbol, Int}})

    MASTER_SOL= array_to_struct_generic(prob, arr, BasecaseSolution, array_lengths)
    return MASTER_SOL

end


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

function define_array_lengths(prob::Ref{SCACOPFdata})

    FIELD_SIZES_DICT= Ref(Dict(
        :v_n => size(prob[].N, 1),  # Same size as the number of rows in `prob.N`
        :theta_n => size(prob[].N, 1),  # Same size as the number of rows in `prob.N`
        :b_s => size(prob[].SSh, 1),  # Same size as the number of rows in `prob.SSh`
        :p_g => size(prob[].G, 1),  # Same size as the number of rows in `prob.N`
        :q_g => size(prob[].G, 1),  # Same size as the number of rows in `prob.N`
        :base_cost => 1,  # Assuming scalar size for `base_cost`
        :recourse_cost => 1,  # Assuming scalar size for `recourse_cost`
        :total_objective => 1  # Assuming scalar size for `recourse_cost`
    ))

    return FIELD_SIZES_DICT
end


#function full_solution_dim(fieldsizes)
function full_solution_dim(fieldsizes::Ref{Dict{Symbol, Int}})

    return sum(values(fieldsizes[]))

end

function getSolution(ptr, x)

   for (i,v) in enumerate(ptr[].p_g)
       x[i] = v
   end

end


function getObjective(ptr)

   #TODO: return the "whole" objective. see TODO in SCACOPFSubproblems.jl
   return ptr[].total_objective
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


#grad_multiplier and hess_multiplier are uaed for debugging
struct RecourseDerivatives

    gradient::Vector{Float64}
    hessian::Vector{Float64}
   # RecourseDerivativesAlloc(_grad, _hess) = new(_grad, _hess)
    function RecourseDerivatives(_grad, _hess, _len) 

        grad_copy = Vector{Float64}(undef, _len)
        hess_copy = Vector{Float64}(undef, _len)

        grad_copy .= grad_multiplier*_grad[1:_len]
        hess_copy .= hess_multiplier*_hess[1:_len]

        new(grad_copy, hess_copy)
  
    end
end

function get_recourse_derivatives(grad, hess, _len)
    return Ref(RecourseDerivatives(grad, hess, _len))
end

struct SparseMatrixIndexWrap

    rows::Vector{Int64}
    cols::Vector{Int64}

    function SparseMatrixIndexWrap(_rows, _cols, _len) 

        rows_copy = Vector{Int64}(undef, _len)
        cols_copy = Vector{Int64}(undef, _len)

        rows_copy .= _rows[1:_len]
        cols_copy .= _cols[1:_len]

        new(rows_copy, cols_copy)
  
    end
end

function get_sparse_matrix_index_wrap(rows, cols, _len)

    return Ref(SparseMatrixIndexWrap(rows, cols, _len))

end


struct SparseMatrixWrap

    indices::Ref{SparseMatrixIndexWrap}
    values::Vector{Float64}

    function SparseMatrixWrap(ix_ref::Ref{SparseMatrixIndexWrap}, _values, _len) 

        values_copy = Vector{Float64}(undef, _len)
        values_copy .= _values[1:_len]

        new(ix_ref, values_copy)
  
    end
end


function get_sparse_matrix_wrap(ix_ref, _values, _len)
    return Ref(SparseMatrixWrap(ix_ref, _values, _len))

end


struct RecourseSparseHessian
    gradient::Vector{Float64}
    hessian::Ref{SparseMatrixWrap}

    function RecourseSparseHessian(_grad, hess::Ref{SparseMatrixWrap}, _len) 

        grad_copy = Vector{Float64}(undef, _len)
        grad_copy .= _grad[1:_len]

        new(grad_copy, hess)
  
    end

end


function get_recourse_sparse(grad, hess, _len)

    return Ref(RecourseSparseHessian(grad, hess, _len))
end



using LinearAlgebra


function save_sparse_matrix(file_path, sparse_matrix, _iter)
   
    save_opt_sparse_data(file_path, sparse_matrix[], _iter)

end


using CSV, DataFrames, Printf

function save_opt_sparse_data(base_file_path::AbstractString, sparse_matrix::SparseMatrixWrap, iteration::Int)
    indices = sparse_matrix.indices[]
    rows = indices.rows
    cols = indices.cols
    vals = sparse_matrix.values

    n = length(vals)
    @assert length(rows) == n
    @assert length(cols) == n

    # Prepare DataFrame with id, row, col, value (5 decimals)
    new_rows = DataFrame(
        id = 0:(n-1),
        row = rows,
        col = cols,
        value = [@sprintf("%.5f", v) for v in vals]
    )

    # Build the new filename with iteration number
    parts = Base.splitext(base_file_path)
    new_file = string(parts[1], "_iter_", iteration, parts[2])

    CSV.write(new_file, new_rows)
end


function save_array(file_path, ptr, grad)
   
    save_opt_data(file_path, ("norm"=>norm(grad)), [round(x, digits=5) for x in grad], 1)

end


function save_cont_solution(file_path, ptr, prev_sol)

  try
    save_opt_data(file_path, ("objective"=>prev_sol[].cont_cost), [round(x, digits=5) for x in prev_sol[].p_g])

# Given file_path
    dir = dirname(file_path)
    base = basename(file_path)

    stem = splitext(base)[1]  # "solution_10"
    number = split(stem, "_")[end]  # "10"
    new_name = "iterations_$(number).csv" 

    new_path = joinpath(dir, new_name)

    save_opt_iterations(new_path, ("objective"=>prev_sol[].cont_cost))

  catch e
    println(" *** Exception ocorred saving coningency solution: ", e, " ***")
    flush(stdout)
  end

end

function save_solution(file_path, ptr, prev_sol)

  try

    save_opt_data(file_path, ("objective"=>prev_sol[].base_cost), [round(x, digits=5) for x in prev_sol[].p_g])

# Given file_path
    dir = dirname(file_path)
    base = basename(file_path)
    new_name = "iterations.csv" 
    new_path = joinpath(dir, new_name)

    save_opt_iterations(new_path, ("objective"=>prev_sol[].base_cost))
  catch e
    println(" *** Exception ocorred saving base case solution: ", e," ***")
    flush(stdout)
  end

end

function save_opt_data(file_path, kval::Pair{String, Float64}, sol, beg=0)

  fd_name = kval.first
  value = kval.second
  
  if isfile(file_path)
        # If the file exists, read it
        existing_data = CSV.read(file_path, DataFrame)
        
        # Determine the number of rows in the existing file
        num_rows = size(existing_data, 1)
        
        columns = Dict(
        :iteration => num_rows+beg,
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
        :iteration => beg,
        Symbol(fd_name) => value,
        [Symbol("x_$i") => sol[i] for i in 1:length(sol)]...)

        new_data = DataFrame(columns)
        CSV.write(file_path, new_data)
    end

end


function save_opt_iterations(file_path, kval::Pair{String, Float64})

  global start_time
  exec_time = time() - start_time  # End timer

  fd_name = kval.first
  value = kval.second

  # Define the column order explicitly
  col_names = [:iteration, Symbol(fd_name), :execution_time, :time_stamp]

  if isfile(file_path)
        # If the file exists, read it
        existing_data = CSV.read(file_path, DataFrame)
        
        # Determine the number of rows in the existing file
        num_rows = size(existing_data, 1)

        # Create a NamedTuple to ensure column order
        new_row = (; iteration=num_rows, Symbol(fd_name)=>value, execution_time=exec_time, time_stamp=now())
        new_data = DataFrame([new_row], col_names)

        # Append the new row to the existing data
        updated_data = vcat(existing_data, new_data)
        
        # Write the updated data back to the file
        CSV.write(file_path, updated_data)
    else
        # If the file does not exist, create it with iteration set to 0 and the given objective

# Create a single row as a NamedTuple in a vector
        row = (; iteration=0, Symbol(fd_name)=>value, execution_time=exec_time, time_stamp=now())
        new_data = DataFrame([row], col_names)
        CSV.write(file_path, new_data)
    end

end


function solve_base_case_recourse(ptr, prev_sol, ptr_rderivaties)

   global start_time
   start_time = time()

   G = ptr_rderivaties[].gradient
   H = ptr_rderivaties[].hessian

   n = length(G)

   # Add the quadratic approximation to the master problem, namely
   # q(x) = G'*(x-xk) + 0.5* (x-xk)'*H*(x-xk), where xk is the previous master solution
   # \nabla_q (x) = G + H(x-xk)
   # Hessian_q(x) = H

   xk = prev_sol[].p_g  #generation

   # CP? I do not know how to get xk: looks like it is in prev_sol. 

   # Create n scalar functions for each calculation

   ## the index 1 (in argH[1]) is because there each variable is entered separately because of the hessian
   recourse_fx = [ (x) -> begin  G[i]*(x-xk[i]) + 0.5*(x-xk[i])'*H[i]*(x-xk[i]) end for i in 1:n ]
   recourse_gx = [ (argG, x) -> begin argG[1] = G[i] + H[i]*(x-xk[i]) end for i in 1:n ]
   recourse_Hx = [ (argH, x) -> begin argH[1] = H[i]; end for i in 1:n ]


   SOLUTION_WITH_RECOURSE=
              Ref(solve_basecase(ptr[], get_optimimizer_base_case_recourse(), 
              recourse_f=recourse_fx, recourse_g=recourse_gx, recourse_H=recourse_Hx,
              previous_solution=prev_sol[])[1])

   # allocated_bytes = Base.gc_bytes() 
   #println("SOLUTION_WITH_RECOURSE Memory allocated: ", allocated_bytes, " bytes")

   return SOLUTION_WITH_RECOURSE
end

function solve_base_case_recourse_sparse(ptr, prev_sol, ptr_rderivaties)

   global start_time
   start_time = time()

   G = ptr_rderivaties[].gradient
   ptr_hess = ptr_rderivaties[].hessian
   idex_arrays = ptr_hess[].indices

   rows = idex_arrays[].rows
   cols = idex_arrays[].cols
   H = ptr_hess[].values

   n = length(G)

   xk = prev_sol[].p_g  #generation

   # Create n scalar functions for each calculation

   ## the index 1 (in argH[1]) is because there each variable is entered separately because of the hessian
   recourse_fx = [ (x) -> begin  G[i]*(x-xk[i]) + 0.5*(x-xk[i])'*H[i]*(x-xk[i]) end for i in 1:n ]
   recourse_gx = [ (argG, x) -> begin argG[1] = G[i] + H[i]*(x-xk[i]) end for i in 1:n ]
   recourse_Hx = [ (argH, x) -> begin argH[1] = H[i]; end for i in 1:n ]

   SOLUTION_WITH_RECOURSE=
              Ref(solve_basecase(ptr[], get_optimimizer_base_case_recourse(), 
              recourse_f=recourse_fx, recourse_g=recourse_gx, recourse_H=recourse_Hx,
              previous_solution=prev_sol[])[1])

   # allocated_bytes = Base.gc_bytes() 
   #println("SOLUTION_WITH_RECOURSE Memory allocated: ", allocated_bytes, " bytes")

   return SOLUTION_WITH_RECOURSE
end


function solve_base_case(ptr)

   global start_time
   start_time = time()

   SOLUTION_WITH_RECOURSE_BASE= Ref(solve_basecase(ptr[], get_optimimizer_base_case())[1])
  # SOLUTION_WITH_RECOURSE_BASE= Ref(solve_basecase(ptr[], get_optimimizer()))

   #allocated_bytes = Base.gc_bytes()
   #println("BASE Memory allocated: ", allocated_bytes, " bytes")

   return SOLUTION_WITH_RECOURSE_BASE

end


function get_instance_files(case)

    if haskey(ENV, "PATH_TO_INSTANCES")
        path_to_instances = ENV["PATH_TO_INSTANCES"]

    else
        println("'PATH_TO_INSTANCES' not set!")
        exit(1)
    end

    example_path = joinpath(path_to_instances, case, "")
  
    raw_file = example_path * "case.raw"
    rop_file = example_path * "case.rop"
    if haskey(ENV, "CONTINGENCY_FILE")
       confile = ENV["CONTINGENCY_FILE"]

    else
       confile = "case"
    end
    con_file = example_path * confile *".con"

    return raw_file, rop_file, con_file
end

function get_number_of_contingencies(case)

    return number_of_contingencies(load_ACOPF(get_instance_files(case)...))

end


function number_of_contingencies(ptr)

   return length(ptr[].cont_labels)

end

function number_of_columns(ptr)

   return nrow(ptr[].G)

end


function solve_contingency_pridec(ptr, i::Int64, ptr_basesol)
 
   global start_time
   start_time = time()

   ptr_basesol[].psd_hash = hash(ptr[])

   CONT_SOL = Ref(solve_contingency(ptr[], i, ptr_basesol[], get_optimimizer_contingency()))
   #debug: println(CONT_SOL[])
   return CONT_SOL

end


