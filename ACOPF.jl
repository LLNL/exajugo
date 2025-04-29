# modules
using Pkg;
if dirname(PROGRAM_FILE) == ""
    Pkg.activate(".")
	push!(LOAD_PATH, "./modules")
else
    Pkg.activate(dirname(PROGRAM_FILE))
	push!(LOAD_PATH, string(dirname(PROGRAM_FILE), "/modules"))
end
using Ipopt, JuMP, Printf
using SCACOPFSubproblems

# function to read, solve rectangular OPF and write solution to a given location

function ACOPF(instance_dir::String, solution_dir::String;
				system_dir::Union{Nothing, String} = nothing)
	println("Reading ACOPF instance from "*instance_dir*" ... ")
    psd = SCACOPFdata(instance_dir)
	println("Done reading data.")
    solve_and_save_OPF(psd, solution_dir, system_dir = system_dir)
end

function ACOPF(raw_filename::String, rop_filename::String, solution_dir::String;
				system_dir::Union{Nothing, String} = nothing)
	println("Reading ACOPF instance from " * raw_filename * " and " * rop_filename * " ... ")
    psd = SCACOPFdata(raw_filename=raw_filename, rop_filename=rop_filename)
	println("Done reading data.")
    solve_and_save_OPF(psd, solution_dir, system_dir = system_dir)
end

function solve_and_save_OPF(psd::SCACOPFdata, solution_dir::String;
							system_dir::Union{Nothing, String} = nothing)
	println("Solving basecase using sparse OPF ...")
	opt = optimizer_with_attributes(Ipopt.Optimizer,
		                            #"linear_solver" => "ma57",
		                            "sb" => "yes")
	solution = solve_basecase(psd, opt, output_dir = system_dir)
	println("Done solving OPF. Objective value: \$", round(solution.base_cost, digits=1),
		".\nWriting solution to "*solution_dir*" ... ")
	if !ispath(solution_dir)
		mkpath(solution_dir)
	end
	write_solution(solution_dir, psd, solution)
	println("done.")
	return nothing
end

# if instancedir and solutiondir are given from cmd line -> run rectangular OPF
if length(ARGS) == 2
	ACOPF(ARGS[1], ARGS[2]);
elseif length(ARGS) == 3
	# if raw file, rop file and solutiondir are given from cmd line -> run rectangular OPF
	if endswith(ARGS[2], ".rop")
		ACOPF(ARGS[1], ARGS[2], ARGS[3]);
	else
		# if instancedir, solutiondir and systemdir are given from cmd line -> run rectangular OPF
		ACOPF(ARGS[1], ARGS[2], system_dir =  ARGS[3]);
	end
# if raw file, rop file, solutiondir and systemdir are given from cmd line -> run rectangular OPF
elseif length(ARGS) == 4
	ACOPF(ARGS[1], ARGS[2], ARGS[3], system_dir = ARGS[4]);
else
	error("Received ", length(ARGS), " input arguments, but expected 2-4.")
end
