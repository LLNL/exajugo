# modules
using Pkg;

PROGRAM_path = "/g/g91/roruiz/WF_ACOPF/exajugo_github/Project.toml"

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

function SCACOPF(instance_dir::String, solution_dir::String)
	println("Reading instance from "*instance_dir*" ... ")
    psd = SCACOPFdata(instance_dir)
	println("Done reading data.")
    solve_and_save_OPF(psd, solution_dir)
end

function SCACOPF(raw_filename::String, rop_filename::String, con_filename::String,
                     solution_dir::String)
	println("Reading instance from " * raw_filename * ", " * rop_filename *
			"and" * rop_filename *  " ... ")
    psd = SCACOPFdata(raw_filename=raw_filename, rop_filename=rop_filename, 
						con_filename = con_filename)
	println("Done reading data.")
    solve_and_save_OPF(psd, solution_dir)
end

function solve_and_save_OPF(psd::SCACOPFdata, solution_dir::String)
	if length(psd.K[:, :IDout]) == 0
		error("Either no contingencies were provided, or all contingencies were ignored.")
	end
	println("Solving SCACOPF using sparse OPF ...")
	opt = optimizer_with_attributes(Ipopt.Optimizer,
		                            "linear_solver" => "ma27",
		                            "sb" => "yes"
									)
    basecase_solution = solve_basecase(psd, opt)
    contingency = Vector{ContingencySolution}(undef, size(psd.K, 1))
    for k = 1:size(psd.K, 1)
        contingency[k] = solve_contingency(psd, k, basecase_solution, opt)
        contingency[k].cont_id = k
    end
	println("done.")
	return nothing
end

# if instancedir and solutiondir are given from cmd line -> run rectangular OPF
if length(ARGS) == 2
	SCACOPF(ARGS[1]);
elseif length(ARGS) == 4
	SCACOPF(ARGS[1], ARGS[2], ARGS[3]);
else
	error("Received ", length(ARGS), " input arguments, but expected 1 or 3.")
end
