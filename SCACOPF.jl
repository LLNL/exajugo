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
	solution = solve_SC_ACOPF(psd, opt)
	println("Done solving SCACOPF. \nWriting solution to "*solution_dir*" ... ")
	if !ispath(solution_dir)
		mkpath(solution_dir)
	end
	write_solution(solution_dir, psd, solution)
	println("done.")
	return nothing
end

# if instancedir and solutiondir are given from cmd line -> run rectangular OPF
if length(ARGS) == 2
	SCACOPF(ARGS[1], ARGS[2]);
elseif length(ARGS) == 4
	SCACOPF(ARGS[1], ARGS[2], ARGS[3], ARGS[4]);
else
	error("Received ", length(ARGS), " input arguments, but expected 2 or 4.")
end
