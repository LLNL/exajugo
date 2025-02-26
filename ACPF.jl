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

# function to read, solve power flow and write solution to a given location
# NOTE: power flow is solved as an optimization problem, setting the objective to penalize
#       deviations from the dispatch setpoints of generators, except for the generator at 
#       the swing bus

function ACPF(raw_filename::String, solution_dir::String)
	print("Reading instance at ", raw_filename, " ... ")
    psd = SCACOPFdata(raw_filename=raw_filename)
	println("done.\nSolving power flow problem ...")
	opt = optimizer_with_attributes(Ipopt.Optimizer,
		                            "sb" => "no")
	solution, summary = solve_base_power_flow(psd, opt)
	println("Power flow solved.")
    println("Total deviation from dispatch outside SWING buses: ", 
            round(summary[:p_deviations] * psd.MVAbase, digits=2), "MW")
    println("Total active nodal imbalance: ",
            round(summary[:active_nodal_imbalance] * psd.MVAbase, digits=2), "MW")
    println("Total reactive nodal imbalance: ",
            round(summary[:reactive_nodal_imbalance] * psd.MVAbase, digits=2), "MW")
    println("Total branch overloads: ",
            round(summary[:branch_overloads] * psd.MVAbase, digits=2), "MW")
    print("Writing solution to ", solution_dir, " ... ")
	if !ispath(solution_dir)
		mkpath(solution_dir)
	end
	write_solution(solution_dir, psd, solution)
	println("done.")
	return 
end

# if instancedir and solutiondir are given from cmd line -> run rectangular OPF
if length(ARGS) == 2
	ACPF(ARGS[1], ARGS[2]);
end
