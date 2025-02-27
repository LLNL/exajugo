# read command line parameters

INSTDIR = ARGS[1]
PYEVALDIR = ARGS[2]
TMPDIR = ARGS[3]
EMPTYCONFILE = joinpath(dirname(PROGRAM_FILE), "case.con")
EMPTYSOL2FILE = joinpath(dirname(PROGRAM_FILE), "solution2.txt")

# load modules

println("\nLoading modules ...")
clock = @elapsed begin
	# external modules
	using DataFrames, CSV, Ipopt, JuMP
	# internal modules
	push!(LOAD_PATH, joinpath(dirname(PROGRAM_FILE), "../modules"))
	using InstanceReader, SCACOPF, SolutionEvaluator, SolutionWriter
end
println("All modules loaded in ", round(clock, digits=2), " secs.")
flush(stdout)

println("\nDefining dependencies ...")
clock = @elapsed begin

# define function to execute evaluation script of the competition

function PyEvaluation(RAWFILE::String, ROPFILE::String, INLFILE::String, SOL1FILE::String)::Real

	# run python code for evaluation
	run(`python -W ignore $PYEVALDIR/test.py $RAWFILE $ROPFILE $EMPTYCONFILE $INLFILE $SOL1FILE $EMPTYSOL2FILE $TMPDIR/.summary.csv $TMPDIR/.details.csv`)
	
	# read solution value
	details = CSV.read(TMPDIR*"/.details.csv")
	
	# delete temporal file
	rm(TMPDIR*"/.details.csv")

	# return
	return details[!, :obj][1]

end

# define function to test correctness of ACOPF

function TestACOPF(RAWFILE::String, ROPFILE::String, INLFILE::String)

	# read instance
	println("\tReading and parsing instance data ...")
	clock = @elapsed begin
		MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches,
			switchedshunts, generatordsp, activedsptables, costcurves, governorresponse,
			contingencies = readinstance(RAWFILE, ROPFILE, INLFILE, EMPTYCONFILE)
		N, L, T, SSh, G, K, P = GOfmt2params(MVAbase, buses, loads, fixedbusshunts, generators,
			ntbranches, tbranches, switchedshunts, generatordsp, activedsptables, costcurves,
			governorresponse, contingencies)
	end
	println("\tData prepared for problem formulation in ", round(clock, digits=2), " secs.")
	flush(stdout)
	
        o = GoOptions()
	# solve AC OPF
	println("\tSolving AC OPF problem for the selected instance ...")
	clock = @elapsed begin
		SCACOPFSol =
			try
				#change the type of K.IDout from Vector{Union{Nothing, Int64}} to Vector{Int64}
				K.IDout = something.(K.IDout)
				SolveSCACOPF(o, N, L, T, SSh, G, K, P, SCACOPF.ipopt_opt())
			catch
				return NaN, NaN
			end
		implementationvalue = SCACOPFSol.prodcost
	end
	println("\tAC OPF problem solved in ", round(clock, digits=2), " secs.")
	flush(stdout)
	
	# write base case solution
	println("\tWritting AC OPF solution ...")
	clock = @elapsed begin
		writesolution(TMPDIR, SCACOPFSol.base_state.v_n, SCACOPFSol.base_state.theta_n, SCACOPFSol.base_state.b_s, 
		SCACOPFSol.base_state.p_g, SCACOPFSol.base_state.q_g, MVAbase, N, SSh, G, generators)
	end
	println("\tAC OPF solution written in ", round(clock, digits=2), " secs.")
	flush(stdout)
	
	# compute solution using GO competition evaluation code
	competitionvalue = PyEvaluation(RAWFILE, ROPFILE, INLFILE, TMPDIR*"/solution1.txt")
	
	# delete intermediate files
	rm(TMPDIR*"/solution1.txt")
	
	# return 
	return implementationvalue, competitionvalue
end

end
println("All dependencies defined in ", round(clock, digits=2), " secs.")
flush(stdout)

# run test over all instances within INSTDIR

println("\nTesting AC OPF for all instances within ", INSTDIR, " ...")
juliavalues = Float64[]
pythonvalues = Float64[]
for (root, dirs, files) in walkdir(INSTDIR)
	for file in files
		if file == "case.raw"
			println("\nTesting for instance at ", root)
			rawfile = joinpath(root, "case.raw")
			ropfile = joinpath(root, "case.rop")
			if !isfile(ropfile)
				ropfile = joinpath(root*"/..", "case.rop")
			end
			inlfile = joinpath(root, "case.inl")
			if !isfile(inlfile)
				inlfile = joinpath(root*"/..", "case.inl")
			end
                        # solve 
			jlvalue, pyvalue = TestACOPF(rawfile, ropfile, inlfile)

			println(root, "  Our objective: ", jlvalue, "    GO objective:  ", pyvalue)
			push!(juliavalues, jlvalue)
			push!(pythonvalues, pyvalue)
			flush(stdout)
		end
	end
end

# compute success ratio

reldiff = abs.(juliavalues .- pythonvalues) ./ max.(1., abs.(pythonvalues))
numsuccess = 0
numundetermined = 0
numfailed = 0
for i = 1:length(juliavalues)
	if isnan(juliavalues[i])
		global numundetermined += 1
	elseif reldiff[i] < 1E-3
		global numsuccess += 1
	else
		global numfailed += 1
	end
end
println("Success: ", numsuccess, "/", length(juliavalues))
println("Failed: ", numfailed, ", Undetermined: ", numundetermined)
