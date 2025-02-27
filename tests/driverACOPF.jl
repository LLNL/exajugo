#Driver for ACOPF without contingency using solver ma57 from Ipopt
# read command line parameters

#Directory that includes all required input files(raw,rop,inl)
#for more information on the input files, check InstanceReader module 
INSTDIR = ARGS[1] 
#Directory that output files will be (tbd)
TMPDIR = ARGS[2]
#load empty con file since con file is required for some modules but contingency is not considered here
EMPTYCONFILE = joinpath(dirname(PROGRAM_FILE), "case.con")

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

#define function to solve ACOPF
function SolveACOPF(RAWFILE::String, ROPFILE::String, INLFILE::String)

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
	
    o = SCACOPF.GoOptions()
	# solve AC OPF
	println("\tSolving AC OPF problem for the selected instance ...")
	clock = @elapsed begin
		SCACOPFSol =
			try
                #change the type of K.IDout from Vector{Union{Nothing, Int64}} to Vector{Int64}
				K.IDout = something.(K.IDout)
				#input 2: "quadratic" for approximated quadratic penalties, "piecewise_linear" for GO piecewise linear penalties
				# SCACOPF.SolveSCACOPF(o, "piecewise_linear", N, L, T, SSh, G, K, P, Ipopt.Optimizer)
				SolveSCACOPF(o, N, L, T, SSh, G, K, P, SCACOPF.ipopt_opt())
			catch err
				println(err)
				rethrow(err)
				return NaN, NaN
			end
	end
	println("\tAC OPF problem solved in ", round(clock, digits=2), " secs.")
	flush(stdout)
	
end

#run the test with required input files within INSTDIR
obj_value = Float64[]
rawfile = joinpath(INSTDIR, "case.raw")
ropfile = joinpath(INSTDIR, "case.rop")
inlfile = joinpath(INSTDIR, "case.inl")
#solve 
obj_value = SolveACOPF(rawfile, ropfile, inlfile)