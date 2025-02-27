#Driver for ACOPF without contingency using solver ma57 from Ipopt
# read command line parameters

#Directory that includes all required input files(raw,rop,inl)
#for more information on the input files, check InstanceReader module 
INSTDIR = ARGS[1] 
#Directory that output files will be 
TMPDIR = ARGS[2]
# Percent of bidding load out of total load
 RatioIn = ARGS[3]
 RatioBidding = tryparse(Float64, RatioIn)
 # number of timesteps
 TimeIn = ARGS[4]
 Tend = tryparse(Int64, TimeIn)
  # day index
  Day = tryparse(Int64, ARGS[5])
#  # method 
# MethodSol = ARGS[5]
#load empty con file since con file is required for some modules but contingency is not considered here
EMPTYCONFILE = joinpath(dirname(PROGRAM_FILE), "case.con")

# load modules

println("\nLoading modules ...")
clock = @elapsed begin
	# external modules
	using DataFrames, CSV, Ipopt, JuMP
	# internal modules
	push!(LOAD_PATH, joinpath(dirname(PROGRAM_FILE), "../modules"))
	using InstanceReaderMPbus, MPACOPFbus, SolutionWriterMPbus
	using SCACOPF, SolutionEvaluator, MPACOPF
end
println("All modules loaded in ", round(clock, digits=2), " secs.")
flush(stdout)

println("\nDefining dependencies ...")

#define function to solve ACOPF
function callMPACOPF(RAWFILE::String, ROPFILE::String, INLFILE::String, CTfile::String, plantfile::String)

	# read instance
	println("\tReading and parsing instance data ...")
	clock = @elapsed begin
   		MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches,
        	switchedshunts, generatordsp, activedsptables, costcurves, governorresponse,
        	contingencies, CT, plant = readinstance(RAWFILE, ROPFILE, INLFILE, EMPTYCONFILE, CTfile, plantfile)
	end
	 
	vN, vL, vT, vSSh, vG, vK, vP, timestep, vLoad, plant, PlantLoadall, DBratio, numPlants = GOfmt2paramsMPbus(MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches,
   	 switchedshunts, generatordsp, activedsptables, costcurves, governorresponse,contingencies,CT, plant, RatioBidding, Tend, Day)
	# return 
	println("\tData prepared for problem formulation in ", round(clock, digits=2), " secs.")
	flush(stdout)
	# return
    o = MPACOPFbus.GoOptions()
	# solve AC OPF
	println("\tSolving AC OPF problem for the selected instance ...")
	clock = @elapsed begin
		MPACOPFbusSol =
			try
                SolveMPACOPFbus(o, vN, vL, vT, vSSh, vG, vK, vP, timestep, vLoad, plant, MVAbase, MPACOPF.ipopt_optPrint())
			catch err
				println(err)
				rethrow(err)
				return NaN
			end
	end
	println("\tAC OPF problem solved in ", round(clock, digits=2), " secs.")
	flush(stdout)

    # return
	# write solution
	println("\tWritting AC OPF solution ...")
	clock = @elapsed begin
		 plantpower, LMPs, Gen, GenCost = MPsolstate(MPACOPFbusSol)
		writesolutionAllBuses(TMPDIR, timestep, plantpower, LMPs, Gen, GenCost, numPlants, Day, RatioBidding)

	end
	println("\tAC OPF solution written in ", round(clock, digits=2), " secs.")
	flush(stdout)
	
end

## function to get outputs from the solution

function MPsolstate(MPACOPFbusSol::MPACOPFbus.MPbusSol)

	plantpower = MPACOPFbusSol.plantpowerAll
	LMPs = MPACOPFbusSol.LMPsbusSingle

	Gtot = MPACOPFbusSol.GenAll
	gencoststep = MPACOPFbusSol.GenCostStep

	return  plantpower, LMPs, Gtot, gencoststep
end

#run the test with required input files within INSTDIR
obj_value = Float64[]
rawfile = joinpath(INSTDIR, "case.raw")
ropfile = joinpath(INSTDIR, "case.rop")
inlfile = joinpath(INSTDIR, "case.inl")
CTfile = joinpath(INSTDIR, "Load.CSV")
# CTfile = joinpath(INSTDIR, "case.CSV")
plantfile = joinpath(INSTDIR, "plant.csv")
#solve 
obj_value = callMPACOPF(rawfile, ropfile, inlfile, CTfile, plantfile)