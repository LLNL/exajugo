# read command line parameters

#Directory that includes all required file
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
    using InstanceReaderMP, MPACOPF, SCACOPF
end
println("All modules loaded in ", round(clock, digits=2), " secs.")
flush(stdout)

CTfile = joinpath(INSTDIR, "case.CSV")
plantfile = joinpath(INSTDIR, "plant.CSV")
RAWFILE = joinpath(INSTDIR, "case.raw")
ROPFILE = joinpath(INSTDIR, "case.rop")
INLFILE = joinpath(INSTDIR, "case.inl")

# read instance
println("\tReading and parsing instance data ...")
clock = @elapsed begin
    MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches,
        switchedshunts, generatordsp, activedsptables, costcurves, governorresponse,
        contingencies, CT, plant = readinstance(RAWFILE, ROPFILE, INLFILE, EMPTYCONFILE, CTfile, plantfile)
end
vN, vL, vT, vSSh, vG, vK, vP, timestep, plant = GOfmt2paramsMP(MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches,
    switchedshunts, generatordsp, activedsptables, costcurves, governorresponse,contingencies,CT, plant)


println("\tSolving AC OPF problem for the selected instance ...")
o = MPACOPF.GoOptions()

SolveMPACOPF(o, vN, vL, vT, vSSh, vG, vK, vP, timestep, plant, SCACOPF.ipopt_opt())

