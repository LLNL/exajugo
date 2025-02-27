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
 # method 
MethodSol = ARGS[5]
#load empty con file since con file is required for some modules but contingency is not considered here
EMPTYCONFILE = joinpath(dirname(PROGRAM_FILE), "case.con")

# load modules

println("\nLoading modules ...")
clock = @elapsed begin
	# external modules
	using DataFrames, CSV, Ipopt, JuMP
	# internal modules
	push!(LOAD_PATH, joinpath(dirname(PROGRAM_FILE), "../modules"))
	using InstanceReaderMP, MPACOPF, SolutionWriterMP
	using SCACOPF, SolutionEvaluator
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
	 
	vN, vL, vT, vSSh, vG, vK, vP, timestep, plant, plant_unsort, DBratio = GOfmt2paramsMP(MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches,
   	 switchedshunts, generatordsp, activedsptables, costcurves, governorresponse,contingencies,CT, plant, RatioBidding, Tend)
	#  writepara(TMPDIR, vG[1],plant_unsort)
	# return 
	println("\tData prepared for problem formulation in ", round(clock, digits=2), " secs.")
	flush(stdout)
	# return
    o = MPACOPF.GoOptions()
	# solve AC OPF
	println("\tSolving AC OPF problem for the selected instance ...")
	clock = @elapsed begin
		MPACOPFSol =
			try
				if MethodSol == "PSQP"
					SolveMasterSQP(o, vN, vL, vT, vSSh, vG, vK, vP, timestep, plant, MVAbase, MPACOPF.ipopt_optPrint(), MPACOPF.ipopt_opt())
				elseif MethodSol == "Preg"
					SolveMaster(o, vN, vL, vT, vSSh, vG, vK, vP, timestep, plant, MVAbase, MPACOPF.ipopt_QN(), MPACOPF.ipopt_opt())
				else
					SolveMPACOPFsingle(o, vN, vL, vT, vSSh, vG, vK, vP, timestep, plant, MVAbase, MPACOPF.ipopt_optPrint())
				end

			catch err
				println(err)
				rethrow(err)
				return NaN
			end
	end
	println("\tAC OPF problem solved in ", round(clock, digits=2), " secs.")
	flush(stdout)

    return
	# write base case solution
	println("\tWritting AC OPF solution ...")
	clock = @elapsed begin
		v_n, theta_n, b_s, p_g, q_g, prodcost, basepen, plantrev, plantpower, LMPs, LMPsN,
		pmslack, ppslack, qmslack, qpslack, lslack, tslack, Lloss, Gtot, eta, pline, ptran = MPsolstate(MPACOPFSol)
		writesolutionBid(TMPDIR, MVAbase, vN, vSSh, vG, generators, plant, timestep, v_n, 
			theta_n, b_s, p_g, q_g, prodcost, basepen, plantrev, plantpower, LMPs, LMPsN, DBratio, RatioBidding) 
		writesolutionPlant(TMPDIR, vN, timestep, plantpower, plant, LMPs, DBratio, RatioBidding)
		# writesolutionSlack(TMPDIR, pmslack, ppslack, qmslack, qpslack, lslack, tslack, vN, vL, vT, timestep, RatioBidding)
		writeloss(TMPDIR, Lloss, Gtot, eta, pline, ptran, vL, vT, timestep, RatioBidding)

	end
	println("\tAC OPF solution written in ", round(clock, digits=2), " secs.")
	flush(stdout)
	
end

## function to get outputs from the solution

function MPsolstate(MPACOPFSol::MPACOPF.MPACOPFBaseState)
	v_n = MPACOPFSol.base_results.v_n
	theta_n = MPACOPFSol.base_results.theta_n
	b_s = MPACOPFSol.base_results.b_s
	p_g = MPACOPFSol.base_results.p_g
	q_g = MPACOPFSol.base_results.q_g

	prodcost = MPACOPFSol.bidding_results.prodcost
	basepen = MPACOPFSol.bidding_results.penalty_base
	plantrev = MPACOPFSol.bidding_results.plantrev
	plantpower = MPACOPFSol.bidding_results.plantpower
	LMPs = MPACOPFSol.bidding_results.LMPs
	LMPsN = MPACOPFSol.bidding_results.LMPsN

	pmslack = MPACOPFSol.slackn.spm
	ppslack = MPACOPFSol.slackn.spp
	qmslack = MPACOPFSol.slackn.sqm
	qpslack = MPACOPFSol.slackn.sqp
	lslack = MPACOPFSol.slackn.sl
	tslack = MPACOPFSol.slackn.st

	Lloss = MPACOPFSol.TLoss.LineLoss
	Gtot = MPACOPFSol.TLoss.GenTotal
	eta = MPACOPFSol.TLoss.Eta
	pline = MPACOPFSol.TLoss.p_line
	ptran = MPACOPFSol.TLoss.p_tran

	return v_n, theta_n, b_s, p_g, q_g, prodcost, basepen, plantrev, plantpower, LMPs, LMPsN,
		pmslack, ppslack, qmslack, qpslack, lslack, tslack, Lloss, Gtot, eta, pline, ptran
end

## function to plot 



#run the test with required input files within INSTDIR
obj_value = Float64[]
rawfile = joinpath(INSTDIR, "case.raw")
ropfile = joinpath(INSTDIR, "case.rop")
inlfile = joinpath(INSTDIR, "case.inl")
CTfile = joinpath(INSTDIR, "case.CSV")
plantfile = joinpath(INSTDIR, "plant.csv")
#solve 
obj_value = callMPACOPF(rawfile, ropfile, inlfile, CTfile, plantfile)