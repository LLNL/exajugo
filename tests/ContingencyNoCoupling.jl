# mimic current MyJulia2:
# 1. Solve ACOPF to obtain a base case solution
# 2. For each contingency (separately), with CouplingMode = :ignore
# 	2.a Solve contingency ignoring coupling
# 	2.b Treat contingency solution as base case solution and post-compute objective value
# 	2.c Verify that objective observed by Ipopt is the same as the post-computed

# read command line parameters

INSTDIR = ARGS[1]

# load modules

println("\nLoading modules ...")
clock = @elapsed begin
	# external modules
	using DataFrames, CSV, Ipopt, JuMP
	# internal modules
	push!(LOAD_PATH, joinpath(dirname(PROGRAM_FILE), "../modules"))
	using InstanceReader, SCACOPF, SolutionEvaluator
end
println("All modules loaded in ", round(clock, digits=2), " secs.")
flush(stdout)

println("\nDefining dependencies ...")
clock = @elapsed begin

# define function to test correctness of ACOPF

function TestContingencies(RAWFILE::String, ROPFILE::String, INLFILE::String, CONFILE::String)

	# read instance
	println("\tReading and parsing instance data ...")
	clock = @elapsed begin
		MVAbase, buses, loads, fixedbusshunts, generators, ntbranches, tbranches,
			switchedshunts, generatordsp, activedsptables, costcurves, governorresponse,
			contingencies = readinstance(RAWFILE, ROPFILE, INLFILE, CONFILE)
		N, L, T, SSh, G, K, P = GOfmt2params(MVAbase, buses, loads, fixedbusshunts, generators,
			ntbranches, tbranches, switchedshunts, generatordsp, activedsptables, costcurves,
			governorresponse, contingencies)
	end
	println("\tData prepared for problem formulation in ", round(clock, digits=2), " secs.")
	flush(stdout)

        o = GoOptions()
        #o.CouplingMode=:ignore

	# solve AC OPF
        println("\tSolving AC OPF problem for the selected instance ...")
        clock = @elapsed begin
	    SCACOPFSol =
		try
                    #solve SCACOPF with each type of contingencies
                    numOfEachCont=1
                    Ksubset = first(K[(K.ConType.==:Generator),:], numOfEachCont)
                    #append!(Ksubset, first(K[(K.ConType.==:Transformer),:], numOfEachCont))
                    #append!(Ksubset, first(K[(K.ConType.==:Line),:],        numOfEachCont))

                    # for no contingencies 
                    #Ksubset = K[Int[],:]
		    SCACOPFSol = SolveSCACOPF(o, N, L, T, SSh, G, Ksubset, P,
			SCACOPF.ipopt_opt())
		catch
		    error("failed to solve AC OPF for instance at ", RAWFILE)
		end
	end
	println("\tAC OPF problem solved in ", round(clock, digits=2), " secs.")
	flush(stdout)
	
	# solve all contingencies 
	#_, v_n0, theta_n0, b_s0, p_g0, q_g0, _ = SCACOPFSol
	numfailedipopt = 0
	numfailedtest = 0
	for k = 1:size(K, 1)
		clock = @elapsed begin
		println("\tTesting contingency ", k, " out of ", size(K, 1), "... ")
		
                o.SmoothingParamCoupling = 1e-2
		# solve contingency problem 
		#_, v_nk, theta_nk, b_sk, p_gk, q_gk, contingencypenalty =
                ContSol =
		    try
			SCACOPF.SolveContingencyBackend(o, N, L, T, SSh, G,
					                (K[!, :ConType][k], K[!, :IDout][k], K[!, :Contingency][k]),
					                P, SCACOPFSol, SCACOPF.ipopt_opt())
		    catch
			numfailedipopt += 1
			print("\tIpopt failed to solve contingency ", k, ".")
			flush(stdout)
			continue
		    end
		    
		# evaluate penalties using function for the base case
		if K[!, :ConType][k] == :Generator
			idxvec = collect(1:size(G,1))
			deleteat!(idxvec, Int(findfirst(G[!, :Generator] .== K[!, :IDout][k])))
			Gk = view(G, idxvec, :)
			p_gk_eval = view(ContSol.p_g, idxvec)
			q_gk_eval = view(ContSol.q_g, idxvec)
			Lk = L
			Tk = T
		elseif K[!, :ConType][k] == :Line
			Gk = G
			p_gk_eval = ContSol.p_g
			q_gk_eval = ContSol.q_g
			idxvec = collect(1:size(L,1))
			deleteat!(idxvec, Int(findfirst(L[!, :Line] .== K[!, :IDout][k])))
			Lk = view(L, idxvec, :)
			Tk = T
		elseif K[!, :ConType][k] == :Transformer
			Gk = G
			p_gk_eval = ContSol.p_g
			q_gk_eval = ContSol.q_g
			Lk = L
			idxvec = collect(1:size(T,1))
			deleteat!(idxvec, Int(findfirst(T[!, :Transformer] .== K[!, :IDout][k])))
			Tk = view(T, idxvec, :)
		else
			error("unrecognized contingency type ", K[!, :ConType][k])
		end
                contingencypenalty = ContSol.penalty
		#_, _, contingencypenaltyeval = BaseCaseValue(v_nk, theta_nk, b_sk,
		#	p_gk_eval, q_gk_eval, N, Lk, Tk, SSh, Gk, P)
  		_, _, contingencypenaltyeval = BaseCaseValue(ContSol.v_n, ContSol.theta_n, ContSol.b_s,
			p_gk_eval, q_gk_eval, N, Lk, Tk, SSh, Gk, P)
                    
		
		# test that penalties are equal
		if abs(contingencypenalty - contingencypenaltyeval)/max(1.0, contingencypenaltyeval) > 1E-8
			numfailedtest += 1
		end
		
		end
		println("\tTest for contingency ", k, " completed in ",
			round(clock, digits=2), " secs.\n",
			"\tIpopt penalty : ", contingencypenalty, "\n",
			"\tEval. penalty : ", contingencypenaltyeval)
		flush(stdout)
	end
	
	# return number of fails and total contingencies
	return size(K, 1), numfailedipopt, numfailedtest
	
end

end
println("All dependencies defined in ", round(clock, digits=2), " secs.")
flush(stdout)

# run test over all instances within INSTDIR

println("\nTesting contingency solve for all instances within ", INSTDIR, " ...")
numcon = Int[]
numfailedipopt = Int[]
numfailedtest = Int[]
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
			confile = joinpath(root, "case.con")
			numconinstance, numfailedipoptinstance, numfailedtestinstance =
				TestContingencies(rawfile, ropfile, inlfile, confile)
			println("\nFinished testing instance at ", root, "\n",
				"# contingencies : ", numconinstance, "\n",
				"# Ipopt failed  : ", numfailedipoptinstance, "\n",
				"# test failed   : ", numfailedtestinstance)
			push!(numcon, numconinstance)
			push!(numfailedipopt, numfailedipoptinstance)
			push!(numfailedtest, numfailedtestinstance)
			flush(stdout)
		end
	end
end

# compute success ratio
numcontingencies = sum(numcon)
numipoptfails = sum(numfailedipopt)
numtestfails = sum(numfailedtest)
numsuccess = numcontingencies - numipoptfails - numtestfails
println("Success      : ", numsuccess, "/", numcontingencies, "\n",
	"Failed Ipopt : ", numipoptfails, "/", numcontingencies, "\n",
	"Failed test  : ", numtestfails, "/", numcontingencies)
