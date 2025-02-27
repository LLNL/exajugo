# read command line parameters

TESTDIR = ARGS[1]
OUTFNAME = ARGS[2]

# load modules

println("\nLoading modules ...")
clock = @elapsed begin
	# external modules
	using LinearAlgebra
	# internal modules
	push!(LOAD_PATH, joinpath(dirname(PROGRAM_FILE), "../modules"))
	using InstanceReader, SolutionReader, SolutionEvaluator
end
println("All modules loaded in ", round(clock, digits=2), " secs.")
flush(stdout)

# function to compute key performance indicators of test

function computeindicators(testdir::String)
	N, L, T, SSh, G, K, P = ParseInstance(testdir)
	L[:RateBase] .*= .75
	T[:RateBase] .*= 1.
	v_n, theta_n, b_s, p_g, q_g = ReadBaseSolution(testdir, 100.0, N, SSh, G)
	p_li, q_li, p_ti, q_ti, pslackm_n, pslackp_n, qslackm_n, qslackp_n, sslack_li, sslack_ti =
		BaseCaseSolution(v_n, theta_n, b_s, p_g, q_g, N, L, T, SSh, G)
	_, prodcost, penalties = BaseSolutionValue(p_g, pslackm_n, pslackp_n, qslackm_n, qslackp_n,
		sslack_li, sslack_ti, G, P)
	p_li_cpp, q_li_cpp, p_ti_cpp, q_ti_cpp, pslackm_n_cpp, pslackp_n_cpp,
		qslackm_n_cpp, qslackp_n_cpp, sslack_li_cpp, sslack_ti_cpp =
		ReadBaseSolutionExtras(testdir, 100.0, N, L, T)
	deltaN = (norm(pslackm_n - pslackm_n_cpp) + norm(pslackp_n - pslackp_n_cpp) +
		norm(qslackm_n - qslackm_n_cpp) + norm(qslackp_n - qslackp_n_cpp))/size(N, 1)
	deltaL = (norm(p_li - p_li_cpp) + norm(q_li - q_li_cpp) + norm(sslack_li - sslack_li_cpp))/size(L, 1)
	deltaT = (norm(p_ti - p_ti_cpp) + norm(q_ti - q_ti_cpp) + norm(sslack_ti - sslack_ti_cpp))/size(T, 1)
	return prodcost, penalties, deltaN, deltaL, deltaT
end

# create execution scripts for all instances within INSTDIR

println("\nComputing perfomance indicators for all tests within ", TESTDIR, ".")
f = open(OUTFNAME, "w")
println(f, "testdir,basecost,basepenalty,deltaN,deltaL,deltaT")
for (root, dirs, files) in walkdir(TESTDIR)
	for file in files
		if file == "solution1_extras.txt"
			prodcost, penalties, deltaN, deltaL, deltaT = computeindicators(root)
			println(f, root, ",", prodcost, ",", penalties, ",", deltaN, ",",
				deltaL, ",", deltaT)
		end
	end
end
close(f)
println("All done.")
