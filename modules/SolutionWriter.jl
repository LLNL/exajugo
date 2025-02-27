#__precompile__()

module SolutionWriter

## external modules used

using SparseArrays, Printf, DataFrames

## elements to be exported

export writesolution

## function to write solution block

function writesolutionblock(io::IO, v_n::T, theta_n::T, b_s::T, p_g::T, q_g::T,
	MVAbase::Float64, N::DataFrame, SSh::DataFrame, G::DataFrame,
	generators::DataFrame)::Nothing where {T <: AbstractVector{Float64}}
	
	# compute reactive compensaton at each bus
	SSh_Nidx = convert(Vector{Int}, indexin(SSh[!, :Bus], N[!, :Bus]))
	bcsn = spzeros(size(N, 1))
	for ssh = 1:size(SSh, 1)
		bcsn[SSh_Nidx[ssh]] += b_s[ssh]
	end
	bcsn *= MVAbase
	
	# write bus section
	@printf(io, "--bus section\ni, v(p.u.), theta(deg), bcs(MVAR at v = 1 p.u.)\n")
	for n = 1:size(N, 1)
		@printf(io, "%d, %.20f, %.20f, %.20f\n", N[!, :Bus][n], v_n[n], 180/pi*theta_n[n], bcsn[n])
	end
	
	# write generator section
	gmap = zeros(Int, size(generators, 1))
	for g = 1:size(G, 1)
		gmap[G[!, :Generator][g]] = g
	end
	@printf(io, "--generator section\ni, id, p(MW), q(MW)\n")
	for gi = 1:size(generators, 1)
		if gmap[gi] == 0
			@printf(io, "%d, \'%s\', 0, 0\n", generators[!, :I][gi], generators[!, :ID][gi])
		else
			g = gmap[gi]
			@printf(io, "%d, \'%s\', %.12f, %.12f\n", G[!, :Bus][g], G[!, :BusUnitNum][g],
				MVAbase*p_g[g], MVAbase*q_g[g])
		end
	end
	
	return
	
end

## function to write solution (solution1.txt and solution2.txt)

function writesolution(OutDir::String, v_n::T, theta_n::T, b_s::T, p_g::T, q_g::T,
	MVAbase::Float64, N::DataFrame, SSh::DataFrame, G::DataFrame, generators::DataFrame,
	WriteContingencies::Bool=false, v_nk::U=Array{Float64}(undef, 0, 0),
	theta_nk::U=Array{Float64}(undef, 0, 0), b_sk::U=Array{Float64}(undef, 0, 0),
	p_gk::U=Array{Float64}(undef, 0, 0), q_gk::U=Array{Float64}(undef, 0, 0),
	delta_k::T=Float64[], K::DataFrame=DataFrame(), contingencies::DataFrame=DataFrame())::
	Nothing where {T <: AbstractVector{Float64}, U <: AbstractMatrix{Float64}}
	
	# write base case solution
	f = open(OutDir * "/solution1.txt", "w")
	writesolutionblock(f, v_n, theta_n, b_s, p_g, q_g, MVAbase, N, SSh, G, generators)
	close(f)
	
	# return here if there is nothing else to do
	if !WriteContingencies
		return
	end
	
	# write contingency solutions
	kmap = zeros(Int, size(contingencies, 1))
	for k = 1:size(K, 1)
		kmap[K[!, :Contingency][k]] = k
	end
	f = open(OutDir * "/solution2.txt", "w")
	for ki = 1:size(contingencies, 1)
		if kmap[ki] != 0
			@printf(f, "--contingency\nlabel\n")
			@printf(f, "\'%s\'\n", contingencies[!, :LABEL][ki])
			k = kmap[ki]
			writesolutionblock(f, view(v_nk, :, k), view(theta_nk, :, k),
				view(b_sk, :, k), view(p_gk, :, k), view(q_gk, :, k),
				MVAbase, N, SSh, G, generators)
			@printf(f, "--delta section\ndelta(MW)\n%g\n", MVAbase*delta_k[k])
		else
			error("no matching result for contingency ", contingencies[!, :LABEL][ki])
		end
	end
	close(f)
	
end

end
