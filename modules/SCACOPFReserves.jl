# function to compute ad-hoc generation reserves for base case
# NOTE: this is very conservative, it guarantees almost no saturation in AGC

function GenerationReserves(G::DataFrame, K::DataFrame)::Vector{Float64}
	
	# get list of generators that fail
	ConGen = Int[]
	for k = 1:size(K, 1)
		if K[!, :ConType][k] == :Generator
			push!(ConGen, K[!, :IDout][k])
		end
	end
	ConGidx = convert(Vector{Int}, indexin(ConGen, G[!, :Generator]))
	
	# compute reserve coefficient
	alphap = max.(0.0, G[!, :alpha])
	Pubp = max.(0.0, G[!, :Pub])
	sumalphap = sum(alphap)
	rescoeff = 0.0
	for g = ConGidx
		rcg = Pubp[g]/(sumalphap - alphap[g])
		if rcg > rescoeff
			rescoeff = rcg
		end
	end
	
	# return vector of reserve margin per generator
	return min.(rescoeff*alphap, G[!, :Pub]-G[!, :Plb])
	
end
