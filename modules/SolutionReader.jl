#__precompile__()

module SolutionReader

## external modules used

using DataFrames

## internal modules used

## elements to be exported

export ReadBaseSolution, ReadBaseSolutionExtras

## function to read solution block

function readsolutionblock(io::IO)
	
	# check that first line is bus header
	@assert readline(io) == "--bus section"
	readline(io);	# ignoring headers
	
	# read bus section
	I_n = Int[]
	v_n = Float64[]
	theta_n = Float64[]
	b_n = Float64[]
	l = readline(io)
	while l != "--generator section"
		fields = split(l, ", ")
		push!(I_n, parse(Int, fields[1]))
		push!(v_n, parse(Float64, fields[2]))
		push!(theta_n, parse(Float64, fields[3]))
		push!(b_n, parse(Float64, fields[4]))
		l = readline(io)
	end
	readline(io);	# ignoring headers
	
	# read generator section
	I_g = Int[]
	ID_g = String[]
	p_g = Float64[]
	q_g = Float64[]
	l = readline(io)
	while l != "" && l[1:2] != "--"
		fields = split(l, ", ")
		push!(I_g, parse(Int, fields[1]))
		push!(ID_g, replace(fields[2], "'" => ""))
		push!(p_g, parse(Float64, fields[3]))
		push!(q_g, parse(Float64, fields[4]))
		l = readline(io)
	end

	return I_n, v_n, theta_n, b_n, I_g, ID_g, p_g, q_g
	
end

## function to read base case solution

function ReadBaseSolution(OutDir::String, MVAbase::Float64, N::DataFrame, SSh::DataFrame, G::DataFrame)

	# read base solution file
	io = open(joinpath(OutDir, "solution1.txt"), "r")
	I__n, v__n, theta__n, b__n, I__g, ID__g, p__g, q__g = readsolutionblock(io)
	close(io)

	# parse solution v__n, theta__n to vectors v_n, theta_n
	theta__n *= pi/180
	Nidx = convert(Vector{Int}, indexin(N[:Bus], I__n))
	v_n = Vector{Float64}(undef, size(N, 1))
	theta_n = Vector{Float64}(undef, size(N, 1))
	for n = 1:size(N, 1)
		v_n[n] = v__n[Nidx[n]]
		theta_n[n] = theta__n[Nidx[n]]
	end

	# parse solution b__n to vector b_s
	b__n /= MVAbase
	SSh_Nidx = convert(Vector{Int}, indexin(SSh[:Bus], N[:Bus]))
	b_s = Vector{Float64}(undef, size(SSh, 1))
	for ssh = 1:size(SSh, 1)
		bn = b__n[SSh_Nidx[ssh]]
		if bn < SSh[:Blb][ssh]
			b_s[ssh] = SSh[:Blb][ssh]
		elseif bn > SSh[:Bub][ssh]
			b_s[ssh] = SSh[:Bub][ssh]
		else
			b_s[ssh] = bn
		end
		b__n[SSh_Nidx[ssh]] -= b_s[ssh]
	end
	if sum(abs.(b__n)) > 1E-4
		@warn(string("there are ", round(MVAbase*sum(abs.(b__n)), digits=4),
			"MVAR unassigned to shunts."))
	end

	# parse solution p__g, q__g to vectors p_g, q_g
	p__g /= MVAbase
	q__g /= MVAbase
	Gidx = convert(Vector{Int}, indexin(string.(G[:Bus], ":", G[:BusUnitNum]),
		string.(I__g, ":", ID__g)))
	p_g = Vector{Float64}(undef, size(G, 1))
	q_g = Vector{Float64}(undef, size(G, 1))
	for g = 1:size(G, 1)
		p_g[g] = p__g[Gidx[g]]
		q_g[g] = q__g[Gidx[g]]
	end
	
	return v_n, theta_n, b_s, p_g, q_g
	
end

## function to read solution extras block

function readsolutionextrasblock(io::IO)
	
	# check that first line is bus header
	@assert readline(io) == "--bus section"
	
	# read bus section
	readline(io);	# ignoring headers
	I_n = Int[]
	pslack_n = Float64[]
	qslack_n = Float64[]
	l = readline(io)
	while l != "" && l[1:2] != "--"
		fields = split(l, ",")
		push!(I_n, parse(Int, fields[1]))
		push!(pslack_n, parse(Float64, fields[2]))
		push!(qslack_n, parse(Float64, fields[3]))
		l = readline(io)
	end
	
	# read line section
	readline(io);	# ignoring headers
	I_l = Int[]
	J_l = Int[]
	CKT_l = String[]
	p1_l = Float64[]
	q1_l = Float64[]
	sslack1_l = Float64[]
	p2_l = Float64[]
	q2_l = Float64[]
	sslack2_l = Float64[]
	l = readline(io)
	while l != "" && l[1:2] != "--"
		fields = split(l, ",")
		push!(I_l, parse(Int, fields[1]))
		push!(J_l, parse(Int, fields[2]))
		push!(CKT_l, fields[3])
		push!(p1_l, parse(Float64, fields[4]))
		push!(q1_l, parse(Float64, fields[5]))
		push!(sslack1_l, parse(Float64, fields[6]))
		push!(p2_l, parse(Float64, fields[7]))
		push!(q2_l, parse(Float64, fields[8]))
		push!(sslack2_l, parse(Float64, fields[9]))
		l = readline(io)
	end
	
	# read transformer section
	readline(io);	# ignoring headers
	I_t = Int[]
	J_t = Int[]
	CKT_t = String[]
	p1_t = Float64[]
	q1_t = Float64[]
	sslack1_t = Float64[]
	p2_t = Float64[]
	q2_t = Float64[]
	sslack2_t = Float64[]
	l = readline(io)
	while l != "" && l[1:2] != "--"
		fields = split(l, ",")
		push!(I_t, parse(Int, fields[1]))
		push!(J_t, parse(Int, fields[2]))
		push!(CKT_t, fields[3])
		push!(p1_t, parse(Float64, fields[4]))
		push!(q1_t, parse(Float64, fields[5]))
		push!(sslack1_t, parse(Float64, fields[6]))
		push!(p2_t, parse(Float64, fields[7]))
		push!(q2_t, parse(Float64, fields[8]))
		push!(sslack2_t, parse(Float64, fields[9]))
		l = readline(io)
	end
	
	return (I_n, pslack_n, qslack_n),
		(I_l, J_l, CKT_l, p1_l, q1_l, sslack1_l, p2_l, q2_l, sslack2_l),
		(I_t, J_t, CKT_t, p1_t, q1_t, sslack1_t, p2_t, q2_t, sslack2_t)
	
end

## function to read base case solution extras

function ReadBaseSolutionExtras(OutDir::String, MVAbase::Float64,
	N::DataFrame, L::DataFrame, T::DataFrame)
	
	# read base solution extras file
	io = open(joinpath(OutDir, "solution1_extras.txt"), "r")
	busdata, linedata, trafodata = readsolutionextrasblock(io)
	close(io)
	
	# parse pslack_n, qslack_n
	Nidx = indexin(N[:Bus], busdata[1])
	@assert all(Nidx .!= nothing)
	pslackm_n = Vector{Float64}(undef, size(N, 1))
	pslackp_n = Vector{Float64}(undef, size(N, 1))
	qslackm_n = Vector{Float64}(undef, size(N, 1))
	qslackp_n = Vector{Float64}(undef, size(N, 1))
	for n = 1:size(N, 1)
		if busdata[2][Nidx[n]] > 0
			pslackp_n[n] = busdata[2][Nidx[n]]
			pslackm_n[n] = 0.0
		else
			pslackp_n[n] = 0.0
			pslackm_n[n] = -busdata[2][Nidx[n]]
		end
		if busdata[3][Nidx[n]] > 0
			qslackp_n[n] = busdata[3][Nidx[n]]
			qslackm_n[n] = 0.0
		else
			qslackp_n[n] = 0.0
			qslackm_n[n] = -busdata[3][Nidx[n]]
		end
	end
	
	# parse p_li, q_li, sslack_li
	Lidx = indexin(string.(L[:From], ":", L[:To], ":", L[:CktID]),
		string.(linedata[1], ":", linedata[2], ":", linedata[3]))
	@assert all(Lidx .!= nothing)
	p_li = Array{Float64, 2}(undef, size(L, 1), 2)
	q_li = Array{Float64, 2}(undef, size(L, 1), 2)
	sslack_li = Array{Float64, 2}(undef, size(L, 1), 2)
	for l = 1:size(L, 1), i=1:2
		p_li[l,i] = linedata[4+(i-1)*3][Lidx[l]]
		q_li[l,i] = linedata[5+(i-1)*3][Lidx[l]]
		sslack_li[l,i] = linedata[6+(i-1)*3][Lidx[l]]
	end
	
	# parse p_ti, q_ti, sslack_ti
	Tidx = indexin(string.(T[:From], ":", T[:To], ":", T[:CktID]),
		string.(trafodata[1], ":", trafodata[2], ":", trafodata[3]))
	@assert all(Tidx .!= nothing)
	p_ti = Array{Float64, 2}(undef, size(T, 1), 2)
	q_ti = Array{Float64, 2}(undef, size(T, 1), 2)
	sslack_ti = Array{Float64, 2}(undef, size(T, 1), 2)
	for t = 1:size(T, 1), i=1:2
		p_ti[t,i] = trafodata[4+(i-1)*3][Tidx[t]]
		q_ti[t,i] = trafodata[5+(i-1)*3][Tidx[t]]
		sslack_ti[t,i] = trafodata[6+(i-1)*3][Tidx[t]]
	end
	
	return p_li, q_li, p_ti, q_ti,
		pslackm_n, pslackp_n, qslackm_n, qslackp_n,
		sslack_li, sslack_ti
	
end

end
