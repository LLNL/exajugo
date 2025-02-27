#__precompile__()

module SolutionEvaluator

## elements to be exported

export BaseCaseSolution, BaseSolutionValue, BaseCaseValue,
	InitialBaseCaseSolution, InitialBaseCaseValue, InitialContingenciesSolutionsFromBase

## load external modules

using DataFrames, JuMP, GLPK
const MOI = JuMP.MOI
const MOIU = JuMP.MOIU

## load internal modules and functions
Base.include(@__MODULE__,"GoUtils.jl")
Base.include(@__MODULE__,"InstanceReader.jl")
using .GoUtils
import .InstanceReader: DELTA

## constants

## function to compute a full base-case solution (post-computing flows and slacks)

function powerflow(vi::Real, vj::Real, deltaij::Real, A::Real, B::Real, C::Real, Theta::Real=0)::Real
	return A*vi^2 + B*vi*vj*cos(deltaij + Theta) + C*vi*vj*sin(deltaij + Theta)
end

function BaseCaseSolution(v_n::V, theta_n::V, b_s::V, p_g::W, q_g::W,
	N::DataFrame, L::AbstractDataFrame, T::AbstractDataFrame, SSh::DataFrame, G::AbstractDataFrame,
	K::Union{Tuple{Symbol, Int}, Nothing}=nothing; IndexSets::Union{Tuple, Nothing}=nothing) where
	{V <: Vector{Float64}, W <: AbstractVector{Float64}}
	
	# compute index and sets for performant formulation
	if K == nothing
		K = DataFrame(Any[Int[], Symbol[], Int[]], [:Contingency,:ConType,:IDout])
	end
        L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn, K_outidx =
            build_or_split_indexsets(IndexSets, N, L, T, SSh, G, K)

	RefBus = G_Nidx[argmax(G[!, :Pub])]	# bus with the biggest generator is reference
	theta_n .-= theta_n[RefBus]
        isKaTuple = ( typeof(K) == Tuple{Symbol, Int} )
        if isKaTuple; @assert 1==size(K_outidx,1) "size of K_outidx not as expected"; end

	# compute line flows and slacks
	p_li = Array{Float64}(undef, size(L, 1), 2)
	q_li = Array{Float64}(undef, size(L, 1), 2)
	sslack_li = Array{Float64}(undef, size(L, 1), 2)
	for l=1:size(L, 1), i=1:2
		if isKaTuple && K[1] == :Line && K_outidx == l
			p_li[l,i] = 0.0
			q_li[l,i] = 0.0
			sslack_li[l,i] = 0.0
			continue
		end
		vi = v_n[L_Nidx[l,i]]
		vj = v_n[L_Nidx[l,3-i]]
		deltaij = theta_n[L_Nidx[l,i]] - theta_n[L_Nidx[l,3-i]]
		p_li[l,i] = powerflow(vi, vj, deltaij, L[!, :G][l], -L[!, :G][l], -L[!, :B][l])
		q_li[l,i] = powerflow(vi, vj, deltaij, -L[!, :B][l]-L[!, :Bch][l]/2, L[!, :B][l], -L[!, :G][l])
		sslack_li[l,i] = max(0.0, sqrt(p_li[l,i]^2 + q_li[l,i]^2) - L[!, :RateBase][l]*vi)
	end

	# compute transformer flows and slacks
	p_ti = Array{Float64}(undef, size(T, 1), 2)
	q_ti = Array{Float64}(undef, size(T, 1), 2)
	sslack_ti = Array{Float64}(undef, size(T, 1), 2)
	for t=1:size(T, 1)
                if isKaTuple &&  K[1] == :Transformer && K_outidx == t
			for i = 1:2
				p_ti[t,i] = 0.0
				q_ti[t,i] = 0.0
				sslack_ti[t,i] = 0.0
			end
			continue
                end
		v1 = v_n[T_Nidx[t,1]]
		v2 = v_n[T_Nidx[t,2]]
		delta12 = theta_n[T_Nidx[t,1]] - theta_n[T_Nidx[t,2]]
		p_ti[t,1] = powerflow(v1, v2, delta12, T[!, :G][t]/T[!, :Tau][t]^2+T[!, :Gm][t],
			-T[!, :G][t]/T[!, :Tau][t], -T[!, :B][t]/T[!, :Tau][t], -T[!, :Theta][t])
		q_ti[t,1] = powerflow(v1, v2, delta12, -T[!, :B][t]/T[!, :Tau][t]^2-T[!, :Bm][t],
			T[!, :B][t]/T[!, :Tau][t], -T[!, :G][t]/T[!, :Tau][t], -T[!, :Theta][t])
		p_ti[t,2] = powerflow(v2, v1, -delta12, T[!, :G][t],
			-T[!, :G][t]/T[!, :Tau][t], -T[!, :B][t]/T[!, :Tau][t], T[!, :Theta][t])
		q_ti[t,2] = powerflow(v2, v1, -delta12, -T[!, :B][t],
			T[!, :B][t]/T[!, :Tau][t], -T[!, :G][t]/T[!, :Tau][t], T[!, :Theta][t])
		for i=1:2
			sslack_ti[t,i] = max(0.0, sqrt(p_ti[t,i]^2 + q_ti[t,i]^2) - T[!, :RateBase][t])
		end
	end
	
	# compute power balance slacks
	pslackm_n = Vector{Float64}(undef, size(N, 1))
	pslackp_n = Vector{Float64}(undef, size(N, 1))
	qslackm_n = Vector{Float64}(undef, size(N, 1))
	qslackp_n = Vector{Float64}(undef, size(N, 1))
	for n=1:size(N, 1)
		pslack = 0.0
		qslack = 0.0
		if length(Gn[n]) > 0
			pslack += sum(p_g[g] for g=Gn[n])
			qslack += sum(q_g[g] for g=Gn[n])
		end
		pslack -= N[!, :Pd][n] + N[!, :Gsh][n]*v_n[n]^2
		qslack -= N[!, :Qd][n] - N[!, :Bsh][n]*v_n[n]^2
		if length(SShn[n]) > 0
			qslack -= sum(-b_s[s] for s=SShn[n])*v_n[n]^2
		end
		if length(Lidxn[n]) > 0
			pslack -= sum(p_li[Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n]))
			qslack -= sum(q_li[Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n]))
		end
		if length(Tidxn[n]) > 0
			pslack -= sum(p_ti[Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n]))
			qslack -= sum(q_ti[Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n]))
		end
		pslackm_n[n] = max(0.0, -pslack)
		pslackp_n[n] = max(0.0, pslack)
		qslackm_n[n] = max(0.0, -qslack)
		qslackp_n[n] = max(0.0, qslack)
	end
	
	## return base case solution
	return p_li, q_li, p_ti, q_ti, pslackm_n, pslackp_n, qslackm_n, qslackp_n,
		sslack_li, sslack_ti
	
end

## function to compute solution value given production and slacks

function productioncost(p::Real, p_h::T, c_h::T, LPSolver)::Real where {T <: Vector{V} where V <: Real}
	if length(p_h) != length(c_h)
		DimensionMismatch()
	end
	npairs = length(p_h)
	m = Model(LPSolver)
	@variable(m, t_h[1:npairs] >= 0)
	@constraint(m, sum(t_h[h] for h=1:npairs) == 1)
	@constraint(m, sum(p_h[h]*t_h[h] for h=1:npairs) == p)
	@objective(m, Min, sum(c_h[h]*t_h[h] for h=1:npairs))
	JuMP.optimize!(m)
	if JuMP.termination_status(m) != MOI.OPTIMAL
		error("error computing production cost, termination status: ", JuMP.termination_status(m))
	end
	return JuMP.objective_value(m)
end

function penaltycost(slack::Real, quantities_h::T, penalties_h::T, LPSolver)::Real where {T <: Vector{V} where V <: Real}
	if length(quantities_h) != length(penalties_h)
		DimensionMismatch()
	end
	ntranches = length(quantities_h)
	m = Model(LPSolver)
	@variable(m, 0 <= sigma_h[h=1:ntranches] <= quantities_h[h])
	@constraint(m, sum(sigma_h[h] for h=1:ntranches) == slack)
	@objective(m, Min, sum(penalties_h[h]*sigma_h[h] for h=1:ntranches))
	JuMP.optimize!(m)
	if JuMP.termination_status(m) != MOI.OPTIMAL
		error("error computing penalty cost, termination status: ", JuMP.termination_status(m))
	end
	return JuMP.objective_value(m)
end

function BaseSolutionValue(p_g::W, pslackm_n::T, pslackp_n::T, qslackm_n::T, qslackp_n::T,
	sslack_li::V, sslack_ti::V, G::AbstractDataFrame, P::DataFrame,
	LPSolver=with_optimizer(GLPK.Optimizer, msg_lev=GLPK.MSG_OFF);
	ReturnDetailedComposition::Bool=false) where {W <: AbstractVector{Float64},
	T <: Vector{Float64}, V <: Array{Float64, 2}}
	
	# compute production cost
	c_g = Vector{Float64}(undef, length(p_g))
	for g = 1:length(p_g)
		c_g[g] = productioncost(p_g[g], G[!, :CostPi][g], G[!, :CostCi][g], LPSolver)
	end
	
	# compute penalties on power balance
	csigmap_n = Vector{Float64}(undef, length(pslackm_n))
	csigmaq_n = Vector{Float64}(undef, length(pslackm_n))
	pidx = findfirst(P[!, :Slack] .== :P)
	qidx = findfirst(P[!, :Slack] .== :Q)
	for n = 1:length(pslackm_n)
		csigmap_n[n] = penaltycost(pslackm_n[n] + pslackp_n[n],
			P[!, :Quantities][pidx], P[!, :Penalties][pidx], LPSolver)
		csigmaq_n[n] = penaltycost(qslackm_n[n] + qslackp_n[n],
			P[!, :Quantities][qidx], P[!, :Penalties][qidx], LPSolver)
	end
	
	# compute penalties on thermal limits
	sidx = findfirst(P[!, :Slack] .== :S)
	csigmas_li = Array{Float64}(undef, size(sslack_li))
	for l=1:size(sslack_li,1), i=1:2
		csigmas_li[l,i] = penaltycost(sslack_li[l,i],
			P[!, :Quantities][sidx], P[!, :Penalties][sidx], LPSolver)
	end
	csigmas_ti = Array{Float64}(undef, size(sslack_ti))
	for t=1:size(sslack_ti,1), i=1:2
		csigmas_ti[t,i] = penaltycost(sslack_ti[t,i],
			P[!, :Quantities][sidx], P[!, :Penalties][sidx], LPSolver)
	end
	
	# compute objective
	objproduction = sum(c_g[g] for g=1:length(p_g))
	objpenalties = DELTA*(sum((csigmap_n[n] + csigmaq_n[n]) for n=1:length(pslackm_n)) +
		sum(csigmas_li[l,i] for l=1:size(sslack_li, 1), i=1:2) +
		sum(csigmas_ti[t,i] for t=1:size(sslack_ti, 1), i=1:2))
	
	# return objective
	if ReturnDetailedComposition
		return objproduction+objpenalties, c_g, csigmap_n, csigmaq_n, csigmas_li, csigmas_ti
	else
		return objproduction+objpenalties, objproduction, objpenalties
	end
	
end

## function to compute base-case value of a base-case solution

function BaseCaseValue(v_n::V, theta_n::V, b_s::V, p_g::W, q_g::W,
	N::DataFrame, L::AbstractDataFrame, T::AbstractDataFrame, SSh::DataFrame, G::AbstractDataFrame,
	P::DataFrame, LPSolver=GLPK.Optimizer;
	# P::DataFrame, LPSolver=with_optimizer(GLPK.Optimizer, msg_lev=GLPK.MSG_OFF);
	ReturnDetailedComposition::Bool=false, IndexSets::Union{Tuple, Nothing}=nothing) where
	{V <: Vector{Float64}, W <: AbstractVector{Float64}}
	
	# compute slacks
	_, _, _, _, pslackm_n, pslackp_n, qslackm_n, qslackp_n,
		sslack_li, sslack_ti = BaseCaseSolution(v_n, theta_n, b_s, p_g, q_g,
		N, L, T, SSh, G, IndexSets=IndexSets)
	
	# return solution value
	return BaseSolutionValue(p_g, pslackm_n, pslackp_n, qslackm_n, qslackp_n,
		sslack_li, sslack_ti, G, P, LPSolver,
		ReturnDetailedComposition=ReturnDetailedComposition)
	
end

## function to compute full initial base-case solution

function InitialBaseCaseSolution(N::DataFrame, L::DataFrame, T::DataFrame, SSh::DataFrame,
	G::DataFrame; IndexSets::Union{Tuple, Nothing}=nothing)
	
	# collect initial solution provided by the competition
	v0_n = N[!, :v0]
	theta0_n = N[!, :theta0]
	b0_s = SSh[!, :B0]
	p0_g = G[!, :p0]
	q0_g = G[!, :q0]
	
	# call base solution computation function
	p0_li, q0_li, p0_ti, q0_ti, pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n,
		sslack0_li, sslack0_ti = BaseCaseSolution(v0_n, theta0_n, b0_s, p0_g, q0_g,
		N, L, T, SSh, G, IndexSets=IndexSets)
	
	## return initial base-case solution
	return v0_n, theta0_n, p0_li, q0_li, p0_ti, q0_ti, b0_s, p0_g, q0_g, pslack0m_n, pslack0p_n,
		qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti
	
end

## function to compute initial base-case solution value

function InitialBaseCaseValue(N::DataFrame, L::DataFrame, T::DataFrame, SSh::DataFrame,
	G::DataFrame, P::DataFrame, LPSolver=with_optimizer(GLPK.Optimizer, msg_lev=GLPK.MSG_OFF);
	ReturnDetailedComposition::Bool=false, IndexSets::Union{Tuple, Nothing}=nothing)
	
	# collect initial solution provided by the competition
	v0_n = N[!, :v0]
	theta0_n = N[!, :theta0]
	b0_s = SSh[!, :B0]
	p0_g = G[!, :p0]
	q0_g = G[!, :q0]
	
	## return inital base-case solution value
	return BaseCaseValue(v0_n, theta0_n, b0_s, p0_g, q0_g, N, L, T, SSh, G, P,
		LPSolver, ReturnDetailedComposition=ReturnDetailedComposition,
		IndexSets=IndexSets)
	
end

function 
    InitialContingenciesSolutionsFromBase(N, L, T, SSh, G, K,
                                          v0_n, theta0_n, p0_li, q0_li, p0_ti, q0_ti, b0_s, p0_g, q0_g, 
                                          pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti;
                                          IndexSets)
    nK = size(K, 1)

    v0_nk       = repeat(v0_n, 1, nK)
    theta0_nk   = repeat(theta0_n, 1, nK)
    p0_lik      = repeat(p0_li, inner=(1, 1, nK))
    q0_lik      = repeat(q0_li, inner=(1, 1, nK))
    p0_tik      = repeat(p0_ti, inner=(1, 1, nK))
    q0_tik      = repeat(q0_ti, inner=(1, 1, nK))
    b0_sk       = repeat(b0_s, 1, nK)
    p0_gk       = repeat(p0_g, 1, nK)
    q0_gk       = repeat(q0_g, 1, nK)
    pslack0m_nk = repeat(pslack0m_n, 1, nK)
    pslack0p_nk = repeat(pslack0p_n, 1, nK)
    qslack0m_nk = repeat(qslack0m_n, 1, nK)
    qslack0p_nk = repeat(qslack0p_n, 1, nK)
    sslack0_lik = repeat(sslack0_li, inner=( 1, 1, nK))
    sslack0_tik = repeat(sslack0_ti, inner=(1, 1, nK))

    return v0_nk, theta0_nk, p0_lik, q0_lik, p0_tik, q0_tik, b0_sk, p0_gk, q0_gk, pslack0m_nk, pslack0p_nk, qslack0m_nk, qslack0p_nk, sslack0_lik, sslack0_tik
end

end
