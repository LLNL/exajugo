#__precompile__()

module MPACOPF

## elements to be exported

export SolveMPACOPF, SolveMaster, SolveMasterSQP, SolveMPACOPFsingle

## load external modules

using DataFrames, Roots, JuMP, Printf, SparseArrays, LinearAlgebra, MathOptInterface
const MOI = JuMP.MOI
const MOIU = JuMP.MOIU
import JuMP: add_to_expression!

## load internal modules and functions
# Base.include(@__MODULE__,"InstanceReader.jl")
Base.include(@__MODULE__,"InstanceReaderMP.jl")
Base.include(@__MODULE__,"SmoothApproximations.jl")
Base.include(@__MODULE__,"GoUtils.jl")
Base.include(@__MODULE__,"SolutionEvaluator.jl")
using  .GoUtils, .SmoothApproximations, .InstanceReaderMP
import .SolutionEvaluator: InitialBaseCaseSolution, BaseCaseSolution, InitialContingenciesSolutionsFromBase

import DataStructures: OrderedDict

## GoOptions struct
include("GoOptions.jl")
include("MPACOPFStates.jl")
include("SCACOPFStates.jl")
## Reserves
include("SCACOPFReserves.jl")

## constants

const MAXANGLEDIFF = pi/2	# radians
const EPSILON = 1E-12		# -
scalSubObj = 1

## function to get the start value of a variable or fixed variable (fixed before JuMP)

function getstartvalue(x::Union{JuMP.VariableRef, T})::Real where {T <: Real}
	if typeof(x) <: JuMP.VariableRef
		value = JuMP.start_value(x)
		if value == nothing
			return T(0)
		else
			return value
		end
	else
		return x
	end
end

## function to add a rectangular power flow constraint to a model

function addpowerflowcon!(m::JuMP.Model, pq::Union{JuMP.VariableRef, Real},
	vi::Union{JuMP.VariableRef, Real}, vj::Union{JuMP.VariableRef, Real},
	thetai::Union{JuMP.VariableRef, Real}, thetaj::Union{JuMP.VariableRef, Real},
	A::Real, B::Real, C::Real, Theta::Real=0)
	return @NLconstraint(m, pq == A*vi^2 + B*vi*vj*cos(thetai - thetaj + Theta) +
		C*vi*vj*sin(thetai - thetaj + Theta))
end

## function to add all power flow constraints to a model

function addpowerflowcons!(m::JuMP.Model, v_n::AbstractVector, theta_n::AbstractVector,
	p_li::AbstractArray, q_li::AbstractArray, p_ti::AbstractArray, q_ti::AbstractArray,
	b_s::AbstractVector, p_g::AbstractVector, q_g::AbstractVector,
	pslackm_n::AbstractVector, pslackp_n::AbstractVector,
	qslackm_n::AbstractVector, qslackp_n::AbstractVector,
	sslack_li::AbstractArray, sslack_ti::AbstractArray,
	N::DataFrame, L::DataFrame, T::DataFrame, 
	# IndexSets::Tuple, ThermalLimits::Symbol;
	IndexSets::Tuple, ThermalLimits::Symbol, ts::Int64, Linfo::Vector{Any}, Tinfo::Vector{Any};
	SysCond::Symbol=:BaseCase, ConType::Symbol=:None, outidx::Int=0)

	# parse necessary index sets
	L_Nidx = IndexSets[1]
	T_Nidx = IndexSets[2]
	Lidxn = IndexSets[5]
	Lin = IndexSets[6]
	Tidxn = IndexSets[7]
	Tin = IndexSets[8]
	SShn = IndexSets[9]
	Gn = IndexSets[10]
	
	# thermal rating to be used
	if SysCond == :BaseCase
		RateSymb = :RateBase
	else
		RateSymb = :RateEmer
	end
	
	# thermal limits and power flows -- lines
	for l=1:size(L, 1), i=1:2
		if SysCond == :Contingency && ConType == :Line && outidx == l
			continue
		end
		if ThermalLimits == :quadr
			@constraint(m, p_li[l,i]^2 + q_li[l,i]^2 <=
				(L[!, RateSymb][l]*v_n[L_Nidx[l,i]] + sslack_li[l,i])^2)
                push!(Linfo,(ts=ts,l=l,i=i))
		else
			@NLconstraint(m, (p_li[l,i]^2 + q_li[l,i]^2+1e-6)^0.5 <=
				L[!, RateSymb][l]*v_n[L_Nidx[l,i]]+ sslack_li[l,i])
		end

		addpowerflowcon!(m, p_li[l,i], v_n[L_Nidx[l,i]], v_n[L_Nidx[l,3-i]],
			theta_n[L_Nidx[l,i]], theta_n[L_Nidx[l,3-i]],
			L[!,:G][l], -L[!,:G][l], -L[!,:B][l])
		addpowerflowcon!(m, q_li[l,i], v_n[L_Nidx[l,i]], v_n[L_Nidx[l,3-i]],
			theta_n[L_Nidx[l,i]], theta_n[L_Nidx[l,3-i]],
			-L[!,:B][l]-L[!,:Bch][l]/2, L[!,:B][l], -L[!,:G][l])
	end
	
    	# thermal limits and power flows -- transformers
	for t=1:size(T, 1), i=1:2
		if SysCond == :Contingency && ConType == :Transformer && outidx == t
			println("transformer continue time ", ts, "trans", t, "i ", i)
			continue
		end
		if ThermalLimits == :quadr
			@constraint(m, p_ti[t,i]^2 + q_ti[t,i]^2 <=
				(T[!, RateSymb][t]+ sslack_ti[t,i])^2 )
                push!(Tinfo,(ts=ts,t=t,i=i))
		else
			@NLconstraint(m, (p_ti[t,i]^2 + q_ti[t,i]^2+1e-6)^0.5 <=
				T[!, RateSymb][t]+ sslack_ti[t,i])
		end

        if i == 1
            addpowerflowcon!(m, p_ti[t,1], v_n[T_Nidx[t,1]], v_n[T_Nidx[t,2]],
                theta_n[T_Nidx[t,1]], theta_n[T_Nidx[t,2]],
			    T[!,:G][t]/T[!,:Tau][t]^2+T[!,:Gm][t], -T[!,:G][t]/T[!,:Tau][t],
			    -T[!,:B][t]/T[!,:Tau][t], -T[!,:Theta][t])
		    addpowerflowcon!(m, q_ti[t,1], v_n[T_Nidx[t,1]], v_n[T_Nidx[t,2]],
			    theta_n[T_Nidx[t,1]], theta_n[T_Nidx[t,2]],
			    -T[!,:B][t]/T[!,:Tau][t]^2-T[!,:Bm][t], T[!,:B][t]/T[!,:Tau][t],
			    -T[!,:G][t]/T[!,:Tau][t], -T[!,:Theta][t])
        else
            addpowerflowcon!(m, p_ti[t,2], v_n[T_Nidx[t,2]], v_n[T_Nidx[t,1]],
			    theta_n[T_Nidx[t,2]], theta_n[T_Nidx[t,1]],
			    T[!,:G][t], -T[!,:G][t]/T[!,:Tau][t],
			    -T[!,:B][t]/T[!,:Tau][t], T[!,:Theta][t])
		    addpowerflowcon!(m, q_ti[t,2], v_n[T_Nidx[t,2]], v_n[T_Nidx[t,1]],
			    theta_n[T_Nidx[t,2]], theta_n[T_Nidx[t,1]],
			    -T[!,:B][t], T[!,:B][t]/T[!,:Tau][t],
			    -T[!,:G][t]/T[!,:Tau][t], T[!,:Theta][t])
        end
	end

	for n = 1:size(N, 1)
		@NLconstraint(m, sum(q_g[g] for g=Gn[n]) - N[!,:Qd][n] 
			- (-N[!,:Bsh][n] - sum(b_s[s] for s=SShn[n]))*v_n[n]^2
			- sum(q_li[Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n])) 
			- sum(q_ti[Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n])) 
			== qslackp_n[n] - qslackm_n[n])
	end

end

## function to add piecewise linear penalty function model

function startvalueauxvarspenaltyfuncs(slackstart::Real, quantities::Vector{Float64},
	penalties::Vector{Float64})::Vector{Float64}
	@assert slackstart >= 0 "start value for slack should be non-negative"
	idx = sortperm(penalties)
	startvalue = zeros(Float64, length(quantities))
	for i = idx
		startvalue[i] = min(quantities[i], slackstart)
		slackstart -= startvalue[i]
		if slackstart < 1e-12
			break
		end
	end
	return startvalue
end

function addpenaltyfunction!(m::JuMP.Model, slack::Union{JuMP.GenericAffExpr, JuMP.VariableRef},
	quantities::Vector{Float64}, penalties::Vector{Float64}, SetStart::Bool=true,
	slackstart::Union{Real, Nothing}=nothing)
	if length(quantities) != length(penalties)
		DimensionMismatch("quantities and penalties should have the same length")
	end
	ntranches = length(quantities)
	sigma_h = @variable(m, [h=1:ntranches], lower_bound=0.0, upper_bound=quantities[h])
	if SetStart
		sigma0_h = startvalueauxvarspenaltyfuncs(slackstart, quantities, penalties)
		JuMP.set_start_value.(sigma_h, sigma0_h)
	end
	@constraint(m, slack == sum(sigma_h[h] for h=1:ntranches))
	return @expression(m, sum(penalties[h]*sigma_h[h] for h=1:ntranches))
end

## function to add piecewise linear penalty for a case 

function addpenaltyfunctions!(m::JuMP.Model,
	pslackm_n::AbstractVector, pslackp_n::AbstractVector,
	qslackm_n::AbstractVector, qslackp_n::AbstractVector,
	sslack_li::AbstractArray, sslack_ti::AbstractArray,
	N::DataFrame, L::DataFrame, T::DataFrame, P::DataFrame,
	SetStart::Bool=true)
	
	# define expression to collect all penalty terms
	casepenalty = zero(JuMP.AffExpr)
	
	# shorthand for penalty functions
	penaltyfun(x, idx, startvalue) = addpenaltyfunction!(m, x, P[!,:Quantities][idx],
		P[!,:Penalties][idx], SetStart, startvalue)
	penaltyfun(x, idx) = addpenaltyfunction!(m, x, P[!,:Quantities][idx], P[!,:Penalties][idx],
		SetStart, if SetStart getstartvalue(x) else nothing end)

	# imbalance penalty
	idxP = findfirst(P[!,:Slack] .== :P)
	idxQ = findfirst(P[!,:Slack] .== :Q)
	if SetStart
		for n = 1:size(N, 1)
			add_to_expression!(casepenalty, penaltyfun(pslackm_n[n] + pslackp_n[n], idxP,
				getstartvalue(pslackm_n[n]) + getstartvalue(pslackp_n[n])))
			add_to_expression!(casepenalty, penaltyfun(qslackm_n[n] + qslackp_n[n], idxQ,
				getstartvalue(qslackm_n[n]) + getstartvalue(qslackp_n[n])))
		end
	else
		for n = 1:size(N, 1)
			add_to_expression!(casepenalty, penaltyfun(pslackm_n[n] + pslackp_n[n], idxP))
			add_to_expression!(casepenalty, penaltyfun(qslackm_n[n] + qslackp_n[n], idxQ))
		end
	end
	
	# thermal overload
	idxS = findfirst(P[!,:Slack] .== :S)
	for l = 1:size(L, 1), i=1:2
		add_to_expression!(casepenalty, penaltyfun(sslack_li[l,i], idxS))
	end
	for t = 1:size(T, 1), i=1:2
		add_to_expression!(casepenalty, penaltyfun(sslack_ti[t,i], idxS))
	end
	# return expression containing all penalty terms
	return casepenalty
	
end

function getstartcostcoefficients(p_g0::Float64, costsamplepoints::Vector{Float64}, lbg::Float64, ubg::Float64)
    #
    # #! some consistency checks 
    #
    if findmin(costsamplepoints)[1] > lbg + 1e-8
        @warn "lower bounds inconsistent with generator cost sample points!?!"
        println("p_g0=", p_g0, " ", 
                "Pi pts:", costsamplepoints, 
                "  lower bound:", lbg, 
                " and min cost point > generator lower bound (violation is ", lbg-findmin(costsamplepoints)[1], ")")
    end
    
    if findmax(costsamplepoints)[1] < ubg - 1e-8
        @warn "upper bounds inconsistent with generator cost sample points!?!"
        println("p_g0=", p_g0, " ", 
                "Pi pts:", costsamplepoints, 
                "  upper bound:", ubg, 
                " and min cost point > generator lower bound (violation is ", -ubg + findmax(costsamplepoints)[1], ")")
            end
    # #! end of checks

    costcoeff=zeros(length(costsamplepoints))
    i = findfirst(costsamplepoints .> p_g0)
    if i==1
        mincostpoint = findmin(costsamplepoints)[1]
        if mincostpoint-p_g0 > 1e-8
            println("mincostpoint=", mincostpoint, " but p_g0=", p_g0, " violation=",  mincostpoint-p_g0)
            @warn false "p_g0 is too much outside the min of generation cost points"
            p_g0 = mincostpoint
        else
            p_g0 = mincostpoint
        end
        i = findfirst(costsamplepoints .> p_g0)
    end
    if i==nothing
        #p_g0 >= max of costsample point
        (maxcostpoint, imax) = findmax(costsamplepoints)
        if p_g0>maxcostpoint+1e-8
            println("maxcostpoint=", maxcostpoint, " but p_g0=", p_g0, " violation=",  p_g0-maxcostpoint)
            @warn false "p_g0 is too much outside the max of generation cost points"
            p_g0 = maxcostpoint
        else
            p_g0 = maxcostpoint
        end
        i = imax
    end
    @assert i!= nothing
    @assert i>=2
    @assert i<= length(costsamplepoints)

    @assert costsamplepoints[i] >= p_g0
    @assert costsamplepoints[i-1] <= p_g0
    @assert costsamplepoints[i]-costsamplepoints[i-1] >= 1e-10

    costcoeff[i-1] = 1 - (p_g0-costsamplepoints[i-1])/(costsamplepoints[i]-costsamplepoints[i-1])
    costcoeff[i]   = 1 - costcoeff[i-1]
    return costcoeff
end

## function to construct quadratic penalty terms for a case 

function addquadpenaltyblock!(expr::JuMP.GenericQuadExpr,
	slack::AbstractArray, a::Float64, b::Float64)::Nothing
	# a*sum(x^2) + b*sum(x)
	for i in eachindex(slack)
		add_to_expression!(expr, a, slack[i], slack[i])
	end
	for i in eachindex(slack)
		add_to_expression!(expr, b, slack[i])
	end
	return nothing
end

function quadpenaltyfunctions(
	pslackm_n::AbstractVector, pslackp_n::AbstractVector,
	qslackm_n::AbstractVector, qslackp_n::AbstractVector,
	sslack_li::AbstractArray, sslack_ti::AbstractArray,
	P::DataFrame)
	# plantslack::AbstractVector, P::DataFrame)
	
	# define expression to collect all penalty terms
	casepenalty = zero(JuMP.QuadExpr)
	
	# compute coefficients (a*x^2 + b*x)
	function quadcoeffs(idx)
		perm = sortperm(P[!,:Penalties][idx])
		quantities = view(P[!,:Quantities][idx], perm)
		penalties = view(P[!,:Penalties][idx], perm)
		b = penalties[1]
		a = (penalties[2] - b)/(2*quantities[1])
		@assert b > 0 && a > 0
		return a*100, b*100
	end
	aP, bP = quadcoeffs(findfirst(P[!,:Slack] .== :P))
	aQ, bQ = quadcoeffs(findfirst(P[!,:Slack] .== :Q))
	aS, bS = quadcoeffs(findfirst(P[!,:Slack] .== :S))
	
	# println("aP: ", aP, "bP: ", bP)
	# println("aQ: ", aQ, "bQ: ", bQ)
	# println("aS: ", aS, "bS: ", bS)
	
	# imbalance penalties
	addquadpenaltyblock!(casepenalty, pslackm_n, aP, bP)
	addquadpenaltyblock!(casepenalty, pslackp_n, aP, bP)
	addquadpenaltyblock!(casepenalty, qslackm_n, aQ, bQ)
	addquadpenaltyblock!(casepenalty, qslackp_n, aQ, bQ)
	
	# thermal overload
	addquadpenaltyblock!(casepenalty, sslack_li, aS, bS)
	addquadpenaltyblock!(casepenalty, sslack_ti, aS, bS)

	# # zero out plant power on other buses
	# addquadpenaltyblock!(casepenalty, plantslack, aP, bP)

	# return expression containing all penalty terms
	return casepenalty
	
end

function Ylinkpen(YUslack::AbstractVector, YZslack::AbstractVector, YVslack::AbstractVector,
	 YUmslack::AbstractVector, YVmslack::AbstractVector)

	 # define expression to collect all penalty terms
	 Ylinkpenalty = zero(JuMP.QuadExpr)

	 # compute coefficients (a*x^2 + b*x)
	 a = 1e8
	#  b = 1e8
	#  a = 0.0
	 b = 0.0

	 # infeasibility penalties
	addquadpenaltyblock!(Ylinkpenalty, YUslack, a, b)
	addquadpenaltyblock!(Ylinkpenalty, YZslack, a, b)
	addquadpenaltyblock!(Ylinkpenalty, YVslack, a, b)
	addquadpenaltyblock!(Ylinkpenalty, YUmslack, a, b)
	addquadpenaltyblock!(Ylinkpenalty, YVmslack, a, b)

	# return expression containing all penalty terms
	return Ylinkpenalty

end

## function for solving each subproblem
function SolveSub(o::GoOptions, N::DataFrame, L::DataFrame, T::DataFrame,
	SSh::DataFrame, G::DataFrame, K::DataFrame, P::DataFrame, plant::DataFrame, numP::Int64, MVAbase::Float64, Ylink::Ylinkdecomts, timestep::Vector{Int64}, NLSolver; 
		IndexSets::Union{Tuple, Nothing}=nothing,
		StartingPoint::Union{SCACOPFState,Nothing}=nothing)

	if numP > 3
		println("U[1:2],V[1:3],Z[1:2] ", Ylink.YUk[1:2], " ", Ylink.YVk[1:3], " ", Ylink.YZk[1:2])
	else
		println("U, V[1:3], Z ", Ylink.YUk, " ", Ylink.YVk[1:3], " ", Ylink.YZk)
	end
	flush(stdout)
	# compute index and sets for performant formulation
	L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn, K_outidx =
	build_or_split_indexsets(IndexSets, N, L, T, SSh, G, K)

	RefBus = G_Nidx[argmax(G[!,:Pub])]	# bus with the biggest generator is reference

	# (some) output for options used
	println("SolveSCACOPF with Options: CouplingMode=", o.CouplingMode, 
  	" Smoothing=", o.SmoothingParamCoupling,
	" ThermalLimits=", o.ThermalLimits)
	

	if StartingPoint==nothing
		# compute SC AC OPF starting point 
		v0_n, theta0_n, p0_li, q0_li, p0_ti, q0_ti, b0_s, p0_g, q0_g, 
		pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti =
 	 InitialBaseCaseSolution( N, L, T, SSh, G, 
						   IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
								SShn, Gn, K_outidx))
		theta0_n .-= theta0_n[RefBus]
	else
		## base case
		v0_n, theta0_n, b0_s, p0_g, q0_g = getBaseStates(StartingPoint)
		p0_li, q0_li, p0_ti, q0_ti, pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti = 
 	 	BaseCaseSolution(v0_n, theta0_n, b0_s, p0_g, q0_g,
				   N, L, T, SSh, G, 
				   IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
						SShn, Gn, K_outidx))
		println("HOT Started")
	end

	numG = size(G,1)
	# create model
	m = Model(NLSolver)

	# ACOPF variables
	@variable(m, N[!,:Vlb][n] <= v_n[n=1:size(N, 1)] <= N[!,:Vub][n], start = v0_n[n])
	@variable(m, theta_n[n=1:size(N, 1)], start = theta0_n[n])
	@variable(m, p_li[l=1:size(L, 1), i=1:2], start = p0_li[l,i])
	@variable(m, q_li[l=1:size(L, 1), i=1:2], start = q0_li[l,i])
	@variable(m, p_ti[t=1:size(T, 1), i=1:2], start = p0_ti[t,i])
	@variable(m, q_ti[t=1:size(T, 1), i=1:2], start = q0_ti[t,i])
	@variable(m, SSh[!,:Blb][s] <= b_s[s=1:size(SSh, 1)] <= SSh[!,:Bub][s], start = b0_s[s])
	@variable(m, G[!,:Plb][g] <= p_g[g=1:size(G, 1)] <= G[!,:Pub][g], start = p0_g[g])
	@variable(m, G[!,:Qlb][g] <= q_g[g=1:size(G, 1)] <= G[!,:Qub][g], start = q0_g[g])
	@variable(m, pslackm_n[n=1:size(N, 1)] >= 0, start = pslack0m_n[n])
	@variable(m, pslackp_n[n=1:size(N, 1)] >= 0, start = pslack0p_n[n])
	@variable(m, qslackm_n[n=1:size(N, 1)] >= 0, start = qslack0m_n[n])
	@variable(m, qslackp_n[n=1:size(N, 1)] >= 0, start = qslack0p_n[n])
	@variable(m, sslack_li[l=1:size(L, 1), i=1:2] >= 0, start = sslack0_li[l,i])
	@variable(m, sslack_ti[t=1:size(T, 1), i=1:2] >= 0, start = sslack0_ti[t,i])

	# slack variables for Ylink 
	@variable(m, YUslack[vs=1:numP*2] >= 0, start = 0.0)
	@variable(m, YZslack[vs=1:numP*2] >= 0, start = 0.0)
	@variable(m, YVslack[vs=1:numG*2] >= 0, start = 0.0)
	@variable(m, YUmslack[vs=1:numP*2] >= 0, start = 0.0)
	@variable(m, YVmslack[vs=1:numG*2] >= 0, start = 0.0)

	# bidding entity variables
	@variable(m, plant[!,:pplantlb][np] <= pplant[np=1:numP] <= plant[!,:pplantub][np])
	@variable(m, plprod[np=1:numP] >= 0, start = 0.0)
	
	# bidding plant constraints
	# calculate plant prod from energy consumption at corresponding bus
	@constraint(m, [np=1:numP], plprod[np] == plant[!,:plantprod][np].*pplant[np])
	# ramping constraint for plant production 
	if Ylink.ts < timestep[end]
		# @constraint(m, CUU[np=1:numP], pplant[np] <= Ylink.YUk[np]+plant[!,:plantRRL][np]/2)
		# @constraint(m, CUD[np=1:numP], Ylink.YUk[np] - plant[!,:plantRRL][np]/2 <= pplant[np])
		@constraint(m, CUU[np=1:numP], pplant[np] <= Ylink.YUk[np]+plant[!,:plantRRL][np]/2 + YUslack[np])
		@constraint(m, CUD[np=1:numP], Ylink.YUk[np] - plant[!,:plantRRL][np]/2 - YUslack[np+numP] <= pplant[np])
	end
	if Ylink.ts > timestep[1]
		# @constraint(m, CUUE[np=1:numP], pplant[np] <= Ylink.YUkm[np]+plant[!,:plantRRL][np]/2)
		# @constraint(m, CUDE[np=1:numP], Ylink.YUkm[np]/2 - plant[!,:plantRRL][np] <= pplant[np])
		@constraint(m, CUUE[np=1:numP], pplant[np] <= Ylink.YUkm[np]+plant[!,:plantRRL][np]/2 + YUmslack[np])
		@constraint(m, CUDE[np=1:numP], Ylink.YUkm[np] - plant[!,:plantRRL][np]/2 - YUmslack[np+numP] <= pplant[np])
	end
	# temporal constraint: plant production within demand and flex*demand 
	Du = (plant[!,:plantflex].+1).*(plant[!,:plantdemand]./24)
	Dl = plant[!,:plantdemand]./24
	@constraint(m, CZU[np=1:numP], plprod[np] - Ylink.YZk[np] <= Du[np] + YZslack[np])
	@constraint(m, CZL[np=1:numP], Dl[np] - YZslack[np+numP] <= plprod[np] - Ylink.YZk[np])
	# @constraint(m, CZU[np=1:numP], plprod[np] - Ylink.YZk[np] <= (1+plant[!,:plantflex][1]).*plant[!,:plantdemand][1]/24)
	# @constraint(m, CZL[np=1:numP], plant[!,:plantdemand][1]/24 <= plprod[np] - Ylink.YZk[np])

	# ramping constraints for power generation
	if Ylink.ts < timestep[end]
		@constraint(m, CVU[g=1:numG], p_g[g] <= Ylink.YVk[g] + 1*G[!,:Pub][g]/2 + YVslack[g])
		@constraint(m, CVD[g=1:numG], Ylink.YVk[g] - 1*G[!,:Pub][g]/2 - YVslack[g+numG] <= p_g[g])
		# @constraint(m, CVU[g=1:numG], p_g[g] <= Ylink.YVk[g] + 1*G[!,:Pub][g]/2)
		# @constraint(m, CVD[g=1:numG], Ylink.YVk[g] - 1*G[!,:Pub][g]/2 <= p_g[g])
	end
	if Ylink.ts > timestep[1]
		@constraint(m, CVUE[g=1:numG], p_g[g] <= Ylink.YVkm[g] + 1*G[!,:Pub][g]/2 + YVmslack[g])
		@constraint(m, CVDE[g=1:numG], Ylink.YVkm[g] - 1*G[!,:Pub][g]/2 - YVmslack[g+numG] <= p_g[g])
		# @constraint(m, CVUE[g=1:numG], p_g[g] <= Ylink.YVkm[g] + 1*G[!,:Pub][g]/2)
		# @constraint(m, CVDE[g=1:numG], Ylink.YVkm[g] - 1*G[!,:Pub][g]/2 <= p_g[g])
	end

	#active power balance with LMP=dual
	@constraint(m, PBPC[np=1:numP], sum(p_g[g] for g=Gn[plant[!,:busidx][np]]) - N[!,:Pd][plant[!,:busidx][np]] - N[!,:Gsh][plant[!,:busidx][np]]*v_n[plant[!,:busidx][np]]^2 -
			pplant[np]- 
			sum(p_li[Lidxn[plant[!,:busidx][np]][lix],Lin[plant[!,:busidx][np]][lix]] for lix=1:length(Lidxn[plant[!,:busidx][np]])) -
			sum(p_ti[Tidxn[plant[!,:busidx][np]][tix],Tin[plant[!,:busidx][np]][tix]] for tix=1:length(Tidxn[plant[!,:busidx][np]])) ==
			pslackp_n[plant[!,:busidx][np]] - pslackm_n[plant[!,:busidx][np]])

	nbusn = findall( x -> !(x in plant[!,:plantid]), N[!, :Bus])
	@constraint(m, PB[n=nbusn], sum(p_g[g] for g=Gn[n]) - N[!,:Pd][n] - N[!,:Gsh][n]*v_n[n]^2 -
	sum(p_li[Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n])) -
	sum(p_ti[Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n])) ==
	pslackp_n[n] - pslackm_n[n])

	# fix angle at reference bus to zero
	JuMP.fix(theta_n[RefBus], 0.0, force=true)

	# add power flow constraints
	addpowerflowcons!(m, v_n, theta_n, p_li, q_li, p_ti, q_ti, b_s, p_g, q_g,
	pslackm_n, pslackp_n, qslackm_n, qslackp_n, sslack_li, sslack_ti,
	N, L, T, (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
	SShn, Gn, K_outidx), o.ThermalLimits)

	# register approximation functions
	sclamp = nothing
	sstep = nothing

	## objective function
	
	# production cost
	prodcost = zero(JuMP.AffExpr)
	# basepenalty = zero(JuMP.QuadExpr)
	plantrev = zero(JuMP.AffExpr)

	for g = 1:size(G, 1)
	    npairs = length(G[!,:CostPi][g])
            @assert sort(G[!,:CostPi][g]) == G[!,:CostPi][g]
 
            t_h0 = getstartcostcoefficients(JuMP.start_value(p_g[g]), G[!,:CostPi][g], G[!,:Plb][g], G[!,:Pub][g])
            t_h = @variable(m, [i=1:npairs], lower_bound=0.0, start=t_h0[i])
            
 	    @constraint(m, sum(t_h[h] for h=1:npairs) == 1.0)
	    @constraint(m, p_g[g] == sum(G[!,:CostPi][g][h]*t_h[h] for h=1:npairs))
            
	    add_to_expression!(prodcost, @expression(m, sum(G[!,:CostCi][g][h]*t_h[h] for h=1:npairs)))
	end
	
	# base case penalties
	basepenalty = quadpenaltyfunctions(pslackm_n, pslackp_n, qslackm_n, qslackp_n,
				sslack_li, sslack_ti, P)

	# plant prod revenue
	for np=1:numP
		add_to_expression!(plantrev, pplant[np].*plant[!,:plantrev][np])
	end

	# Ylink penalties
	Ylinkpennalty = Ylinkpen(YUslack, YZslack, YVslack, YUmslack, YVmslack)

	# declare objective
	# @objective(m, Min, (prodcost - plantrev) + DELTA*basepenalty + DELTA*Ylinkpennalty)
	@objective(m, Min, scalSubObj*( (prodcost - plantrev) + DELTA*basepenalty + DELTA*Ylinkpennalty) )
	
	## attempt to solve ACOPF
	JuMP.optimize!(m)

	    # termination_status, primal_status, and dual_status: http://www.juliaopt.org/JuMP.jl/dev/solutions/
        solvestatusprimal = JuMP.primal_status(m)
        solvestatus = JuMP.termination_status(m)
        println("JuMP optimize MOI termination_status: ", solvestatus, " ", solvestatusprimal)
	if (solvestatusprimal != MOI.FEASIBLE_POINT) && (solvestatusprimal != MOI.NEARLY_FEASIBLE_POINT)
		error("solver failed to find a feasible solution")
	end
	
	## objective breakdown
	prodcostvalue = scalSubObj*JuMP.value(prodcost)

	penaltybasevalue = scalSubObj*DELTA*JuMP.value(basepenalty)
  
	plantrevvalue = scalSubObj*JuMP.value(plantrev)
  
	penaltyYlinkvalue = scalSubObj*DELTA*JuMP.value(Ylinkpennalty)

	# prodcostvalue = JuMP.value(prodcost)
	# penaltybasevalue = DELTA*JuMP.value(basepenalty)
	# plantrevvalue = JuMP.value(plantrev)
	# penaltyYlinkvalue = DELTA*JuMP.value(Ylinkpennalty)
	
	objectivevalue = prodcostvalue - plantrevvalue + penaltybasevalue + penaltyYlinkvalue

	println("Objective:           \$", round(objectivevalue, digits=1), "\n",
	"Generation cost:     \$", round(prodcostvalue, digits=1), "\n",
	"Penalty base:        \$", round(penaltybasevalue, digits=1), "\n",
	"Penalty Ylink:      \$", round(penaltyYlinkvalue, digits=1), "\n",
	"Plant revenue:        \$", round(plantrevvalue, digits=1), "\n")
	flush(stdout)

	function getDual(dualgk::Containers.DenseAxisArray{Float64}, numV::Int64, CU::T, 
		CD::T) where {T <: Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}}
		for vs = 1:numV
			dualgk[vs] = dual(CU[vs]) - dual(CD[vs])
			# dualgk[vs] = dual(CU[vs]) + dual(CD[vs])
		end
		return dualgk
	end

	# initialize duals/gradient
	dualgkU =  Containers.DenseAxisArray{Float64}(undef, 1:numP)
	dualgkV =  Containers.DenseAxisArray{Float64}(undef, 1:numG)
	dualgkZ =  Containers.DenseAxisArray{Float64}(undef, 1:numP)
	dualgkUm =  Containers.DenseAxisArray{Float64}(undef, 1:numP)
	dualgkVm =  Containers.DenseAxisArray{Float64}(undef, 1:numG)
	# get dDt/dYt
	if Ylink.ts < timestep[end]
		dualgkU = getDual(dualgkU, numP, CUU, CUD)
		dualgkV = getDual(dualgkV, numG, CVU, CVD)
		# dualgkV .= 0.0
	else
		dualgkU .= 0.0
		dualgkV .= 0.0
	end

	if Ylink.ts > timestep[1]
		dualgkUm = getDual(dualgkUm, numP, CUUE, CUDE)
		dualgkVm = getDual(dualgkVm, numG, CVUE, CVDE)
		# dualgkVm .= 0.0
	else
		dualgkUm .= 0.0
		dualgkVm .= 0.0
	end
	dualgkZ = getDual(dualgkZ, numP, CZU, CZL)
	# combine gradients
	gk = vcat(dualgkU, dualgkZ, dualgkV, dualgkUm, dualgkVm)

	return objectivevalue, gk

end

## function for solving master problem with user defined grad & hess
function SolveMaster(o::GoOptions, N::Vector{Any}, L::Vector{Any}, T::Vector{Any},
	SSh::Vector{Any}, G::Vector{Any}, K::Vector{Any}, P::Vector{Any}, timestep::Vector{Int64}, plant::DataFrame, MVAbase::Float64, NLSolver, subNLSolver; 
		IndexSets::Union{Tuple, Nothing}=nothing,
		StartingPoint::Union{SCACOPFState,Nothing}=nothing)

	println("decomposed problem Master")
	flush(stdout)
	# check # of plants
	nplant = findall( x -> x in plant[!,:plantid], N[1][!, :Bus])
	numP = size(nplant,1)
	if numP != 0
		plant[!,:busidx] .= nplant
	end
	println("# plants: ", numP)
	flush(stdout)

	# number of generators
	numG = size(G[1], 1)
	# some commonly used numbers: 
	# end of time, number of elements in Yts and Y
	tsend = timestep[end]
	nele = 2*numP+numG
	nelet = tsend*nele

	# hessian parameter
	betak = 2.0

	# Ylink variable bounds and initial values
	function Yint()
		startv = zeros(tsend, 2*numP+numG)
		Ylb = zeros(tsend, 2*numP+numG)
		Yub = zeros(tsend, 2*numP+numG)
		for vs = 1:(2*numP+numG)
			if vs <= numP
				# YU init = pplantlb
				# startv[:,vs] .= (plant[!,:pplantub][1] - plant[!,:pplantlb][1] + plant[!,:plantRRL][1])/2
				startv[:,vs] .= 0.45
				# Ylb[:,vs] .= startv[:,vs]*0.99
				# Yub[:,vs] .= startv[:,vs]*1.01
				Ylb[:,vs] .= (plant[!,:pplantlb][1] - plant[!,:plantRRL][1]/2)
				Yub[:,vs] .= (plant[!,:pplantub][1] + plant[!,:plantRRL][1]/2)
	
			elseif vs > numP && vs <= 2*numP
				# YZ init = (demand+flex)/24
				# startv[:,vs] .= plant[!,:plantdemand][1]/24
				startv[:,vs] .= 0.0
				# Ylb[:,vs] .= -0.01+startv[:,vs]
				# Yub[:,vs] .= 0.01+startv[:,vs]
				Ylb[:,vs] .= -1*(1+plant[!,:plantflex][1]).*plant[!,:plantdemand][1]/24
				Yub[:,vs] .= (1+plant[!,:plantflex][1]).*plant[!,:plantdemand][1]/24
			else
				# YV
				startv[:,vs] .= 1*G[1][!,:Pub][vs-numP*2]/2
				# Ylb[:,vs] .= 0.99*startv[:,vs]
				# Yub[:,vs] .= 1.01*startv[:,vs]
				Ylb[:,vs] .= G[1][!,:Plb][vs-numP*2] - 1*G[1][!,:Pub][vs-numP*2]/2
				Yub[:,vs] .= G[1][!,:Pub][vs-numP*2]/2 + G[1][!,:Pub][vs-numP*2]
			end

		end

		# println("start ", startv[:,1:(2*numP+3)], " Ylb ", Ylb[:,1:(2*numP+3)], " Yub ", Yub[:,1:(2*numP+3)])
		return startv, Ylb, Yub
	end

	# function to call solver for subproblems and save outputs
	function fsub(Ylink)
		# initialize obj 
		objk = 0.0
		# reshape Ylink to 2D
		YlinkA = collect(Ylink)
		Ylink2D = reshape(YlinkA, (tsend,(2*numP+numG)))
		# initialize U,V,Z from Y
		YUlink = Ylink2D[:,1:numP]
		YZlink = Ylink2D[:,(numP+1):2*numP]
		YVlink = Ylink2D[:,(2*numP+1):end]
		# println("U inti ", YUlink, " Z init ", YZlink, " V[1:4] init ", YVlink[:,1:4])
		# initialization to store gradients
		gk = Vector{Float64}(undef,nelet)
		gkts = zeros(tsend,(2*numP+numG))
		gkmts = zeros(tsend,(numP+numG))
		
		for ts = timestep
			# get current values of Ylink 
			if ts == 1
				Ylinkts = Ylinkdecomts(YUlink[ts,:], YZlink[ts,:], YVlink[ts,:], YUlink[ts,:], YVlink[ts,:], ts)
			else
				Ylinkts = Ylinkdecomts(YUlink[ts,:], YZlink[ts,:], YVlink[ts,:], YUlink[ts-1,:], YVlink[ts-1,:], ts)
			end
			println("------------subproblem ts ", ts, "-----------------")
			flush(stdout)
			# call solver for subproblems and output objk & gradients
			objts, gts = SolveSub(o, N[ts], L[ts], T[ts], SSh[ts], G[ts], K[ts], P[ts], plant, numP, MVAbase, Ylinkts, timestep, subNLSolver)
			# add all obj from subproblems
			objk += objts
			# store dDt/dYt
			gkts[ts,:] .= gts[1:(2*numP+numG)]
			# store dDt/dUtm1 and dVt/dUtm1
			gkmts[ts,:] .= gts[(2*numP+numG+1):end]
		end
		# gUkt = dDt/dUt + dDt+1/dUt, gVkt = dDt/dVt + dDt+1/dVt, i.e., dD/dU1 = dD2/dU1 + dD1/dU1
		# gkmts[ts+1] = dDt+1/dYt, gkts[ts] = dDt/dYt
		for ts = timestep[1:(end-1)]
			gkts[ts,1:numP] += gkmts[(ts+1),1:numP]
			gkts[ts,(numP*2+1):end] += gkmts[(ts+1),(numP+1):end]
		end
		# store gradients (U,Z,V) in gk 
		gk[1:numP*tsend] = gkts[1:numP*tsend]
		gk[(numP*tsend+1):numP*2*tsend] = gkts[(numP*tsend+1):numP*2*tsend]
		gk[(numP*2*tsend+1):end] = gkts[(numP*2*tsend+1):end]
		# println("objk ", objk)
		# println("grad[1:10] ", gk[1:10])

		return objk, gk
	end

	# user defined function f objective
	function f(fgrad, Ylink...)
		println("enter f")
		flush(stdout)
		# println("fgrad[1:10] before sub ", fgrad[1:10])
		# println("Ylink ", Ylink[1:9])
		# call to get obj and grad
		fobjk, gra = fsub(Ylink)
		fgrad .= gra
		# println("fgrad[1:9] after sub ", fgrad[1:9])
		println("f objk ", fobjk)
		flush(stdout)
		return fobjk
	end
	# user defined function gradient of f 
	function ∇f(fgrad, g, Ylink...)
		println("fgrad[1:9] in ∇f ", fgrad[1:9])	
		flush(stdout)
		# println("g[5] g[9] before ", g[5], " ", g[9])	
		g .= fgrad
		println("g[15] g[20] ", g[15], " ", g[20])
		flush(stdout)
		return g
	end

	# initialize gradient g 
	fgrad = ones(nelet)
	# define functions to allow extra inputs 
	ff(Ylink...) = f(fgrad, Ylink...)
	∇ff(g, Ylink...) = ∇f(fgrad, g, Ylink...)

	# create model
	m = Model(NLSolver)

	# Master variables Y
	startv, Ylb, Yub = Yint()
	@variable(m, Ylb[ts,vs] <= Ylink[ts=1:tsend, vs=1:(2*numP+numG)] <= Yub[ts,vs], start = startv[ts,vs])

	println("start[1,1:8] ", startv[1,1:8])
	# println(" Ylb ", Ylb[:,1:(2*numP+3)]) 
	# println(" Yub ", Yub[:,1:(2*numP+3)])
	# println("start values ", start_value.(Ylink[:,1:(2*numP+3)]))
	flush(stdout)
	
	# sum_t(Z) = 0 
	for vs = 1:numP
		@constraint(m, sum(Ylink[ts, vs+numP] for ts=timestep) == 0)
	end
	# Ut and Zt: demand/(T*unitProd) - RRL/2 <= Ut - Zt/unitProd <= (demand+flex)/(T*unitProd) + RRL/2
	for ts = timestep
		@constraint(m, [vs=1:numP], Ylink[ts,vs] - Ylink[ts,numP+vs]/plant[!,:plantprod][vs] <= plant[!,:plantRRL][vs]/2 + (1+plant[!,:plantflex][vs])*plant[!,:plantdemand][vs]/24/plant[!,:plantprod][vs])
		@constraint(m, [vs=1:numP], Ylink[ts,vs] - Ylink[ts,numP+vs]/plant[!,:plantprod][vs] >= -plant[!,:plantRRL][vs]/2 + plant[!,:plantdemand][vs]/24/plant[!,:plantprod][vs])
	end
	# Ut-1 and Zt: same relationship
	for ts = timestep[1:(end-1)]
		@constraint(m, [vs=1:numP], Ylink[ts,vs] - Ylink[ts+1,numP+vs]/plant[!,:plantprod][vs] <= plant[!,:plantRRL][vs]/2 + (1+plant[!,:plantflex][vs])*plant[!,:plantdemand][vs]/24/plant[!,:plantprod][vs])
		@constraint(m, [vs=1:numP], Ylink[ts,vs] - Ylink[ts+1,numP+vs]/plant[!,:plantprod][vs] >= -plant[!,:plantRRL][vs]/2 + plant[!,:plantdemand][vs]/24/plant[!,:plantprod][vs])
	end
	# -RRL <= Ut - Ut-1 <= RRL
	for ts = timestep[1:(end-1)]
		@constraint(m, [vs=1:numP], -plant[!,:plantRRL][vs] <= Ylink[ts+1,vs] - Ylink[ts,vs] <= plant[!,:plantRRL][vs])
	end

	# -RRL_gen <= Vt - Vt-1 <= RRL_gen
	for ts = timestep[1:(end-1)]
		@constraint(m, [vs=1:numG], -1*G[1][!,:Pub][vs] <= Ylink[ts+1,vs+2*numP] - Ylink[ts,vs+2*numP] <= 1*G[1][!,:Pub][vs])
	end

	# register(m, :my_obj, nelet, fs, ∇fs, ∇²f)
	# @NLobjective(m, Min, my_obj(Ylink...))

	println("before register ")
	flush(stdout)
	# register functions 
	register(m, :my_objf, nelet, ff, ∇ff)
	objtest = zero(JuMP.QuadExpr)
	for ts = timestep
		add_to_expression!(objtest, sum(Ylink[ts,vs]^2 for vs=1:(2*numP+numG)))
	end
	# declare objective
	@NLobjective(m, Min, my_objf(Ylink...))
	# @NLobjective(m, Min, objtest)
	println("after declare obj ")
	flush(stdout)
	## attempt to solve SCACOPF
	JuMP.optimize!(m)
 
        # termination_status, primal_status, and dual_status: http://www.juliaopt.org/JuMP.jl/dev/solutions/
        solvestatusprimal = JuMP.primal_status(m)
        solvestatus = JuMP.termination_status(m)
        println("JuMP optimize MOI termination_status: ", solvestatus, " ", solvestatusprimal)
	if (solvestatusprimal != MOI.FEASIBLE_POINT) && (solvestatusprimal != MOI.NEARLY_FEASIBLE_POINT)
		error("solver failed to find a feasible solution")
	end

end

## function for solving each subproblem with regularization terms 
function SolveSubr(o::GoOptions, N::DataFrame, L::DataFrame, T::DataFrame,
	SSh::DataFrame, G::DataFrame, K::DataFrame, P::DataFrame, plant::DataFrame, numP::Int64, MVAbase::Float64, 
	Ylink::Ylinkdecomts, Lk::Vector{Float64} , timestep::Vector{Int64}, NLSolver; 
		IndexSets::Union{Tuple, Nothing}=nothing,
		StartingPoint::Union{SCACOPFState,Nothing}=nothing)

	if numP > 3
		println("U[1:2],V[1:3],Z[1:2] ", Ylink.YUk[1:2], " ", Ylink.YVk[1:3], " ", Ylink.YZk[1:2])
	else
		println("U, V[1:3], Z ", Ylink.YUk, " ", Ylink.YVk[1:3], " ", Ylink.YZk)
	end
	flush(stdout)
	# compute index and sets for performant formulation
	L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn, K_outidx =
	build_or_split_indexsets(IndexSets, N, L, T, SSh, G, K)

	RefBus = G_Nidx[argmax(G[!,:Pub])]	# bus with the biggest generator is reference

	# (some) output for options used
	println("SolveSCACOPF with Options: CouplingMode=", o.CouplingMode, 
  	" Smoothing=", o.SmoothingParamCoupling,
	" ThermalLimits=", o.ThermalLimits)
	

	if StartingPoint==nothing
		# compute SC AC OPF starting point 
		v0_n, theta0_n, p0_li, q0_li, p0_ti, q0_ti, b0_s, p0_g, q0_g, 
		pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti =
 	 InitialBaseCaseSolution( N, L, T, SSh, G, 
						   IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
								SShn, Gn, K_outidx))
		theta0_n .-= theta0_n[RefBus]
	else
		## base case
		v0_n, theta0_n, b0_s, p0_g, q0_g = getBaseStates(StartingPoint)
		p0_li, q0_li, p0_ti, q0_ti, pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti = 
 	 	BaseCaseSolution(v0_n, theta0_n, b0_s, p0_g, q0_g,
				   N, L, T, SSh, G, 
				   IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
						SShn, Gn, K_outidx))
		println("HOT Started")
	end

	numG = size(G,1)
	# create model
	m = Model(NLSolver)

	# ACOPF variables
	@variable(m, N[!,:Vlb][n] <= v_n[n=1:size(N, 1)] <= N[!,:Vub][n], start = v0_n[n])
	@variable(m, theta_n[n=1:size(N, 1)], start = theta0_n[n])
	@variable(m, p_li[l=1:size(L, 1), i=1:2], start = p0_li[l,i])
	@variable(m, q_li[l=1:size(L, 1), i=1:2], start = q0_li[l,i])
	@variable(m, p_ti[t=1:size(T, 1), i=1:2], start = p0_ti[t,i])
	@variable(m, q_ti[t=1:size(T, 1), i=1:2], start = q0_ti[t,i])
	@variable(m, SSh[!,:Blb][s] <= b_s[s=1:size(SSh, 1)] <= SSh[!,:Bub][s], start = b0_s[s])
	@variable(m, G[!,:Plb][g] <= p_g[g=1:size(G, 1)] <= G[!,:Pub][g], start = p0_g[g])
	@variable(m, G[!,:Qlb][g] <= q_g[g=1:size(G, 1)] <= G[!,:Qub][g], start = q0_g[g])
	@variable(m, pslackm_n[n=1:size(N, 1)] >= 0, start = pslack0m_n[n])
	@variable(m, pslackp_n[n=1:size(N, 1)] >= 0, start = pslack0p_n[n])
	@variable(m, qslackm_n[n=1:size(N, 1)] >= 0, start = qslack0m_n[n])
	@variable(m, qslackp_n[n=1:size(N, 1)] >= 0, start = qslack0p_n[n])
	@variable(m, sslack_li[l=1:size(L, 1), i=1:2] >= 0, start = sslack0_li[l,i])
	@variable(m, sslack_ti[t=1:size(T, 1), i=1:2] >= 0, start = sslack0_ti[t,i])

	# slack variables for Ylink 
	@variable(m, YUslack[vs=1:numP*2] >= 0, start = 0.0)
	@variable(m, YZslack[vs=1:numP*2] >= 0, start = 0.0)
	@variable(m, YVslack[vs=1:numG*2] >= 0, start = 0.0)
	@variable(m, YUmslack[vs=1:numP*2] >= 0, start = 0.0)
	@variable(m, YVmslack[vs=1:numG*2] >= 0, start = 0.0)

	# slack variables for regularization
	@variable(m, YUr[vs=1:numP], start = 0.0)
	@variable(m, YZr[vs=1:numP], start = 0.0)
	@variable(m, YVr[vs=1:numG], start = 0.0)
	@variable(m, YUmr[vs=1:numP], start = 0.0)
	@variable(m, YVmr[vs=1:numG], start = 0.0)
	rho_R = 1e-4
	function findrho(LkY::Vector{Float64}, YL::Vector{Float64}, YU::Vector{Float64})
		n = size(LkY,1)
		YLU = zeros(n)
		rhomin = zeros(n)
		for i = 1:n
			YLU[i] = findmax([abs(YL[i]), abs(YU[i])])[1]
			if abs(LkY[i]) > 1e-6
				rhomin[i] = rho_R*YLU[i]/LkY[i]
			else
				rhomin[i] = 1e-8*YLU[i]
			end
		end
		rhobase = YLU.*1e-8
		rmin = findmin(rhomin)[1]
		rhoY = sum(findmax([rhobase; rmin])[1])
		println("rhoY ", rhoY)
		return rhoY
	end

	# bidding entity variables
	@variable(m, plant[!,:pplantlb][np] <= pplant[np=1:numP] <= plant[!,:pplantub][np])
	@variable(m, plprod[np=1:numP] >= 0, start = 0.0)
	
	# bidding plant constraints
	# calculate plant prod from energy consumption at corresponding bus
	@constraint(m, [np=1:numP], plprod[np] == plant[!,:plantprod][np].*pplant[np])
	
	# ramping constraint for plant production 
	rhoU = findrho(Lk[1:numP], -plant[!,:plantRRL]./2, plant[!,:plantRRL]./2)
	rhoUm = findrho(Lk[(numP*2+numG+1):(numP*3+numG)], -plant[!,:plantRRL]./2, plant[!,:plantRRL]./2)
	if Ylink.ts < timestep[end]
		@constraint(m, CUU[np=1:numP], pplant[np] + rhoU*YUr[np] <= Ylink.YUk[np]+plant[!,:plantRRL][np]/2 + YUslack[np])
		@constraint(m, CUD[np=1:numP], Ylink.YUk[np] - plant[!,:plantRRL][np]/2 - YUslack[np+numP] <= pplant[np] + rhoU*YUr[np])
	end
	if Ylink.ts > timestep[1]
		@constraint(m, CUUE[np=1:numP], pplant[np] + rhoUm*YUmr[np] <= Ylink.YUkm[np]+plant[!,:plantRRL][np]/2 + YUmslack[np])
		@constraint(m, CUDE[np=1:numP], Ylink.YUkm[np] - plant[!,:plantRRL][np]/2 - YUmslack[np+numP] <= pplant[np] + rhoUm*YUmr[np] )
	end
	# temporal constraint: plant production within demand and flex*demand 
	Du = (plant[!,:plantflex].+1).*(plant[!,:plantdemand]./24)
	Dl = plant[!,:plantdemand]./24
	rhoZ = findrho(Lk[(numP+1):numP*2], Dl, Du)
	@constraint(m, CZU[np=1:numP], plprod[np] - Ylink.YZk[np] + rhoZ*YZr[np] <=  Du[np] + YZslack[np])
	@constraint(m, CZL[np=1:numP], Dl[np] - YZslack[np+numP] <= plprod[np] - Ylink.YZk[np] + rhoZ*YZr[np])

	# ramping constraints for power generation
	rhoV = findrho(Lk[(numP*2+1):(numP*2+numG)], -G[!,:Pub]./2, G[!,:Pub]./2)
	rhoVm = findrho(Lk[(numP*3+numG+1):end], -G[!,:Pub]./2, G[!,:Pub]./2)
	if Ylink.ts < timestep[end]
		@constraint(m, CVU[g=1:numG], p_g[g] + rhoV*YVr[g] <= Ylink.YVk[g] + 1*G[!,:Pub][g]/2 + YVslack[g])
		@constraint(m, CVD[g=1:numG], Ylink.YVk[g] - 1*G[!,:Pub][g]/2 - YVslack[g+numG] <= p_g[g] + rhoV*YVr[g])
	end
	if Ylink.ts > timestep[1]
		@constraint(m, CVUE[g=1:numG], p_g[g] + rhoVm*YVmr[g] <= Ylink.YVkm[g] + 1*G[!,:Pub][g]/2 + YVmslack[g])
		@constraint(m, CVDE[g=1:numG], Ylink.YVkm[g] - 1*G[!,:Pub][g]/2 - YVmslack[g+numG] <= p_g[g] + rhoVm*YVmr[g])
	end
	
	#active power balance with LMP=dual
	@constraint(m, PBPC[np=1:numP], sum(p_g[g] for g=Gn[plant[!,:busidx][np]]) - N[!,:Pd][plant[!,:busidx][np]] - N[!,:Gsh][plant[!,:busidx][np]]*v_n[plant[!,:busidx][np]]^2 -
			pplant[np] - 
			sum(p_li[Lidxn[plant[!,:busidx][np]][lix],Lin[plant[!,:busidx][np]][lix]] for lix=1:length(Lidxn[plant[!,:busidx][np]])) -
			sum(p_ti[Tidxn[plant[!,:busidx][np]][tix],Tin[plant[!,:busidx][np]][tix]] for tix=1:length(Tidxn[plant[!,:busidx][np]])) ==
			pslackp_n[plant[!,:busidx][np]] - pslackm_n[plant[!,:busidx][np]])

	nbusn = findall( x -> !(x in plant[!,:plantid]), N[!, :Bus])
	@constraint(m, PB[n=nbusn], sum(p_g[g] for g=Gn[n]) - N[!,:Pd][n] - N[!,:Gsh][n]*v_n[n]^2 -
	sum(p_li[Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n])) -
	sum(p_ti[Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n])) ==
	pslackp_n[n] - pslackm_n[n])

	# fix angle at reference bus to zero
	JuMP.fix(theta_n[RefBus], 0.0, force=true)

	# add power flow constraints
	addpowerflowcons!(m, v_n, theta_n, p_li, q_li, p_ti, q_ti, b_s, p_g, q_g,
	pslackm_n, pslackp_n, qslackm_n, qslackp_n, sslack_li, sslack_ti,
	N, L, T, (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
	SShn, Gn, K_outidx), o.ThermalLimits)

	# register approximation functions
	sclamp = nothing
	sstep = nothing

	## objective function
	
	# production cost
	prodcost = zero(JuMP.AffExpr)
	plantrev = zero(JuMP.AffExpr)
	slackpenalty = zero(JuMP.QuadExpr)

	for g = 1:size(G, 1)
	    npairs = length(G[!,:CostPi][g])
            @assert sort(G[!,:CostPi][g]) == G[!,:CostPi][g]
 
            t_h0 = getstartcostcoefficients(JuMP.start_value(p_g[g]), G[!,:CostPi][g], G[!,:Plb][g], G[!,:Pub][g])
            t_h = @variable(m, [i=1:npairs], lower_bound=0.0, start=t_h0[i])
            
 	    @constraint(m, sum(t_h[h] for h=1:npairs) == 1.0)
	    @constraint(m, p_g[g] == sum(G[!,:CostPi][g][h]*t_h[h] for h=1:npairs))
            
	    add_to_expression!(prodcost, @expression(m, sum(G[!,:CostCi][g][h]*t_h[h] for h=1:npairs)))
	end
	
	# base case penalties
	basepenalty = quadpenaltyfunctions(pslackm_n, pslackp_n, qslackm_n, qslackp_n,
				sslack_li, sslack_ti, P)

	# plant prod revenue
	for np=1:numP
		add_to_expression!(plantrev, pplant[np].*plant[!,:plantrev][np])
	end

	# Ylink penalties
	Ylinkpennalty = Ylinkpen(YUslack, YZslack, YVslack, YUmslack, YVmslack)

	# slack penalties
	YUZrs = @expression(m, (YUr-Lk[1:numP]).^2+(YZr-Lk[(numP+1):numP*2]).^2+(YUmr-Lk[(numP*2+numG+1):(numP*3+numG)]).^2)
	YVrs = @expression(m, (YVr-Lk[(numP*2+1):(numP*2+numG)]).^2+(YVmr-Lk[(numP*3+numG+1):end]).^2)

	slackrs = @expression(m, sum(YUZrs) + sum(YVrs))
	add_to_expression!(slackpenalty, 0.5, slackrs)

	# declare objective
	@objective(m, Min, (prodcost + DELTA*basepenalty - plantrev + DELTA*Ylinkpennalty + slackpenalty))
	
	## attempt to solve ACOPF
	JuMP.optimize!(m)

	    # termination_status, primal_status, and dual_status: http://www.juliaopt.org/JuMP.jl/dev/solutions/
        solvestatusprimal = JuMP.primal_status(m)
        solvestatus = JuMP.termination_status(m)
        println("JuMP optimize MOI termination_status: ", solvestatus, " ", solvestatusprimal)
	if (solvestatusprimal != MOI.FEASIBLE_POINT) && (solvestatusprimal != MOI.NEARLY_FEASIBLE_POINT)
		error("solver failed to find a feasible solution")
	end
	
	## objective breakdown
	prodcostvalue = JuMP.value(prodcost)
	penaltybasevalue = DELTA*JuMP.value(basepenalty)
	plantrevvalue = JuMP.value(plantrev)
	penaltyYlinkvalue = DELTA*JuMP.value(Ylinkpennalty)
	slackpenvalue = JuMP.value(slackpenalty)
	objectivevalue = prodcostvalue + penaltybasevalue - plantrevvalue + penaltyYlinkvalue + slackpenvalue

	println("Objective:           \$", round(objectivevalue, digits=1), "\n",
	"Generation cost:     \$", round(prodcostvalue, digits=1), "\n",
	"Penalty base:        \$", round(penaltybasevalue, digits=1), "\n",
	"Penalty Ylink:      \$", round(penaltyYlinkvalue, digits=1), "\n",
	"Penalty Slack:      \$", round(slackpenvalue, digits=1), "\n",
	"Plant revenue:        \$", round(plantrevvalue, digits=1), "\n")
	flush(stdout)

	function getDual(dualgk::Containers.DenseAxisArray{Float64}, numV::Int64, CU::T, 
		CD::T) where {T <: Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}}
		for vs = 1:numV
			dualgk[vs] = dual(CU[vs]) + dual(CD[vs])
		end
		return dualgk
	end

	# initialize duals/gradient
	dualgkU =  Containers.DenseAxisArray{Float64}(undef, 1:numP)
	dualgkV =  Containers.DenseAxisArray{Float64}(undef, 1:numG)
	dualgkZ =  Containers.DenseAxisArray{Float64}(undef, 1:numP)
	dualgkUm =  Containers.DenseAxisArray{Float64}(undef, 1:numP)
	dualgkVm =  Containers.DenseAxisArray{Float64}(undef, 1:numG)
	# get dDt/dYt
	if Ylink.ts < timestep[end]
		dualgkU = getDual(dualgkU, numP, CUU, CUD)
		dualgkV = getDual(dualgkV, numG, CVU, CVD)
		# dualgkV .= 0.0
	else
		dualgkU .= 0.0
		dualgkV .= 0.0
	end

	if Ylink.ts > timestep[1]
		dualgkUm = getDual(dualgkUm, numP, CUUE, CUDE)
		dualgkVm = getDual(dualgkVm, numG, CVUE, CVDE)
		# dualgkVm .= 0.0
	else
		dualgkUm .= 0.0
		dualgkVm .= 0.0
	end
	dualgkZ = getDual(dualgkZ, numP, CZU, CZL)
	# combine gradients
	gk = vcat(dualgkU, dualgkZ, dualgkV, dualgkUm, dualgkVm)

	return objectivevalue, gk

end

## function for solving each subproblem
function SolveSubSQP(o::GoOptions, N::DataFrame, L::DataFrame, T::DataFrame,
	SSh::DataFrame, G::DataFrame, K::DataFrame, P::DataFrame, plant::DataFrame, numP::Int64, MVAbase::Float64, Ylink::Ylinkdecomts, timestep::Vector{Int64}, NLSolver; 
		IndexSets::Union{Tuple, Nothing}=nothing,
		StartingPoint::Union{SCACOPFState,Nothing}=nothing)

	if numP > 3
		println("U[1:2],V[1:3],Z[1:2] ", Ylink.YUk[1:2], " ", Ylink.YVk[1:3], " ", Ylink.YZk[1:2])
	else
		println("U, V[1:3], Z ", Ylink.YUk, " ", Ylink.YVk[1:3], " ", Ylink.YZk)
	end
	flush(stdout)
	# compute index and sets for performant formulation
	L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn, K_outidx =
	build_or_split_indexsets(IndexSets, N, L, T, SSh, G, K)

	RefBus = G_Nidx[argmax(G[!,:Pub])]	# bus with the biggest generator is reference

	# (some) output for options used
	println("SolveSCACOPF with Options: CouplingMode=", o.CouplingMode, 
  	" Smoothing=", o.SmoothingParamCoupling,
	" ThermalLimits=", o.ThermalLimits)
	

	if StartingPoint==nothing
		# compute SC AC OPF starting point 
		v0_n, theta0_n, p0_li, q0_li, p0_ti, q0_ti, b0_s, p0_g, q0_g, 
		pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti =
 	 InitialBaseCaseSolution( N, L, T, SSh, G, 
						   IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
								SShn, Gn, K_outidx))
		theta0_n .-= theta0_n[RefBus]
	else
		## base case
		v0_n, theta0_n, b0_s, p0_g, q0_g = getBaseStates(StartingPoint)
		p0_li, q0_li, p0_ti, q0_ti, pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti = 
 	 	BaseCaseSolution(v0_n, theta0_n, b0_s, p0_g, q0_g,
				   N, L, T, SSh, G, 
				   IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
						SShn, Gn, K_outidx))
		println("HOT Started")
	end

	numG = size(G,1)
	# create model
	m = Model(NLSolver)

	# ACOPF variables
	@variable(m, N[!,:Vlb][n] <= v_n[n=1:size(N, 1)] <= N[!,:Vub][n], start = v0_n[n])
	@variable(m, theta_n[n=1:size(N, 1)], start = theta0_n[n])
	@variable(m, p_li[l=1:size(L, 1), i=1:2], start = p0_li[l,i])
	@variable(m, q_li[l=1:size(L, 1), i=1:2], start = q0_li[l,i])
	@variable(m, p_ti[t=1:size(T, 1), i=1:2], start = p0_ti[t,i])
	@variable(m, q_ti[t=1:size(T, 1), i=1:2], start = q0_ti[t,i])
	@variable(m, SSh[!,:Blb][s] <= b_s[s=1:size(SSh, 1)] <= SSh[!,:Bub][s], start = b0_s[s])
	@variable(m, G[!,:Plb][g] <= p_g[g=1:size(G, 1)] <= G[!,:Pub][g], start = p0_g[g])
	@variable(m, G[!,:Qlb][g] <= q_g[g=1:size(G, 1)] <= G[!,:Qub][g], start = q0_g[g])
	@variable(m, pslackm_n[n=1:size(N, 1)] >= 0, start = pslack0m_n[n])
	@variable(m, pslackp_n[n=1:size(N, 1)] >= 0, start = pslack0p_n[n])
	@variable(m, qslackm_n[n=1:size(N, 1)] >= 0, start = qslack0m_n[n])
	@variable(m, qslackp_n[n=1:size(N, 1)] >= 0, start = qslack0p_n[n])
	@variable(m, sslack_li[l=1:size(L, 1), i=1:2] >= 0, start = sslack0_li[l,i])
	@variable(m, sslack_ti[t=1:size(T, 1), i=1:2] >= 0, start = sslack0_ti[t,i])

	# slack variables for Ylink 
	@variable(m, YUslack[vs=1:numP*2] >= 0, start = 0.0)
	@variable(m, YZslack[vs=1:numP*2] >= 0, start = 0.0)
	@variable(m, YVslack[vs=1:numG*2] >= 0, start = 0.0)
	@variable(m, YUmslack[vs=1:numP*2] >= 0, start = 0.0)
	@variable(m, YVmslack[vs=1:numG*2] >= 0, start = 0.0)

	# bidding entity variables
	@variable(m, plant[!,:pplantlb][np] <= pplant[np=1:numP] <= plant[!,:pplantub][np])
	@variable(m, plprod[np=1:numP] >= 0, start = 0.0)
	
	# bidding plant constraints
	# calculate plant prod from energy consumption at corresponding bus
	@constraint(m, [np=1:numP], plprod[np] == plant[!,:plantprod][np].*pplant[np])
	# ramping constraint for plant production 
	if Ylink.ts < timestep[end]
		# @constraint(m, CUU[np=1:numP], pplant[np] <= Ylink.YUk[np]+plant[!,:plantRRL][np]/2)
		# @constraint(m, CUD[np=1:numP], Ylink.YUk[np] - plant[!,:plantRRL][np]/2 <= pplant[np])
		@constraint(m, CUU[np=1:numP], pplant[np] <= Ylink.YUk[np]+plant[!,:plantRRL][np]/2 + YUslack[np])
		@constraint(m, CUD[np=1:numP], Ylink.YUk[np] - plant[!,:plantRRL][np]/2 - YUslack[np+numP] <= pplant[np])
	end
	if Ylink.ts > timestep[1]
		# @constraint(m, CUUE[np=1:numP], pplant[np] <= Ylink.YUkm[np]+plant[!,:plantRRL][np]/2)
		# @constraint(m, CUDE[np=1:numP], Ylink.YUkm[np]/2 - plant[!,:plantRRL][np] <= pplant[np])
		@constraint(m, CUUE[np=1:numP], pplant[np] <= Ylink.YUkm[np]+plant[!,:plantRRL][np]/2 + YUmslack[np])
		@constraint(m, CUDE[np=1:numP], Ylink.YUkm[np] - plant[!,:plantRRL][np]/2 - YUmslack[np+numP] <= pplant[np])
	end
	# temporal constraint: plant production within demand and flex*demand 
	Du = (plant[!,:plantflex].+1).*(plant[!,:plantdemand]./24)
	Dl = plant[!,:plantdemand]./24
	@constraint(m, CZU[np=1:numP], plprod[np] - Ylink.YZk[np] <= Du[np] + YZslack[np])
	@constraint(m, CZL[np=1:numP], Dl[np] - YZslack[np+numP] <= plprod[np] - Ylink.YZk[np])
	# @constraint(m, CZU[np=1:numP], plprod[np] - Ylink.YZk[np] <= (1+plant[!,:plantflex][1]).*plant[!,:plantdemand][1]/24)
	# @constraint(m, CZL[np=1:numP], plant[!,:plantdemand][1]/24 <= plprod[np] - Ylink.YZk[np])

	# ramping constraints for power generation
	if Ylink.ts < timestep[end]
		@constraint(m, CVU[g=1:numG], p_g[g] <= Ylink.YVk[g] + 1*G[!,:Pub][g]/2 + YVslack[g])
		@constraint(m, CVD[g=1:numG], Ylink.YVk[g] - 1*G[!,:Pub][g]/2 - YVslack[g+numG] <= p_g[g])
		# @constraint(m, CVU[g=1:numG], p_g[g] <= Ylink.YVk[g] + 1*G[!,:Pub][g]/2)
		# @constraint(m, CVD[g=1:numG], Ylink.YVk[g] - 1*G[!,:Pub][g]/2 <= p_g[g])
	end
	if Ylink.ts > timestep[1]
		@constraint(m, CVUE[g=1:numG], p_g[g] <= Ylink.YVkm[g] + 1*G[!,:Pub][g]/2 + YVmslack[g])
		@constraint(m, CVDE[g=1:numG], Ylink.YVkm[g] - 1*G[!,:Pub][g]/2 - YVmslack[g+numG] <= p_g[g])
		# @constraint(m, CVUE[g=1:numG], p_g[g] <= Ylink.YVkm[g] + 1*G[!,:Pub][g]/2)
		# @constraint(m, CVDE[g=1:numG], Ylink.YVkm[g] - 1*G[!,:Pub][g]/2 <= p_g[g])
	end

	#active power balance with LMP=dual
	@constraint(m, PBPC[np=1:numP], sum(p_g[g] for g=Gn[plant[!,:busidx][np]]) - N[!,:Pd][plant[!,:busidx][np]] - N[!,:Gsh][plant[!,:busidx][np]]*v_n[plant[!,:busidx][np]]^2 -
			pplant[np] - 
			sum(p_li[Lidxn[plant[!,:busidx][np]][lix],Lin[plant[!,:busidx][np]][lix]] for lix=1:length(Lidxn[plant[!,:busidx][np]])) -
			sum(p_ti[Tidxn[plant[!,:busidx][np]][tix],Tin[plant[!,:busidx][np]][tix]] for tix=1:length(Tidxn[plant[!,:busidx][np]])) ==
			pslackp_n[plant[!,:busidx][np]] - pslackm_n[plant[!,:busidx][np]])

	nbusn = findall( x -> !(x in plant[!,:plantid]), N[!, :Bus])
	@constraint(m, PB[n=nbusn], sum(p_g[g] for g=Gn[n]) - N[!,:Pd][n] - N[!,:Gsh][n]*v_n[n]^2 -
	sum(p_li[Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n])) -
	sum(p_ti[Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n])) ==
	pslackp_n[n] - pslackm_n[n])

	# fix angle at reference bus to zero
	JuMP.fix(theta_n[RefBus], 0.0, force=true)

	# add power flow constraints
	addpowerflowcons!(m, v_n, theta_n, p_li, q_li, p_ti, q_ti, b_s, p_g, q_g,
	pslackm_n, pslackp_n, qslackm_n, qslackp_n, sslack_li, sslack_ti,
	N, L, T, (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
	SShn, Gn, K_outidx), o.ThermalLimits)

	# register approximation functions
	sclamp = nothing
	sstep = nothing

	## objective function
	
	# production cost
	prodcost = zero(JuMP.AffExpr)
	# basepenalty = zero(JuMP.QuadExpr)
	plantrev = zero(JuMP.AffExpr)
	regpenalty = zero(JuMP.QuadExpr)

	for g = 1:size(G, 1)
	    npairs = length(G[!,:CostPi][g])
            @assert sort(G[!,:CostPi][g]) == G[!,:CostPi][g]
 
            t_h0 = getstartcostcoefficients(JuMP.start_value(p_g[g]), G[!,:CostPi][g], G[!,:Plb][g], G[!,:Pub][g])
            t_h = @variable(m, [i=1:npairs], lower_bound=0.0, start=t_h0[i])
            
 	    @constraint(m, sum(t_h[h] for h=1:npairs) == 1.0)
	    @constraint(m, p_g[g] == sum(G[!,:CostPi][g][h]*t_h[h] for h=1:npairs))
            
	    add_to_expression!(prodcost, @expression(m, sum(G[!,:CostCi][g][h]*t_h[h] for h=1:npairs)))
	end
	
	# base case penalties
	basepenalty = quadpenaltyfunctions(pslackm_n, pslackp_n, qslackm_n, qslackp_n,
				sslack_li, sslack_ti, P)

	# plant prod revenue
	for np=1:numP
		add_to_expression!(plantrev, pplant[np].*plant[!,:plantrev][np])
	end

	# Ylink penalties
	Ylinkpennalty = Ylinkpen(YUslack, YZslack, YVslack, YUmslack, YVmslack)

	# regularization/stablization penalties
	YUrs = @expression(m, (pplant-Ylink.YUk-plant[!,:plantRRL]./2-YUslack[1:numP]).^2 + (Ylink.YUk-plant[!,:plantRRL]./2-YUslack[(1+numP):numP*2]-pplant).^2)
	YUmrs = @expression(m, (pplant-Ylink.YUkm-plant[!,:plantRRL]./2-YUmslack[1:numP]).^2 + (Ylink.YUkm-plant[!,:plantRRL]./2-YUmslack[(1+numP):numP*2]-pplant).^2)
	YZrs = @expression(m, (plprod-Ylink.YZk .-Du-YZslack[1:numP]).^2 + (-YZslack[(1+numP):numP*2].+Dl-plprod+Ylink.YZk).^2)
	Urs = @expression(m, 1e-4/numP*sum(YUrs+YUmrs))
	Zrs = @expression(m, 1e-1/numP*sum(YZrs))

	YVrs = @expression(m, (p_g-Ylink.YVk-1*G[!,:Pub]./2-YVslack[1:numG]).^2 + (Ylink.YVk-1*G[!,:Pub]./2-YVslack[(1+numG):numG*2]-p_g).^2)
	YVmrs = @expression(m, (p_g-Ylink.YVkm-1*G[!,:Pub]./2-YVmslack[1:numG]).^2 + (Ylink.YVkm-1*G[!,:Pub]./2-YVmslack[(1+numG):numG*2]-p_g).^2)
	Vrs = @expression(m, 1e-4*sum(YVrs+YVmrs))

	slackrs = @expression(m, Urs+Zrs+Vrs)
	rho = 1.0
	add_to_expression!(regpenalty, rho, slackrs)



	# prodcost = @expression(m, prodcost/1024)
	# plantrev = @expression(m, plantrev/1024)
	# declare objective
	@objective(m, Min, prodcost - plantrev + DELTA*basepenalty + DELTA*Ylinkpennalty + regpenalty)
	
	## attempt to solve ACOPF
	JuMP.optimize!(m)

	    # termination_status, primal_status, and dual_status: http://www.juliaopt.org/JuMP.jl/dev/solutions/
        solvestatusprimal = JuMP.primal_status(m)
        solvestatus = JuMP.termination_status(m)
        println("JuMP optimize MOI termination_status: ", solvestatus, " ", solvestatusprimal)
	if (solvestatusprimal != MOI.FEASIBLE_POINT) && (solvestatusprimal != MOI.NEARLY_FEASIBLE_POINT)
		error("solver failed to find a feasible solution")
	end
	
	## objective breakdown
	prodcostvalue = JuMP.value(prodcost)
	penaltybasevalue = DELTA*JuMP.value(basepenalty)
	plantrevvalue = JuMP.value(plantrev)
	penaltyYlinkvalue = DELTA*JuMP.value(Ylinkpennalty)
	regpenvalue = JuMP.value(regpenalty)
	objectivevalue = prodcostvalue - plantrevvalue + penaltybasevalue + penaltyYlinkvalue + regpenvalue

	println("Objective:           \$", round(objectivevalue, digits=1), "\n",
	"Generation cost:     \$", round(prodcostvalue, digits=1), "\n",
	"Penalty base:        \$", round(penaltybasevalue, digits=1), "\n",
	"Penalty Ylink:      \$", round(penaltyYlinkvalue, digits=1), "\n",
	"Penalty reg:      \$", round(regpenvalue, digits=1), "\n",
	"Plant revenue:        \$", round(plantrevvalue, digits=1), "\n")
	flush(stdout)

	Upenvalue = JuMP.value(Urs)
	Zpenvalue = JuMP.value(Zrs)
	Vpenvalue = JuMP.value(Vrs)

	println("Ureg pen:      \$", round(Upenvalue, digits=1), "\n",
	"Zreg pen:      \$", round(Zpenvalue, digits=1), "\n",
	"Vreg pen:      \$", round(Vpenvalue, digits=1), "\n")

	function getDual(dualgk::Containers.DenseAxisArray{Float64}, numV::Int64, CU::T, 
		CD::T) where {T <: Vector{ConstraintRef{Model, MathOptInterface.ConstraintIndex{MathOptInterface.ScalarAffineFunction{Float64}, MathOptInterface.LessThan{Float64}}, ScalarShape}}}
		for vs = 1:numV
			dualgk[vs] = dual(CU[vs]) + dual(CD[vs])
		end
		return dualgk
	end

	# initialize duals/gradient
	dualgkU =  Containers.DenseAxisArray{Float64}(undef, 1:numP)
	dualgkV =  Containers.DenseAxisArray{Float64}(undef, 1:numG)
	dualgkZ =  Containers.DenseAxisArray{Float64}(undef, 1:numP)
	dualgkUm =  Containers.DenseAxisArray{Float64}(undef, 1:numP)
	dualgkVm =  Containers.DenseAxisArray{Float64}(undef, 1:numG)
	# get dDt/dYt
	if Ylink.ts < timestep[end]
		dualgkU = getDual(dualgkU, numP, CUU, CUD)
		dualgkV = getDual(dualgkV, numG, CVU, CVD)
		# dualgkV .= 0.0
	else
		dualgkU .= 0.0
		dualgkV .= 0.0
	end

	if Ylink.ts > timestep[1]
		dualgkUm = getDual(dualgkUm, numP, CUUE, CUDE)
		dualgkVm = getDual(dualgkVm, numG, CVUE, CVDE)
		# dualgkVm .= 0.0
	else
		dualgkUm .= 0.0
		dualgkVm .= 0.0
	end
	dualgkZ = getDual(dualgkZ, numP, CZU, CZL)
	# combine gradients
	gk = vcat(dualgkU, dualgkZ, dualgkV, dualgkUm, dualgkVm)

	return objectivevalue, gk

end

## function to slove the master problem via SQP
function SolveMasterSQP(o::GoOptions, N::Vector{Any}, L::Vector{Any}, T::Vector{Any},
	SSh::Vector{Any}, G::Vector{Any}, K::Vector{Any}, P::Vector{Any}, timestep::Vector{Int64}, plant::DataFrame, MVAbase::Float64, SQPNLSolver, NLSolver; 
		IndexSets::Union{Tuple, Nothing}=nothing,
		StartingPoint::Union{SCACOPFState,Nothing}=nothing)

	println("decomposed problem Master SQP")
	# check # of plants
	nplant = findall( x -> x in plant[!,:plantid], N[1][!, :Bus])
	numP = size(nplant,1)
	if numP != 0
		plant[!,:busidx] .= nplant
	end
	println("# plants: ", numP)
	flush(stdout)
	# number of generators
	numG = size(G[1], 1)
	tsend = timestep[end]

	
	# initialize gk=df/dYk, matrix to store Yk
	Yk = zeros(tsend, numP*2+numG)
	gkts = zeros(tsend, 2*numP+numG)
	gktsm = zeros(tsend, numP+numG)
	gkm1 = zeros(tsend, 3*numP+2*numG)

	# SQP parameters
	# hessian constant
	alphak = 2
	# step size
	betak = 1.0
	

	# Ylink variable bounds and initial values
	function Yint()
		startv = zeros(tsend, 2*numP+numG)
		Ylb = zeros(tsend, 2*numP+numG)
		Yub = zeros(tsend, 2*numP+numG)
		for vs = 1:(2*numP+numG)
			if vs <= numP
				# YU init = pplantlb
				startv[:,vs] .= (plant[!,:pplantub][1] - plant[!,:pplantlb][1] + plant[!,:plantRRL][1])/2
				Ylb[:,vs] .= plant[!,:pplantlb][1] - plant[!,:plantRRL][1]/2
				Yub[:,vs] .= plant[!,:pplantub][1] + plant[!,:plantRRL][1]/2
	
			elseif vs > numP && vs <= 2*numP
				# YZ init = (demand+flex)/24
				# startv[:,vs] .= plant[!,:plantdemand][1]/24
				startv[:,vs] .= 0.0
				# Ylb[:,vs] .= -0.1
				# Yub[:,vs] .= 0.1
				# startv = (1+plant[!,:plantflex][1]).*plant[!,:plantdemand][1]/24
				Ylb[:,vs] .= -(1+plant[!,:plantflex][1]).*plant[!,:plantdemand][1]/24
				Yub[:,vs] .= (1+plant[!,:plantflex][1]).*plant[!,:plantdemand][1]/24
			else
				# YV
				startv[:,vs] .= 1*G[1][!,:Pub][vs-numP*2]/2
				Ylb[:,vs] .= G[1][!,:Plb][vs-numP*2] - 1*G[1][!,:Pub][vs-numP*2]/2
				Yub[:,vs] .= 1*G[1][!,:Pub][vs-numP*2]/2 + G[1][!,:Pub][vs-numP*2]
			end

		end

		# println("start ", startv[:,1:(2*numP+3)], " Ylb ", Ylb[:,1:(2*numP+3)], " Yub ", Yub[:,1:(2*numP+3)])
		return startv, Ylb, Yub
	end

	# initialize Yk
	Yk = Yint()[1]
	# initialize to save value of obj
	objvalue = 0.0

	# iteration index and tolerance 
	iteridx = 0
	tol = 1

	# SQP loop 
	while tol > 1e-5 && iteridx < 45
		println("----------------iter ", iteridx, "---------------------")
		flush(stdout)
		# create model
		m = Model(SQPNLSolver)

		# Master variables Y
		Ylb, Yub = Yint()[2:3]
		@variable(m, Ylb[ts,vs] <= Ylink[ts=1:tsend, vs=1:(2*numP+numG)] <= Yub[ts,vs], start = Yk[ts,vs])

		# println("start values ", start_value.(Ylink[:,1:(2*numP+3)]))
		# flush(stdout)

		# sum_t(Z) = 0 
		for vs = 1:numP
			@constraint(m, sum(Ylink[ts, vs+numP] for ts=timestep) == 0)
		end
		# Ut and Zt: demand/(T*unitProd) - RRL/2 <= Ut - Zt/unitProd <= (demand+flex)/(T*unitProd) + RRL/2
		for ts = timestep
			@constraint(m, [vs=1:numP], Ylink[ts,vs] - Ylink[ts,numP+vs]/plant[!,:plantprod][vs] <= plant[!,:plantRRL][vs]/2 + (1+plant[!,:plantflex][vs])*plant[!,:plantdemand][vs]/24/plant[!,:plantprod][vs])
			@constraint(m, [vs=1:numP], -plant[!,:plantRRL][vs]/2 + plant[!,:plantdemand][vs]/24/plant[!,:plantprod][vs] <= Ylink[ts,vs] - Ylink[ts,numP+vs]/plant[!,:plantprod][vs])
		end
		# Ut and Zt+1: same relationship
		for ts = timestep[1:(end-1)]
			@constraint(m, [vs=1:numP], Ylink[ts,vs] - Ylink[ts+1,numP+vs]/plant[!,:plantprod][vs] <= plant[!,:plantRRL][vs]/2 + (1+plant[!,:plantflex][vs])*plant[!,:plantdemand][vs]/24/plant[!,:plantprod][vs])
			@constraint(m, [vs=1:numP], -plant[!,:plantRRL][vs]/2 + plant[!,:plantdemand][vs]/24/plant[!,:plantprod][vs] <= Ylink[ts,vs] - Ylink[ts+1,numP+vs]/plant[!,:plantprod][vs])
		end
		# -RRL <= Ut - Ut-1 <= RRL
		for ts = timestep[1:(end-1)]
			@constraint(m, [vs=1:numP], -plant[!,:plantRRL][vs] <= Ylink[ts+1,vs] - Ylink[ts,vs])
			@constraint(m, [vs=1:numP], Ylink[ts+1,vs] - Ylink[ts,vs] <= plant[!,:plantRRL][vs])
		end
	
		# -RRL_gen <= Vt - Vt-1 <= RRL_gen
		for ts = timestep[1:(end-1)]
			@constraint(m, [vs=1:numG], -1*G[1][!,:Pub][vs] <= Ylink[ts+1,vs+2*numP] - Ylink[ts,vs+2*numP])
			@constraint(m, [vs=1:numG], Ylink[ts+1,vs+2*numP] - Ylink[ts,vs+2*numP] <= 1*G[1][!,:Pub][vs])
		end

		# initialize objk at the beginning of each iteration
		subobjk = 0 
		objreg = zero(JuMP.QuadExpr)
		# get current values for Ylink at ts 
		for ts = timestep
			if ts == 1
				Ylinkts = Ylinkdecomts(Yk[ts,1:numP], Yk[ts,(numP+1):numP*2], Yk[ts,(numP*2+1):end], Yk[ts,1:numP], Yk[ts,(numP*2+1):end], ts)
			else
				Ylinkts = Ylinkdecomts(Yk[ts,1:numP], Yk[ts,(numP+1):numP*2], Yk[ts,(numP*2+1):end], Yk[ts-1,1:numP], Yk[ts-1,(numP*2+1):end], ts)
			end
			println("---------------subproblem ts ", ts, "---------------------")
			flush(stdout)
			# call to solve subproblems
			objts, gts = SolveSubr(o, N[ts], L[ts], T[ts], SSh[ts], G[ts], K[ts], P[ts], plant, numP, MVAbase, Ylinkts, gkm1[ts,:], timestep, NLSolver)
			# add all obj from subproblems
			subobjk += objts
			# save duals from subproblems
			gkm1[ts,:] = gts
			# store dDt/dYt
			gkts[ts,:] .= gts[1:(2*numP+numG)]
			# store dDt/dUtm1 and dVt/dUtm1
			gktsm[ts,:] .= gts[(2*numP+numG+1):end]
		end

		# gUkt = dDt/dUt + dDt+1/dUt, gVkt = dDt/dVt + dDt+1/dVt, i.e., dD/dU1 = dD2/dU1 + dD1/dU1
		# gkmts[ts+1] = dDt+1/dYt, gkts[ts] = dDt/dYt
		for ts = timestep[1:(end-1)]
			gkts[ts,1:numP] += gktsm[(ts+1),1:numP]
			gkts[ts,(numP*2+1):end] += gktsm[(ts+1),(numP+1):end]
		end
		# println("gk[1:2,1:8] ", gkts[1:2,1:8])
		# println("gkm1[1:2,1:8] ", gkm1[1:2,1:8])
		# println("gktsm[1:2,1:8] ", gktsm[1:2,1:8])
		flush(stdout)
		# subobjk *= 1/1024
		# gkts .*= 1/1024
		# alphak *= 1/1024
		# add SQP gradient and hessian terms to objreg
		for ts = timestep
			# g'*(Y-Yk) & 0.5*(Y-Yk)'*Bk*(Y-Yk)
			grad = @expression(m, gkts[ts,:]'*(Ylink[ts,:]-Yk[ts,:]))
			add_to_expression!(objreg, grad)
			hess = @expression(m, 0.5*(Ylink[ts,:]-Yk[ts,:])'*alphak*I*(Ylink[ts,:]-Yk[ts,:]))
			add_to_expression!(objreg, hess)
		end

		# declare objective
		@NLobjective(m, Min, subobjk + objreg)
		## attempt to solve SCACOPF
		JuMP.optimize!(m)
 
        # termination_status, primal_status, and dual_status: http://www.juliaopt.org/JuMP.jl/dev/solutions/
        solvestatusprimal = JuMP.primal_status(m)
        solvestatus = JuMP.termination_status(m)
        println("JuMP optimize MOI termination_status: ", solvestatus, " ", solvestatusprimal)
		if (solvestatusprimal != MOI.FEASIBLE_POINT) && (solvestatusprimal != MOI.NEARLY_FEASIBLE_POINT)
			error("solver failed to find a feasible solution")
		end

		## output Y for next iteration
		dYk = JuMP.value.(Ylink)
		dk = dYk - Yk
		# save objective value fk and reg term
		objvalue = subobjk
		println("Obj value \$", round(objvalue, digits=2))
		# predicted = gk*dk'+0.5*alpha*I*||dk||^2
		# regvalue = sum(gkts*dk')+0.5*alphak*norm(dk)
		regvalue = JuMP.value(objreg)
		println("predicted change ", regvalue)

		# update Yk to Yk+betak*dk
		Yk .= Yk + betak*dk
		# initialize new obj
		objdk = 0
		# check value of f(Yk+betak*dk)
		for ts = timestep
			if ts == 1
				Ylinkts = Ylinkdecomts(Yk[ts,1:numP], Yk[ts,(numP+1):numP*2], Yk[ts,(numP*2+1):end], Yk[ts,1:numP], Yk[ts,(numP*2+1):end], ts)
			else
				Ylinkts = Ylinkdecomts(Yk[ts,1:numP], Yk[ts,(numP+1):numP*2], Yk[ts,(numP*2+1):end], Yk[ts-1,1:numP], Yk[ts-1,(numP*2+1):end], ts)
			end
			println("------------alpha check: subproblem ts ", ts, "---------------------")
			flush(stdout)
			# call to solve subproblems
			objts = SolveSubr(o, N[ts], L[ts], T[ts], SSh[ts], G[ts], K[ts], P[ts], plant, numP, MVAbase, Ylinkts, gkm1[ts,:], timestep, NLSolver)[1]
			# add all obj from subproblems
			objdk += objts
		end
		# objdk *= 1/1024
		Dobj = objdk - objvalue
		println("true change ", Dobj)

		# compute ratio of true change and predicted change
		TPratio = Dobj/regvalue
		# check if current alpha is acceptable
		if TPratio < 0.1 
			# if not, keep Yk=Yk, increase alpha
			println("don't update Yk, change alpha")
			Yk .= Yk - betak*dk
			# alphak *= 10
			if alphak <= 1e6
				alphak *= 10
			elseif alphak < 2e7
				alphak *= 2
				else
					alphak += 100
			end
		end

		tol = norm(dk)
		
		# gradsqr = sum((gkts).^2)
		println("alpha ", alphak)
		println("sum of (Y - Yk).^2 ", tol)
		# tol = norm(dk)*alphak
		# println("sum of (Y - Yk).^2 *alpha ", tol)
		flush(stdout)
		# println("sum of gk.^2 ", gradsqr)
		iteridx += 1

	end
	println("Master Objective:           \$", round(objvalue, digits=2), "\n")
	flush(stdout)
	
end

## function for solving MP AC OPF
function SolveMPACOPF(o::GoOptions, N::Vector{Any}, L::Vector{Any}, T::Vector{Any},
	              SSh::Vector{Any}, G::Vector{Any}, K::Vector{Any}, P::Vector{Any}, timestep::Vector{Int64}, plant::DataFrame, MVAbase::Float64, NLSolver; 
                      IndexSets::Union{Tuple, Nothing}=nothing,
                      StartingPoint::Union{SCACOPFState,Nothing}=nothing)

	# compute index and sets for performant formulation
        L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn, K_outidx =
            build_or_split_indexsets(IndexSets, N[timestep[1]], L[timestep[1]], T[timestep[1]], SSh[timestep[1]], G[timestep[1]], K[timestep[1]])

	RefBus = G_Nidx[argmax(G[timestep[1]][!,:Pub])]	# bus with the biggest generator is reference
	# println("RefBus G_Nidx ", argmax(G[timestep[1]][!,:Pub]))
	NumK = size(K[timestep[1]], 1)

        # (some) output for options used
        println("SolveSCACOPF with ", NumK, " contingencies. \nOptions: CouplingMode=", o.CouplingMode, 
                " Smoothing=", o.SmoothingParamCoupling,
		" ThermalLimits=", o.ThermalLimits)
		flush(stdout)

        if StartingPoint==nothing
	    # compute SC AC OPF starting point 
	    v0_n, theta0_n, p0_li, q0_li, p0_ti, q0_ti, b0_s, p0_g, q0_g, 
            pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti =
                InitialBaseCaseSolution( N[timestep[1]], L[timestep[1]], T[timestep[1]], SSh[timestep[1]], G[timestep[1]], 
                                         IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
		                                      SShn, Gn, K_outidx))
	    theta0_n .-= theta0_n[RefBus]
        else
            #! this is also temporary code
            @assert num_cont(StartingPoint)==NumK

            ## base case
            v0_n, theta0_n, b0_s, p0_g, q0_g = getBaseStates(StartingPoint)
            
            p0_li, q0_li, p0_ti, q0_ti, pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti = 
                BaseCaseSolution(v0_n, theta0_n, b0_s, p0_g, q0_g,
                                 N[timestep[1]], L[timestep[1]], T[timestep[1]], SSh[timestep[1]], G[timestep[1]], 
                                 IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
		                              SShn, Gn, K_outidx))
            println("HOT Started")
        end

	# create model
	m = Model(NLSolver)
	
	# MPACOPF variables
	@variable(m, N[ts][!,:Vlb][n] <= v_n[ts=timestep, n=1:size(N[timestep[1]], 1)] <= N[ts][!,:Vub][n], start = v0_n[n])
	@variable(m, theta_n[ts=timestep, n=1:size(N[timestep[1]], 1)], start = theta0_n[n])
	@variable(m, p_li[ts=timestep, l=1:size(L[timestep[1]], 1), i=1:2], start = p0_li[l,i])
	@variable(m, q_li[ts=timestep, l=1:size(L[timestep[1]], 1), i=1:2], start = q0_li[l,i])
	@variable(m, p_ti[ts=timestep, t=1:size(T[timestep[1]], 1), i=1:2], start = p0_ti[t,i])
	@variable(m, q_ti[ts=timestep, t=1:size(T[timestep[1]], 1), i=1:2], start = q0_ti[t,i])
	@variable(m, SSh[ts][!,:Blb][s] <= b_s[ts=timestep, s=1:size(SSh[timestep[1]], 1)] <= SSh[ts][!,:Bub][s], start = b0_s[s])
	@variable(m, G[ts][!,:Plb][g] <= p_g[ts=timestep, g=1:size(G[timestep[1]], 1)] <= G[ts][!,:Pub][g], start = p0_g[g])
	@variable(m, G[ts][!,:Qlb][g] <= q_g[ts=timestep, g=1:size(G[timestep[1]], 1)] <= G[ts][!,:Qub][g], start = q0_g[g])
	@variable(m, pslackm_n[ts=timestep, n=1:size(N[timestep[1]], 1)] >= 0, start = pslack0m_n[n])
	@variable(m, pslackp_n[ts=timestep, n=1:size(N[timestep[1]], 1)] >= 0, start = pslack0p_n[n])
	@variable(m, qslackm_n[ts=timestep, n=1:size(N[timestep[1]], 1)] >= 0, start = qslack0m_n[n])
	@variable(m, qslackp_n[ts=timestep, n=1:size(N[timestep[1]], 1)] >= 0, start = qslack0p_n[n])
	@variable(m, sslack_li[ts=timestep, l=1:size(L[timestep[1]], 1), i=1:2] >= 0, start = sslack0_li[l,i])
	@variable(m, sslack_ti[ts=timestep, t=1:size(T[timestep[1]], 1), i=1:2] >= 0, start = sslack0_ti[t,i])

	# # check generation cost 
	# t_h01 = [0.2, 0.8]
	# t_h1 = @variable(m, [ts=timestep, i=1:2], lower_bound=0.0, start=t_h01[i])
	# prodcost1 = zero(JuMP.AffExpr)
	# prodcost2 = zero(JuMP.AffExpr)
	# prodcost3 = zero(JuMP.AffExpr)
	# prodcost4 = zero(JuMP.AffExpr)
	# prodcost5 = zero(JuMP.AffExpr)
	# # costh1 = @variable(m, [ts=timestep], start = 0.0)
	# gh1 = @variable(m, [ts=timestep], start = 0.0)

	# check # of plants
	nplant = findall( x -> x in plant[!,:plantid], N[1][!, :Bus])
	if size(nplant,1) != 0
		plant[!,:busidx] .= nplant
	end

	# bidding entity variables at timestep
	@variable(m, plant[!,:pplantlb][np] <= pplant[ts=timestep, np=1:size(nplant,1)] <= plant[!,:pplantub][np])
	@variable(m, plprod[ts=timestep, np=1:size(nplant,1)])
	@variable(m, plprodtot[np=1:size(nplant,1)])
	
	# bidding plant constraints
	# calculate plant prod from energy consumption for each timestep at corresponding bus
	@constraint(m, [ts=timestep, np=1:size(nplant,1)], plprod[ts,np] == plant[!,:plantprod][np].*pplant[ts,np])
	# ramping constraint for plant production
	@NLconstraint(m, [ts=timestep[2:end], np=1:size(nplant,1)], -plant[!,:plantRRL][np] <= pplant[ts,np] - pplant[ts-1,np] <= plant[!,:plantRRL][np])
	# temporal constraint: plant production within demand and flex*demand 
	@constraint(m, [np=1:size(nplant,1)], plprodtot[np] == sum(plprod[ts,np] for ts=timestep))
	@constraint(m, [np=1:size(nplant,1)], plant[!,:plantdemand][np]*timestep[end]/24 <= plprodtot[np] <= (1+plant[!,:plantflex][np]).*plant[!,:plantdemand][np]*timestep[end]/24)

	#active power balance with LMP=dual
	@constraint(m, PBPC[ts = timestep, np=1:size(nplant,1)], sum(p_g[ts,g] for g=Gn[plant[!,:busidx][np]]) - N[ts][!,:Pd][plant[!,:busidx][np]] - N[ts][!,:Gsh][plant[!,:busidx][np]]*v_n[ts,plant[!,:busidx][np]]^2 -
			pplant[ts,np] - 
			sum(p_li[ts,Lidxn[plant[!,:busidx][np]][lix],Lin[plant[!,:busidx][np]][lix]] for lix=1:length(Lidxn[plant[!,:busidx][np]])) -
			sum(p_ti[ts,Tidxn[plant[!,:busidx][np]][tix],Tin[plant[!,:busidx][np]][tix]] for tix=1:length(Tidxn[plant[!,:busidx][np]])) ==
			pslackp_n[ts,plant[!,:busidx][np]] - pslackm_n[ts,plant[!,:busidx][np]])
	
	nbusn = findall( x -> !(x in plant[!,:plantid]), N[1][!, :Bus])
	# println("size of nbusn ", size(nbusn,1), " nbusn[1:10] ", nbusn[1:10])

	@constraint(m, PB[ts=timestep, n=nbusn], sum(p_g[ts,g] for g=Gn[n]) - N[ts][!,:Pd][n] - N[ts][!,:Gsh][n]*v_n[ts,n]^2 -
			sum(p_li[ts,Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n])) -
			sum(p_ti[ts,Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n])) ==
			pslackp_n[ts,n] - pslackm_n[ts,n])


	
	
	# constrain total power generation
	GenTot =  1.4747684463111619e6
	# @constraint(m, sum(p_g[ts,g] for ts=timestep, g=1:size(G[ts], 1)) <= GenTot/MVAbase)

	# Define vectors to save the indexes of inequalities for line and transformer
	Lsta = []
	Tsta = []

	# power flow constraints
	for ts = timestep
		# fix angle at reference bus to zero
		JuMP.fix(theta_n[ts,RefBus], 0.0, force=true)
		
		# add power flow constraints for base case
		addpowerflowcons!(m, v_n[ts,:], theta_n[ts,:], p_li[ts,:,:], q_li[ts,:,:], p_ti[ts,:,:], q_ti[ts,:,:], b_s[ts,:], p_g[ts,:], q_g[ts,:],
			pslackm_n[ts,:], pslackp_n[ts,:], qslackm_n[ts,:], qslackp_n[ts,:], sslack_li[ts,:,:], sslack_ti[ts,:,:],
		 	N[ts], L[ts], T[ts], (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
			#  SShn, Gn, K_outidx), o.ThermalLimits)
			SShn, Gn, K_outidx), o.ThermalLimits, ts, Lsta, Tsta)
	
		# register approximation functions
		sclamp = nothing
		sstep = nothing
	end
	

	println("start ramping constraints")
	flush(stdout)
	for ts = timestep[2:end]
		#ramping constraints
		@NLconstraint(m, [g=1:size(G[ts], 1)], -1*G[ts][!,:Pub][g] <= p_g[ts,g]-p_g[ts-1,g] <= 1*G[ts][!,:Pub][g])
	end


## objective function
prodcost = zero(JuMP.AffExpr)
prodcoststep = zeros(JuMP.AffExpr, timestep[end])
basepenalty = zero(JuMP.QuadExpr)
plantrev = zero(JuMP.AffExpr)

println("start prod cost")
flush(stdout)
for ts = timestep

	# generator production cost
	for g = 1:size(G[ts], 1)
		npairs = length(G[ts][!,:CostPi][g])
			@assert sort(G[ts][!,:CostPi][g]) == G[ts][!,:CostPi][g]

			t_h0 = getstartcostcoefficients(JuMP.start_value(p_g[ts,g]), G[ts][!,:CostPi][g], G[ts][!,:Plb][g], G[ts][!,:Pub][g])
			t_h = @variable(m, [i=1:npairs], lower_bound=0.0, start=t_h0[i])
		
		 @constraint(m, sum(t_h[h] for h=1:npairs) == 1.0)
		@constraint(m, p_g[ts,g] == sum(G[ts][!,:CostPi][g][h]*t_h[h] for h=1:npairs))
		gencost = @expression(m, sum(G[ts][!,:CostCi][g][h]*t_h[h] for h=1:npairs))
		add_to_expression!(prodcost, gencost)
		prodcoststep[ts] += gencost
	end
end

	println("start base case pen")
	flush(stdout)
	for ts = timestep
		# base case penalties
		basecasepenalty = quadpenaltyfunctions(pslackm_n[ts,:], pslackp_n[ts,:], qslackm_n[ts,:], qslackp_n[ts,:],
		sslack_li[ts,:,:], sslack_ti[ts,:,:], P[ts])
		add_to_expression!(basepenalty, basecasepenalty)
		
	end
	
	#plant prod revenue
	println("start plant rev")
	flush(stdout)
	for np=1:size(nplant,1)
		add_to_expression!(plantrev, sum(pplant[ts,np].*plant[!,:plantrev][np] for ts=timestep))
	end
	


	# declare objective
	# @objective(m, Min, (prodcost + DELTA*basepenalty))
	@objective(m, Min, (prodcost + DELTA*basepenalty - plantrev))
	
	## attempt to solve SCACOPF
	JuMP.optimize!(m)
 
        # termination_status, primal_status, and dual_status: http://www.juliaopt.org/JuMP.jl/dev/solutions/
        solvestatusprimal = JuMP.primal_status(m)
        solvestatus = JuMP.termination_status(m)
        println("JuMP optimize MOI termination_status: ", solvestatus, " ", solvestatusprimal)
	if (solvestatusprimal != MOI.FEASIBLE_POINT) && (solvestatusprimal != MOI.NEARLY_FEASIBLE_POINT)
		error("solver failed to find a feasible solution")
	end
	
	## objective breakdown
	#prodcostvalue = JuMP.value(prodcost)
	penaltybasevalue = DELTA*JuMP.value(basepenalty)
	plantrevvalue = JuMP.value(plantrev)
	objectivevalue = objective_value(m)
	# prodcostvalue = objectivevalue - penaltybasevalue 
	prodcostvalue = objectivevalue - penaltybasevalue + plantrevvalue

	
	println("Objective:           \$", round(objectivevalue, digits=1), "\n",
		"Generation cost:     \$", round(prodcostvalue, digits=1), "\n",
		"Penalty base:        \$", round(penaltybasevalue, digits=1), "\n",
		"Plant revenue:        \$", round(plantrevvalue, digits=1), "\n")

	flush(stdout)

	
	## outputs
	plantpower = JuMP.value.(pplant).*MVAbase
	# println("plant power:      ", plantpower, "\n")

	
	LMPs = Containers.DenseAxisArray{Float64}(undef, timestep, 1:size(nplant,1))
	for ts = timestep
		for n = 1:size(nplant,1)
		# for n = nplant
			LMPs[ts,n] = dual(PBPC[ts,n])/MVAbase
		end
	end
	# println("Plant LMPs   ", LMPs, "\n")

	@assert penaltybasevalue <= 0
	
	LMPsN = Containers.DenseAxisArray{Float64}(undef, timestep, 1:size(nbusn,1))
	for ts = timestep
		for n = eachindex(nbusn)
			LMPsN[ts,n] = dual(PB[ts,nbusn[n]])/MVAbase
		end
	end
	println("First bus LMPs (no generator)   ", LMPsN[:,1], "\n")
	println("Bus 1004 LMPs (with generator)    ", LMPsN[:,4], "\n")
	# println("Plant LMPs   ", LMPsN[:,nplant[1]], "\n")

	p_gv = JuMP.value.(p_g)
    Gent = zeros(size(timestep,1))
    Gentot = 0.0 
    for ts = timestep
        Gent[ts] = sum(p_gv[ts,g] for g=1:size(G[ts], 1))
        Gentot += Gent[ts]
    end
    Gentot *= MVAbase
    println("Total generation:  ", Gentot)

	# get values of line losses from line and transformer powers, and power generation sum
	p_liv = JuMP.value.(p_li)
	# println(typeof(p_liv), size(p_liv))
	p_tiv = JuMP.value.(p_ti)
	p_gv = JuMP.value.(p_g)
	Lloss = Containers.DenseAxisArray{Float64}(undef, timestep, 1:size(N[timestep[1]], 1))
	Genb = Containers.DenseAxisArray{Float64}(undef, timestep, 1:size(N[timestep[1]], 1))
	Llosst = zeros(size(timestep,1))
	Gent = zeros(size(timestep,1))
	eta = zeros(size(timestep,1))
	for ts = timestep
		for n = 1:size(N[timestep[1]], 1)
			Genb[ts,n] = sum(p_gv[ts,g] for g=Gn[n]; init=0)
			Lloss[ts,n] = sum(p_liv[ts,Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n]); init=0) +
			sum(p_tiv[ts,Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n]); init=0)
		end
		Llosst[ts] = sum(Lloss[ts,ns] for ns = 1:size(N[timestep[1]], 1); init=0)*MVAbase
		Gent[ts] = sum(Genb[ts,ns] for ns = 1:size(N[timestep[1]], 1); init=0)*MVAbase
		eta[ts] = Llosst[ts]/Gent[ts]
	end
	# println("Percent line loss[1:3] ", eta[1:3])

	# Define empty DataFrames to store constraint information for line and transformer inequalities
	LieqnA = DataFrame(time = Int[], line = Int[], phase = Int[], status = Int[])
	TieqnA = DataFrame(time = Int[], trans = Int[], phase = Int[], status = Int[])

	# calculate constraint values, save status 
	for (ts,l,i) in Lsta
		Lconstv = value(p_li[ts,l,i]^2 + q_li[ts,l,i]^2 - 
		((L[ts][!, :RateBase][l]*v_n[ts,L_Nidx[l,i]]) + sslack_li[ts,l,i])^2)
		Lstatus = isapprox(Lconstv, 0, atol=1e-6) ? 1 : 0  # Change status to 1 for active, 0 for inactive
		if Lstatus == 1
			Lstatus *= ts # relate active status to timestep
			push!(LieqnA, (time=ts, line=l, phase=i, status=Lstatus))
		end
	end
	# println("size of Lieqn: ", size(Lieqn))
	println("active line constraints: ", LieqnA)

	for (ts,t,i) in Tsta
		Tconstv = value(p_ti[ts,t,i]^2 + q_ti[ts,t,i]^2 -
		(T[ts][!, :RateBase][t]+ sslack_ti[ts,t,i])^2)
		Tstatus = isapprox(Tconstv, 0, atol=1e-6) ? 1 : 0  # Change status to 1 for active, 0 for inactive
		if Tstatus == 1
			Tstatus *= ts # relate active status to timestep
			push!(TieqnA, (time=ts, trans=t, phase=i, status=Tstatus))
		end
	end
	# println("size of Tieqn: ", size(Tieqn))
	println("active transformer constraints: ", TieqnA)
	
	GenCostStep = JuMP.value.(prodcoststep)
	println("GenCost @ timestep ", GenCostStep)

	# get values of slack variables
	slackpm = JuMP.value.(pslackm_n)
	slackpp = JuMP.value.(pslackp_n)
	slackqm = JuMP.value.(qslackm_n)
	slackqp = JuMP.value.(qslackp_n)
	slackl = JuMP.value.(sslack_li)
	slackt = JuMP.value.(sslack_ti)


        if false
            kk=2
            vp_k = JuMP.value.(p_gk)
            vp_k = view(vp_k, :, kk)
            vp   = JuMP.value.(p_g)
            vPlb = G[!,:Plb]
            vPub = G[!,:Pub]
            delta= JuMP.value.(delta_k)
            delta= delta[kk]
            alpha = G[!,:alpha]
            for i=1:length(vp_k)
                @printf("g=%4d pg=%12.5e  pg_k=%12.5e diffPlb=%12.5e diffPub=%12.5e  p_k-(p+alpha*delta)=%12.5e fPlb=%12.5e Pub=%12.5e alpha=%12.5e  delta=%12.5e\n", 
                        i, vp[i], vp_k[i], vp_k[i]-vPlb[i], vPub[i]-vp_k[i], vp_k[i]-(vp[i]+alpha[i]*delta), vPlb[i], vPub[i], alpha[i], delta)
                #@printf("alpha[i]=%g  delta=%g\n", alpha[i], delta)
            end
            #println("G Plb", G[!,:Plb])
            #println("G Pub", G[!,:Pub])
        end

        if false
            vq_gk = JuMP.value.(q_gk)
            vq_gk = view(vq_gk,:,1)
            vq_g  = JuMP.value.(q_g)
            vv_nk = JuMP.value.(v_nk)
            vv_nk = view(vv_nk,:,1)
            vv_n  = JuMP.value.(v_n)
            vQlb = G[!,:Qlb]
            vQub = G[!,:Qub]
            
            for i=1:length(vq_gk)
                @printf("PVPQ gen=%3d q_g=%12.5e q_gk=%12.5e v_n=%12.5e v_nk=%12.5e | v-v_k=%12.5e Qub-q_k=%12.5e q_k-Qlb=%12.5e | Qlb=%12.5e Qub=%12.5e\n",
                        i, vq_g[i], vq_gk[i], vv_n[i], vv_nk[i], vv_n[i]- vv_nk[i],  vQub[i]-vq_gk[i],  vq_gk[i]-vQlb[i], vQlb[i], vQub[i])
            end
                    
        end

        if false
            #println(JuMP.value.(v_n))
            println("p_g", JuMP.value.(p_g))
            println("p_")
            println("G Plb", G[!,:Plb])
            println("G Pub", G[!,:Pub])


            println("q_g", JuMP.value.(q_g))
            println("G Qlb", G[!,:Qlb])
            println("G Qub", G[!,:Qub])

            #println("N[!,:PD]", N[!,:Pd])
            #println("N[!,:QD]", N[!,:Qd])
            
            qslackm_nv = JuMP.value.(qslackm_n)
            qslackp_nv = JuMP.value.(qslackp_n)
	    #@variable(m, pslackm_n[n=1:size(N, 1)] >= 0)
            println("pslackm_n", JuMP.value.(pslackm_n))
	    #@variable(m, pslackp_n[n=1:size(N, 1)] >= 0)
            println("pslackp_n", JuMP.value.(pslackp_n))
	    #@variable(m, qslackm_n[n=1:size(N, 1)] >= 0)
            println("qslackm_n", JuMP.value.(qslackm_n))
	    #@variable(m, qslackp_n[n=1:size(N, 1)] >= 0)
            println("qslackm_n", JuMP.value.(qslackp_n))
	    #@variable(m, sslack_li[l=1:size(L, 1), i=1:2] >= 0)
	    #@variable(m, sslack_ti[t=1:size(T, 1), i=1:2] >= 0)
            p_liv = JuMP.value.(p_li)
            p_tiv = JuMP.value.(p_ti)
            q_liv = JuMP.value.(q_li)
            q_tiv = JuMP.value.(q_ti)
            q_gv = JuMP.value.(q_g)
            p_gv = JuMP.value.(p_g)
            b_sv = JuMP.value.(b_s)
            println("b_sv: ", b_sv)
	    for n = 1:size(N, 1)
                #		@NLconstraint(m, sum(q_g[g] for g=Gn[n]) - N[!,:Qd][n] -
		#	(N[!,:Bsh][n] + sum(b_s[s] for s=SShn[n]))*v_n[n]^2 -
		#	sum(q_li[Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n])) -
		#	sum(q_ti[Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n])) ==
		#	qslackp_n[n] - qslackm_n[n])
                sumpg = 0.
                if length(Gn[n])>0
                    sumpg = sum(q_gv[g] for g=Gn[n])
                end
                sumqli=0.;
                if length(Lidxn[n])>0
                    sumqli = sum(q_liv[Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n]))
                end
                sumqti=0.
                if length(Tidxn[n])>0
                    sumqti = sum(q_tiv[Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n]))
                end
		NBsh = N[!,:Bsh][n]
                if length(SShn[n])>0
                    NBsh = NBsh +  sum(b_sv[s] for s=SShn[n])
                end
                if qslackm_nv[n]>1e-6 || qslackp_nv[n]>1e-6
                    @printf("n=%d  sumpg=%12.5e  N[!,:Qd]=%12.5e  N[!,:Bsh]=%12.5e sum_qli=%12.5e sum_qti=%12.5e slackm=%12.5e  slackp=%12.5e\n", n, sumpg, N[!,:Qd][n], NBsh, sumqli, sumqti, qslackm_nv[n], qslackp_nv[n])
                end
            end
        end
	
	## return solution
	return MPACOPFState(MPBaseState(JuMP.value.(v_n), 
                                	JuMP.value.(theta_n), 
                                    JuMP.value.(b_s),
                            JuMP.value.(p_g), 
                                      JuMP.value.(q_g)),
									  biddingRes(prodcostvalue, penaltybasevalue, plantrevvalue, plantpower,LMPs,LMPsN),
									  slackv(slackpm, slackpp, slackqm, slackqp, slackl, slackt),
									  TranLoss(Llosst, Gent, eta, p_liv, p_tiv, Genb, Lloss, GenCostStep, LieqnA, TieqnA))
	
end

## function for solving MP AC OPF
function SolveMPACOPFsingle(o::GoOptions, N::Vector{Any}, L::Vector{Any}, T::Vector{Any},
	              SSh::Vector{Any}, G::Vector{Any}, K::Vector{Any}, P::Vector{Any}, timestep::Vector{Int64}, plant::DataFrame, MVAbase::Float64, NLSolver; 
                      IndexSets::Union{Tuple, Nothing}=nothing,
                      StartingPoint::Union{SCACOPFState,Nothing}=nothing)

	# compute index and sets for performant formulation
        L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn, K_outidx =
            build_or_split_indexsets(IndexSets, N[timestep[1]], L[timestep[1]], T[timestep[1]], SSh[timestep[1]], G[timestep[1]], K[timestep[1]])

	RefBus = G_Nidx[argmax(G[timestep[1]][!,:Pub])]	# bus with the biggest generator is reference
	NumK = size(K[timestep[1]], 1)

        # (some) output for options used
        println("SolveSCACOPF with ", NumK, " contingencies. \nOptions: CouplingMode=", o.CouplingMode, 
                " Smoothing=", o.SmoothingParamCoupling,
		" ThermalLimits=", o.ThermalLimits)
		flush(stdout)

        if StartingPoint==nothing
	    # compute SC AC OPF starting point 
	    v0_n, theta0_n, p0_li, q0_li, p0_ti, q0_ti, b0_s, p0_g, q0_g, 
            pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti =
                InitialBaseCaseSolution( N[timestep[1]], L[timestep[1]], T[timestep[1]], SSh[timestep[1]], G[timestep[1]], 
                                         IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
		                                      SShn, Gn, K_outidx))
	    theta0_n .-= theta0_n[RefBus]
        else
            #! this is also temporary code
            @assert num_cont(StartingPoint)==NumK

            ## base case
            v0_n, theta0_n, b0_s, p0_g, q0_g = getBaseStates(StartingPoint)
            
            p0_li, q0_li, p0_ti, q0_ti, pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti = 
                BaseCaseSolution(v0_n, theta0_n, b0_s, p0_g, q0_g,
                                 N[timestep[1]], L[timestep[1]], T[timestep[1]], SSh[timestep[1]], G[timestep[1]], 
                                 IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
		                              SShn, Gn, K_outidx))
            println("HOT Started")
        end

	# create model
	m = Model(NLSolver)
	
	# MPACOPF variables
	@variable(m, N[ts][!,:Vlb][n] <= v_n[ts=timestep, n=1:size(N[timestep[1]], 1)] <= N[ts][!,:Vub][n], start = v0_n[n])
	@variable(m, theta_n[ts=timestep, n=1:size(N[timestep[1]], 1)], start = theta0_n[n])
	@variable(m, p_li[ts=timestep, l=1:size(L[timestep[1]], 1), i=1:2], start = p0_li[l,i])
	@variable(m, q_li[ts=timestep, l=1:size(L[timestep[1]], 1), i=1:2], start = q0_li[l,i])
	@variable(m, p_ti[ts=timestep, t=1:size(T[timestep[1]], 1), i=1:2], start = p0_ti[t,i])
	@variable(m, q_ti[ts=timestep, t=1:size(T[timestep[1]], 1), i=1:2], start = q0_ti[t,i])
	@variable(m, SSh[ts][!,:Blb][s] <= b_s[ts=timestep, s=1:size(SSh[timestep[1]], 1)] <= SSh[ts][!,:Bub][s], start = b0_s[s])
	@variable(m, G[ts][!,:Plb][g] <= p_g[ts=timestep, g=1:size(G[timestep[1]], 1)] <= G[ts][!,:Pub][g], start = p0_g[g])
	@variable(m, G[ts][!,:Qlb][g] <= q_g[ts=timestep, g=1:size(G[timestep[1]], 1)] <= G[ts][!,:Qub][g], start = q0_g[g])
	@variable(m, pslackm_n[ts=timestep, n=1:size(N[timestep[1]], 1)] >= 0, start = pslack0m_n[n])
	@variable(m, pslackp_n[ts=timestep, n=1:size(N[timestep[1]], 1)] >= 0, start = pslack0p_n[n])
	@variable(m, qslackm_n[ts=timestep, n=1:size(N[timestep[1]], 1)] >= 0, start = qslack0m_n[n])
	@variable(m, qslackp_n[ts=timestep, n=1:size(N[timestep[1]], 1)] >= 0, start = qslack0p_n[n])
	@variable(m, sslack_li[ts=timestep, l=1:size(L[timestep[1]], 1), i=1:2] >= 0, start = sslack0_li[l,i])
	@variable(m, sslack_ti[ts=timestep, t=1:size(T[timestep[1]], 1), i=1:2] >= 0, start = sslack0_ti[t,i])

	# check # of plants
	nplant = findall( x -> x in plant[!,:plantid], N[1][!, :Bus])
	if size(nplant,1) != 0
		plant[!,:busidx] .= nplant
	end

	# bidding entity variables at timestep
	@variable(m, plant[!,:pplantlb][np] <= pplant[ts=timestep, np=1:size(nplant,1)] <= plant[!,:pplantub][np])
	@variable(m, plprod[ts=timestep, np=1:size(nplant,1)])
	@variable(m, plprodtot[np=1:size(nplant,1)])
	# @variable(m, plprod[ts=timestep, np=1:size(nplant,1)], start = 0.0)
	# @variable(m, plprodtot[np=1:size(nplant,1)], start = 0.0)
	
	# bidding plant constraints
	# calculate plant prod from energy consumption for each timestep at corresponding bus
	@constraint(m, [ts=timestep, np=1:size(nplant,1)], plprod[ts,np] == plant[!,:plantprod][np].*pplant[ts,np])
	# ramping constraint for plant production
	@NLconstraint(m, [ts=timestep[2:end], np=1:size(nplant,1)], -plant[!,:plantRRL][np] <= pplant[ts,np] - pplant[ts-1,np] <= plant[!,:plantRRL][np])
	# temporal constraint: plant production within demand and flex*demand 
	@constraint(m, [np=1:size(nplant,1)], plprodtot[np] == sum(plprod[ts,np] for ts=timestep))
	@constraint(m, [np=1:size(nplant,1)], plant[!,:plantdemand][np]*timestep[end]/24 <= plprodtot[np] <= (1+plant[!,:plantflex][np]).*plant[!,:plantdemand][np]*timestep[end]/24)

	#active power balance with LMP=dual
	@constraint(m, PBPC[ts = timestep], sum(p_g[ts,g] for g=1:size(G[ts], 1)) - sum(N[ts][!,:Pd][nd] for nd=1:size(N[ts][!,:Pd],1))
	- sum(pplant[ts,np] for np=1:size(nplant,1))
	# - N[ts][!,:Gsh][plant[!,:busidx][np]]*v_n[ts,plant[!,:busidx][np]]^2 
			# - sum(p_li[ts,Lidxn[plant[!,:busidx][np]][lix],Lin[plant[!,:busidx][np]][lix]] for lix=1:length(Lidxn[plant[!,:busidx][np]])) 
			# - sum(p_ti[ts,Tidxn[plant[!,:busidx][np]][tix],Tin[plant[!,:busidx][np]][tix]] for tix=1:length(Tidxn[plant[!,:busidx][np]])) 
			== sum(pslackp_n[ts,nd] - pslackm_n[ts,nd] for nd=1:size(N[ts][!,:Pd],1)))
		
	#reactive power balance 
	@constraint(m, PBPC2[ts = timestep], sum(q_g[ts,g] for g=1:size(G[ts], 1)) - sum(N[ts][!,:Qd][nd] for nd=1:size(N[ts][!,:Qd],1))
					== sum(qslackp_n[ts,nd] - qslackm_n[ts,nd] for nd=1:size(N[ts][!,:Qd],1)))	

	# constrain total power generation
	
	# @constraint(m, [ts=timestep], sum(p_g[ts,g] for g=1:size(G[ts], 1)) <= GenTot[ts]/MVAbase)
	GenTot = 1474640.6077464963
	# @constraint(m, sum(p_g[ts,g] for ts=timestep, g=1:size(G[ts], 1)) <= GenTot/MVAbase)

	# power flow constraints
	for ts = timestep
		# fix angle at reference bus to zero
		JuMP.fix(theta_n[ts,RefBus], 0.0, force=true)
		
		# add power flow constraints for base case
		addpowerflowcons!(m, v_n[ts,:], theta_n[ts,:], p_li[ts,:,:], q_li[ts,:,:], p_ti[ts,:,:], q_ti[ts,:,:], b_s[ts,:], p_g[ts,:], q_g[ts,:],
			pslackm_n[ts,:], pslackp_n[ts,:], qslackm_n[ts,:], qslackp_n[ts,:], sslack_li[ts,:,:], sslack_ti[ts,:,:],
		 	N[ts], L[ts], T[ts], (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
			SShn, Gn, K_outidx), o.ThermalLimits)
	
		# register approximation functions
		sclamp = nothing
		sstep = nothing
	end
	

	println("start ramping constraints")
	flush(stdout)
	for ts = timestep[2:end]
		#ramping constraints
		@NLconstraint(m, [g=1:size(G[ts], 1)], -1*G[ts][!,:Pub][g] <= p_g[ts,g]-p_g[ts-1,g] <= 1*G[ts][!,:Pub][g])
	end


	## objective function
	prodcost = zero(JuMP.AffExpr)
	basepenalty = zero(JuMP.QuadExpr)
	plantrev = zero(JuMP.AffExpr)

	println("start prod cost")
	for ts = timestep

		# production cost
		for g = 1:size(G[ts], 1)
	    	npairs = length(G[ts][!,:CostPi][g])
            	@assert sort(G[ts][!,:CostPi][g]) == G[ts][!,:CostPi][g]
 
            	t_h0 = getstartcostcoefficients(JuMP.start_value(p_g[ts,g]), G[ts][!,:CostPi][g], G[ts][!,:Plb][g], G[ts][!,:Pub][g])
            	t_h = @variable(m, [i=1:npairs], lower_bound=0.0, start=t_h0[i])
            
 	    	@constraint(m, sum(t_h[h] for h=1:npairs) == 1.0)
	    	@constraint(m, p_g[ts,g] == sum(G[ts][!,:CostPi][g][h]*t_h[h] for h=1:npairs))
            gencost = @expression(m, sum(G[ts][!,:CostCi][g][h]*t_h[h] for h=1:npairs))
	    	add_to_expression!(prodcost, gencost)
		end
	
	end

	println("start base case pen")
	for ts = timestep
		# base case penalties
		basecasepenalty = quadpenaltyfunctions(pslackm_n[ts,:], pslackp_n[ts,:], qslackm_n[ts,:], qslackp_n[ts,:],
		sslack_li[ts,:,:], sslack_ti[ts,:,:], P[ts])
		add_to_expression!(basepenalty, basecasepenalty)
		
	end
	
	#plant prod revenue
	for np=1:size(nplant,1)
		add_to_expression!(plantrev, sum(pplant[ts,np].*plant[!,:plantrev][np] for ts=timestep))
	end
	


	# declare objective
	@objective(m, Min, (prodcost + DELTA*basepenalty - plantrev))
	
	## attempt to solve SCACOPF
	JuMP.optimize!(m)

        # termination_status, primal_status, and dual_status: http://www.juliaopt.org/JuMP.jl/dev/solutions/
        solvestatusprimal = JuMP.primal_status(m)
        solvestatus = JuMP.termination_status(m)
        println("JuMP optimize MOI termination_status: ", solvestatus, " ", solvestatusprimal)
	if (solvestatusprimal != MOI.FEASIBLE_POINT) && (solvestatusprimal != MOI.NEARLY_FEASIBLE_POINT)
		error("solver failed to find a feasible solution")
	end
	
	## objective breakdown
	#prodcostvalue = JuMP.value(prodcost)
	penaltybasevalue = DELTA*JuMP.value(basepenalty)
	plantrevvalue = JuMP.value(plantrev)
	objectivevalue = objective_value(m)
	prodcostvalue = objectivevalue - penaltybasevalue + plantrevvalue

	
	println("Objective:           \$", round(objectivevalue, digits=1), "\n",
		"Generation cost:     \$", round(prodcostvalue, digits=1), "\n",
		"Penalty base:        \$", round(penaltybasevalue, digits=1), "\n",
		"Plant revenue:        \$", round(plantrevvalue, digits=1), "\n")

	flush(stdout)

	## outputs
	LMPs = Containers.DenseAxisArray{Float64}(undef, timestep)
	for ts = timestep
			LMPs[ts] = dual(PBPC[ts])/MVAbase
	end
	println("LMPs:   ", LMPs, "\n")

	p_gv = JuMP.value.(p_g)
    Gent = zeros(size(timestep,1))
    Gentot = 0.0 
    for ts = timestep
        Gent[ts] = sum(p_gv[ts,g] for g=1:size(G[ts], 1))
        Gentot += Gent[ts]
    end
    Gentot *= MVAbase
    println("Total generation:  ", Gentot)

	plantpower = JuMP.value.(pplant).*MVAbase
	println("first plant power:      ", plantpower[:,1], "\n")

	# get values of line losses from line and transformer powers, and power generation sum
	p_liv = JuMP.value.(p_li)
	# println(typeof(p_liv), size(p_liv))
	p_tiv = JuMP.value.(p_ti)
	p_gv = JuMP.value.(p_g)
	Lloss = Containers.DenseAxisArray{Float64}(undef, timestep, 1:size(N[timestep[1]], 1))
	Genb = Containers.DenseAxisArray{Float64}(undef, timestep, 1:size(N[timestep[1]], 1))
	Llosst = zeros(size(timestep,1))
	Gent = zeros(size(timestep,1))
	eta = zeros(size(timestep,1))
	for ts = timestep
		for n = 1:size(N[timestep[1]], 1)
			Genb[ts,n] = sum(p_gv[ts,g] for g=Gn[n]; init=0)
			Lloss[ts,n] = sum(p_liv[ts,Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n]); init=0) +
			sum(p_tiv[ts,Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n]); init=0)
		end
		Llosst[ts] = sum(Lloss[ts,ns] for ns = 1:size(N[timestep[1]], 1); init=0)*MVAbase
		Gent[ts] = sum(Genb[ts,ns] for ns = 1:size(N[timestep[1]], 1); init=0)*MVAbase
		eta[ts] = Llosst[ts]/Gent[ts]
	end
	# println("Percent line loss[1:3] ", eta[1:3])


	# get values of slack variables
	slackpm = JuMP.value.(pslackm_n)
	slackpp = JuMP.value.(pslackp_n)
	slackqm = JuMP.value.(qslackm_n)
	slackqp = JuMP.value.(qslackp_n)
	slackl = JuMP.value.(sslack_li)
	slackt = JuMP.value.(sslack_ti)


        if false
            kk=2
            vp_k = JuMP.value.(p_gk)
            vp_k = view(vp_k, :, kk)
            vp   = JuMP.value.(p_g)
            vPlb = G[!,:Plb]
            vPub = G[!,:Pub]
            delta= JuMP.value.(delta_k)
            delta= delta[kk]
            alpha = G[!,:alpha]
            for i=1:length(vp_k)
                @printf("g=%4d pg=%12.5e  pg_k=%12.5e diffPlb=%12.5e diffPub=%12.5e  p_k-(p+alpha*delta)=%12.5e fPlb=%12.5e Pub=%12.5e alpha=%12.5e  delta=%12.5e\n", 
                        i, vp[i], vp_k[i], vp_k[i]-vPlb[i], vPub[i]-vp_k[i], vp_k[i]-(vp[i]+alpha[i]*delta), vPlb[i], vPub[i], alpha[i], delta)
                #@printf("alpha[i]=%g  delta=%g\n", alpha[i], delta)
            end
            #println("G Plb", G[!,:Plb])
            #println("G Pub", G[!,:Pub])
        end

        if false
            vq_gk = JuMP.value.(q_gk)
            vq_gk = view(vq_gk,:,1)
            vq_g  = JuMP.value.(q_g)
            vv_nk = JuMP.value.(v_nk)
            vv_nk = view(vv_nk,:,1)
            vv_n  = JuMP.value.(v_n)
            vQlb = G[!,:Qlb]
            vQub = G[!,:Qub]
            
            for i=1:length(vq_gk)
                @printf("PVPQ gen=%3d q_g=%12.5e q_gk=%12.5e v_n=%12.5e v_nk=%12.5e | v-v_k=%12.5e Qub-q_k=%12.5e q_k-Qlb=%12.5e | Qlb=%12.5e Qub=%12.5e\n",
                        i, vq_g[i], vq_gk[i], vv_n[i], vv_nk[i], vv_n[i]- vv_nk[i],  vQub[i]-vq_gk[i],  vq_gk[i]-vQlb[i], vQlb[i], vQub[i])
            end
                    
        end

        if false
            #println(JuMP.value.(v_n))
            println("p_g", JuMP.value.(p_g))
            println("p_")
            println("G Plb", G[!,:Plb])
            println("G Pub", G[!,:Pub])


            println("q_g", JuMP.value.(q_g))
            println("G Qlb", G[!,:Qlb])
            println("G Qub", G[!,:Qub])

            #println("N[!,:PD]", N[!,:Pd])
            #println("N[!,:QD]", N[!,:Qd])
            
            qslackm_nv = JuMP.value.(qslackm_n)
            qslackp_nv = JuMP.value.(qslackp_n)
	    #@variable(m, pslackm_n[n=1:size(N, 1)] >= 0)
            println("pslackm_n", JuMP.value.(pslackm_n))
	    #@variable(m, pslackp_n[n=1:size(N, 1)] >= 0)
            println("pslackp_n", JuMP.value.(pslackp_n))
	    #@variable(m, qslackm_n[n=1:size(N, 1)] >= 0)
            println("qslackm_n", JuMP.value.(qslackm_n))
	    #@variable(m, qslackp_n[n=1:size(N, 1)] >= 0)
            println("qslackm_n", JuMP.value.(qslackp_n))
	    #@variable(m, sslack_li[l=1:size(L, 1), i=1:2] >= 0)
	    #@variable(m, sslack_ti[t=1:size(T, 1), i=1:2] >= 0)
            p_liv = JuMP.value.(p_li)
            p_tiv = JuMP.value.(p_ti)
            q_liv = JuMP.value.(q_li)
            q_tiv = JuMP.value.(q_ti)
            q_gv = JuMP.value.(q_g)
            p_gv = JuMP.value.(p_g)
            b_sv = JuMP.value.(b_s)
            println("b_sv: ", b_sv)
	    for n = 1:size(N, 1)
                #		@NLconstraint(m, sum(q_g[g] for g=Gn[n]) - N[!,:Qd][n] -
		#	(N[!,:Bsh][n] + sum(b_s[s] for s=SShn[n]))*v_n[n]^2 -
		#	sum(q_li[Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n])) -
		#	sum(q_ti[Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n])) ==
		#	qslackp_n[n] - qslackm_n[n])
                sumpg = 0.
                if length(Gn[n])>0
                    sumpg = sum(q_gv[g] for g=Gn[n])
                end
                sumqli=0.;
                if length(Lidxn[n])>0
                    sumqli = sum(q_liv[Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n]))
                end
                sumqti=0.
                if length(Tidxn[n])>0
                    sumqti = sum(q_tiv[Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n]))
                end
		NBsh = N[!,:Bsh][n]
                if length(SShn[n])>0
                    NBsh = NBsh +  sum(b_sv[s] for s=SShn[n])
                end
                if qslackm_nv[n]>1e-6 || qslackp_nv[n]>1e-6
                    @printf("n=%d  sumpg=%12.5e  N[!,:Qd]=%12.5e  N[!,:Bsh]=%12.5e sum_qli=%12.5e sum_qti=%12.5e slackm=%12.5e  slackp=%12.5e\n", n, sumpg, N[!,:Qd][n], NBsh, sumqli, sumqti, qslackm_nv[n], qslackp_nv[n])
                end
            end
        end
	
	return
	## return solution
	# return MPACOPFState(MPBaseState(JuMP.value.(v_n), 
    #                             	JuMP.value.(theta_n), 
    #                                 JuMP.value.(b_s),
    #                         JuMP.value.(p_g), 
    #                                   JuMP.value.(q_g)),
	# 								  biddingRes(prodcostvalue, penaltybasevalue, plantrevvalue, plantpower,LMPs,LMPsN),
	# 								  slackv(slackpm, slackpp, slackqm, slackqp, slackl, slackt),
	# 								  TranLoss(Llosst, Gent, eta, p_liv, p_tiv))
	
end

end
