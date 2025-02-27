#__precompile__()

module MPACOPFbus

## elements to be exported

export SolveMPACOPFbus

## load external modules

using DataFrames, Roots, JuMP, Printf, SparseArrays, LinearAlgebra, MathOptInterface
const MOI = JuMP.MOI
const MOIU = JuMP.MOIU
import JuMP: add_to_expression!

## load internal modules and functions
Base.include(@__MODULE__,"MPACOPF.jl")
Base.include(@__MODULE__,"InstanceReaderMPbus.jl")
Base.include(@__MODULE__,"SmoothApproximations.jl")
Base.include(@__MODULE__,"GoUtils.jl")
Base.include(@__MODULE__,"SolutionEvaluator.jl")
using  .GoUtils, .SmoothApproximations, .InstanceReaderMPbus, .MPACOPF
import .SolutionEvaluator: InitialBaseCaseSolution, BaseCaseSolution, InitialContingenciesSolutionsFromBase

import DataStructures: OrderedDict

## GoOptions struct
include("GoOptions.jl")
include("MPACOPFbusStates.jl")
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
	IndexSets::Tuple, ThermalLimits::Symbol;
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
	
	# percent of relaxation of transmission constraint upper bounds
	PerUpper = 1 + 0.01*0
	# thermal limits and power flows -- lines
	for l=1:size(L, 1), i=1:2
		if SysCond == :Contingency && ConType == :Line && outidx == l
			continue
		end
		# if ThermalLimits == :quadr
		# 	@constraint(m, p_li[l,i]^2 + q_li[l,i]^2 <=
		# 		((L[!, RateSymb][l]*v_n[L_Nidx[l,i]])*PerUpper + sslack_li[l,i])^2)
		# else
		# 	@NLconstraint(m, (p_li[l,i]^2 + q_li[l,i]^2+1e-6)^0.5 <=
		# 		(L[!, RateSymb][l]*v_n[L_Nidx[l,i]])*PerUpper + sslack_li[l,i] )
		# end
		addpowerflowcon!(m, p_li[l,i], v_n[L_Nidx[l,i]], v_n[L_Nidx[l,3-i]],
			theta_n[L_Nidx[l,i]], theta_n[L_Nidx[l,3-i]],
			L[!,:G][l], -L[!,:G][l], -L[!,:B][l])
		addpowerflowcon!(m, q_li[l,i], v_n[L_Nidx[l,i]], v_n[L_Nidx[l,3-i]],
			theta_n[L_Nidx[l,i]], theta_n[L_Nidx[l,3-i]],
			-L[!,:B][l]-L[!,:Bch][l]/2, L[!,:B][l], -L[!,:G][l])
	end
	
	# thermal limits and power flows -- transformers
	for t=1:size(T, 1)
		if SysCond == :Contingency && ConType == :Transformer && outidx == t
			continue
		end
		# if ThermalLimits == :quadr
		# 	@constraint(m, [i=1:2], p_ti[t,i]^2 + q_ti[t,i]^2 <=
		# 		((T[!, RateSymb][t])*PerUpper + sslack_ti[t,i])^2 )
		# else
		# 	@NLconstraint(m, [i=1:2], (p_ti[t,i]^2 + q_ti[t,i]^2+1e-6)^0.5 <=
		# 		(T[!, RateSymb][t])*PerUpper + sslack_ti[t,i] )
		# end
		addpowerflowcon!(m, p_ti[t,1], v_n[T_Nidx[t,1]], v_n[T_Nidx[t,2]],
			theta_n[T_Nidx[t,1]], theta_n[T_Nidx[t,2]],
			T[!,:G][t]/T[!,:Tau][t]^2+T[!,:Gm][t], -T[!,:G][t]/T[!,:Tau][t],
			-T[!,:B][t]/T[!,:Tau][t], -T[!,:Theta][t])
		addpowerflowcon!(m, q_ti[t,1], v_n[T_Nidx[t,1]], v_n[T_Nidx[t,2]],
			theta_n[T_Nidx[t,1]], theta_n[T_Nidx[t,2]],
			-T[!,:B][t]/T[!,:Tau][t]^2-T[!,:Bm][t], T[!,:B][t]/T[!,:Tau][t],
			-T[!,:G][t]/T[!,:Tau][t], -T[!,:Theta][t])
		addpowerflowcon!(m, p_ti[t,2], v_n[T_Nidx[t,2]], v_n[T_Nidx[t,1]],
			theta_n[T_Nidx[t,2]], theta_n[T_Nidx[t,1]],
			T[!,:G][t], -T[!,:G][t]/T[!,:Tau][t],
			-T[!,:B][t]/T[!,:Tau][t], T[!,:Theta][t])
		addpowerflowcon!(m, q_ti[t,2], v_n[T_Nidx[t,2]], v_n[T_Nidx[t,1]],
			theta_n[T_Nidx[t,2]], theta_n[T_Nidx[t,1]],
			-T[!,:B][t], T[!,:B][t]/T[!,:Tau][t],
			-T[!,:G][t]/T[!,:Tau][t], T[!,:Theta][t])
	end
	
	# balance
	# for n = 1:size(N, 1)
		# PBPC = @constraint(m, PBPC[n = 1:size(N, 1)], sum(p_g[g] for g=Gn[n]) - N[!,:Pd][n] - N[!,:Gsh][n]*v_n[n]^2 -
		# # @constraint(m, sum(p_g[g] for g=Gn[n]) - N[!,:Pd][n] - N[!,:Gsh][n]*v_n[n]^2 -
		# 	pplant[n] - 
		# 	sum(p_li[Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n])) -
		# 	sum(p_ti[Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n])) .==
		# 	pslackp_n[n] - pslackm_n[n])
			# set_name(PBPC,"PBPC_$(ts)_$(n)")
	for n = 1:size(N, 1)
		@NLconstraint(m, sum(q_g[g] for g=Gn[n]) - N[!,:Qd][n] -
			(-N[!,:Bsh][n] - sum(b_s[s] for s=SShn[n]))*v_n[n]^2 -
			sum(q_li[Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n])) -
			sum(q_ti[Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n])) ==
			qslackp_n[n] - qslackm_n[n])
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
	# qslackm_n::AbstractVector, qslackp_n::AbstractVector,
	# sslack_li::AbstractArray, sslack_ti::AbstractArray,
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
		return a, b
	end
	aP, bP = quadcoeffs(findfirst(P[!,:Slack] .== :P))
	# aQ, bQ = quadcoeffs(findfirst(P[!,:Slack] .== :Q))
	# aS, bS = quadcoeffs(findfirst(P[!,:Slack] .== :S))
	
	# imbalance penalties
	addquadpenaltyblock!(casepenalty, pslackm_n, aP, bP)
	addquadpenaltyblock!(casepenalty, pslackp_n, aP, bP)
	# addquadpenaltyblock!(casepenalty, qslackm_n, aQ, bQ)
	# addquadpenaltyblock!(casepenalty, qslackp_n, aQ, bQ)
	
	# # thermal overload
	# addquadpenaltyblock!(casepenalty, sslack_li, aS, bS)
	# addquadpenaltyblock!(casepenalty, sslack_ti, aS, bS)

	# # zero out plant power on other buses
	# addquadpenaltyblock!(casepenalty, plantslack, aP, bP)

	# return expression containing all penalty terms
	return casepenalty
	
end


## function for solving MP AC OPF
function SolveMPACOPFbus(o::GoOptions, N::Vector{Any}, L::Vector{Any}, T::Vector{Any},
	              SSh::Vector{Any}, G::Vector{Any}, K::Vector{Any}, P::Vector{Any}, timestep::Vector{Int64}, 
                  vLoad::Vector{Any},plant::DataFrame, MVAbase::Float64, NLSolver; 
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
	# @variable(m, N[ts][!,:Vlb][n] <= v_n[ts=timestep, n=1:size(N[timestep[1]], 1)] <= N[ts][!,:Vub][n], start = v0_n[n])
	# @variable(m, theta_n[ts=timestep, n=1:size(N[timestep[1]], 1)], start = theta0_n[n])
	# @variable(m, p_li[ts=timestep, l=1:size(L[timestep[1]], 1), i=1:2], start = p0_li[l,i])
	# @variable(m, q_li[ts=timestep, l=1:size(L[timestep[1]], 1), i=1:2], start = q0_li[l,i])
	# @variable(m, p_ti[ts=timestep, t=1:size(T[timestep[1]], 1), i=1:2], start = p0_ti[t,i])
	# @variable(m, q_ti[ts=timestep, t=1:size(T[timestep[1]], 1), i=1:2], start = q0_ti[t,i])
	# @variable(m, SSh[ts][!,:Blb][s] <= b_s[ts=timestep, s=1:size(SSh[timestep[1]], 1)] <= SSh[ts][!,:Bub][s], start = b0_s[s])
	@variable(m, G[ts][!,:Plb][g] <= p_g[ts=timestep, g=1:size(G[timestep[1]], 1)] <= G[ts][!,:Pub][g], start = p0_g[g])
	@variable(m, G[ts][!,:Qlb][g] <= q_g[ts=timestep, g=1:size(G[timestep[1]], 1)] <= G[ts][!,:Qub][g], start = q0_g[g])
	@variable(m, pslackm_n[ts=timestep] >= 0, start = pslack0m_n[1])
	@variable(m, pslackp_n[ts=timestep] >= 0, start = pslack0p_n[1])
	@variable(m, qslackm_n[ts=timestep] >= 0, start = qslack0m_n[1])
	@variable(m, qslackp_n[ts=timestep] >= 0, start = qslack0p_n[1])
	# @variable(m, qslackm_n[ts=timestep, n=1:size(N[timestep[1]], 1)] >= 0, start = qslack0m_n[n])
	# # @variable(m, qslackp_n[ts=timestep, n=1:size(N[timestep[1]], 1)] >= 0, start = qslack0p_n[n])
	# @variable(m, sslack_li[ts=timestep, l=1:size(L[timestep[1]], 1), i=1:2] >= 0, start = sslack0_li[l,i])
	# @variable(m, sslack_ti[ts=timestep, t=1:size(T[timestep[1]], 1), i=1:2] >= 0, start = sslack0_ti[t,i])

	# bidding entity variables at timestep
	@variable(m, plant[!,:pplantlb][1] <= pplant[ts=timestep] <= plant[!,:pplantub][1])
	@variable(m, plprod[ts=timestep])
	@variable(m, plprodtot)

	
	# bidding plant constraints
	# calculate plant prod from energy consumption for each timestep
	@constraint(m, [ts=timestep], plprod[ts] == plant[!,:plantprod][1]*pplant[ts])
	# ramping constraint for plant production
	@NLconstraint(m, [ts=timestep[2:end]], -plant[!,:plantRRL][1] <= pplant[ts] - pplant[ts-1] <= plant[!,:plantRRL][1])
	# temporal constraint: plant production within demand and flex*demand 
	@constraint(m, plprodtot == sum(plprod[ts] for ts=timestep))
	@constraint(m, plant[!,:plantdemand][1]*timestep[end]/24 <= plprodtot <= (1+plant[!,:plantflex][1])*plant[!,:plantdemand][1]*timestep[end]/24)

    # overall active power balance at each time step 
    @constraint(m, PBPC[ts = timestep], sum(p_g[ts,g] for g=1:size(G[ts], 1)) - vLoad[ts] -
    		pplant[ts] == pslackp_n[ts]-pslackm_n[ts])

	# overall reactive power balance at each time step 
	@constraint(m, [ts = timestep], sum(q_g[ts,g] for g=1:size(G[ts], 1)) - sum(N[ts][!,:Qd][n] for n=1:size(N, 1)) 
			== qslackp_n[ts]-qslackm_n[ts])



	# constrain total power generation
	# GenTot = 1.446655699976021e6
	# @constraint(m, sum(p_g[ts,g] for ts=timestep, g=1:size(G[ts], 1)) <= GenTot/MVAbase)	

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
			prodcoststep[ts] += gencost
		end
	end

    #plant prod revenue
    add_to_expression!(plantrev, sum(pplant[ts]*plant[!,:plantrev][1] for ts=timestep))
	# for np=1:size(nplant,1)
	# 	add_to_expression!(plantrev, sum(pplant[ts,np].*plant[!,:plantrev][np] for ts=timestep))
	# end

	println("start base case pen")
    	# compute coefficients (a*x^2 + b*x)
	function quadcoeffs(idx)
		perm = sortperm(P[1][!,:Penalties][idx])
		quantities = view(P[1][!,:Quantities][idx], perm)
		penalties = view(P[1][!,:Penalties][idx], perm)
		b = penalties[1]
		a = (penalties[2] - b)/(2*quantities[1])
		@assert b > 0 && a > 0
		return a, b
	end
	aP, bP = quadcoeffs(findfirst(P[1][!,:Slack] .== :P))
	aQ, bQ = quadcoeffs(findfirst(P[1][!,:Slack] .== :Q))
    addquadpenaltyblock!(basepenalty, pslackm_n, aP, bP)
    addquadpenaltyblock!(basepenalty, pslackp_n, aP, bP)
	addquadpenaltyblock!(basepenalty, qslackm_n, aQ, bQ)
	addquadpenaltyblock!(basepenalty, qslackp_n, aQ, bQ)
    # add_to_expression!(basepenalty, basecasepenalty)


	# declare objective
	@objective(m, Min, (prodcost + DELTA*basepenalty - plantrev))
	# @objective(m, Min, (prodcost + DELTA*basepenalty))
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
	# prodcostvalue = objectivevalue - penaltybasevalue
	
	println("Objective:           \$", round(objectivevalue, digits=1), "\n",
		"Generation cost:     \$", round(prodcostvalue, digits=1), "\n",
		"Penalty base:        \$", round(penaltybasevalue, digits=1), "\n",
		"Plant revenue:        \$", round(plantrevvalue, digits=1), "\n")

	flush(stdout)


	## outputs
	plantpower = JuMP.value.(pplant).*MVAbase
	println("plant power:      ", plantpower, "\n")

	
	LMPs = Containers.DenseAxisArray{Float64}(undef, timestep)
	for ts = timestep
			LMPs[ts] = dual(PBPC[ts])/MVAbase
	end
	println("LMPs:   ", LMPs, "\n")

	@assert penaltybasevalue <= 0
	
    p_gv = JuMP.value.(p_g)
    Gent = zeros(size(timestep,1))
    Gentot = 0.0 
    for ts = timestep
        Gent[ts] = sum(p_gv[ts,g] for g=1:size(G[ts], 1))*MVAbase
        Gentot += Gent[ts]
    end
    # Gentot *= MVAbase
    println("Total generation:  ", Gentot)

	GenCostStep = JuMP.value.(prodcoststep)
	println("GenCost @ timestep ", GenCostStep)

	# get values of slack variables
	# slackpm = JuMP.value.(pslackm_n)
	slackpp = JuMP.value.(pslackp_n)
	# slackqm = JuMP.value.(qslackm_n)
	slackqp = JuMP.value.(qslackp_n)
	# slackl = JuMP.value.(sslack_li)
	# slackt = JuMP.value.(sslack_ti)


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
	
	# return solution
	return MPbusState(plantpower, LMPs, Gent, GenCostStep)
	
end

end
