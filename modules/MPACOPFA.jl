#__precompile__()

module MPACOPFA

## elements to be exported
export SolveMPACOPFA

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
		return a, b
	end
	aP, bP = quadcoeffs(findfirst(P[!,:Slack] .== :P))
	aQ, bQ = quadcoeffs(findfirst(P[!,:Slack] .== :Q))
	aS, bS = quadcoeffs(findfirst(P[!,:Slack] .== :S))
	
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

## function for solving MP AC OPF
function SolveMPACOPFA(o::GoOptions, N::Vector{Any}, L::Vector{Any}, T::Vector{Any},
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
	GenTot =  1.4769745621171894e6
	@constraint(m, sum(p_g[ts,g] for ts=timestep, g=1:size(G[ts], 1)) <= GenTot/MVAbase)

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

    GenCostStep = JuMP.value.(prodcoststep)
	println("GenCost @ timestep ", GenCostStep)

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
	println("active transformer constraints: ", TieqnA)

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

end