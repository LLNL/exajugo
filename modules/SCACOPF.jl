#__precompile__()

module SCACOPF

## elements to be exported

export SolveSCACOPF, SolveContingencyFastBackend, SolveContingencyTwoStep, SolveContingencyBackend #,SolveAllContingencies

## load external modules

using DataFrames, Roots, JuMP, Printf
const MOI = JuMP.MOI
const MOIU = JuMP.MOIU
import JuMP: add_to_expression!

## load internal modules and functions
Base.include(@__MODULE__,"InstanceReader.jl")
Base.include(@__MODULE__,"SmoothApproximations.jl")
Base.include(@__MODULE__,"GoUtils.jl")
Base.include(@__MODULE__,"SolutionEvaluator.jl")
using .InstanceReader, .GoUtils, .SmoothApproximations
import .SolutionEvaluator: InitialBaseCaseSolution, BaseCaseSolution, InitialContingenciesSolutionsFromBase

import DataStructures: OrderedDict

## GoOptions struct
include("GoOptions.jl")
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
	
	# thermal limits and power flows -- lines
	for l=1:size(L, 1), i=1:2
		if SysCond == :Contingency && ConType == :Line && outidx == l
			continue
		end
		if ThermalLimits == :quadr
			@constraint(m, p_li[l,i]^2 + q_li[l,i]^2 <=
				(L[!, RateSymb][l]*v_n[L_Nidx[l,i]] + sslack_li[l,i])^2)
		else
			@NLconstraint(m, (p_li[l,i]^2 + q_li[l,i]^2+1e-6)^0.5 <=
				L[!, RateSymb][l]*v_n[L_Nidx[l,i]] + sslack_li[l,i])
		end
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
		if ThermalLimits == :quadr
			@constraint(m, [i=1:2], p_ti[t,i]^2 + q_ti[t,i]^2 <=
				(T[!, RateSymb][t] + sslack_ti[t,i])^2)
		else
			@NLconstraint(m, [i=1:2], (p_ti[t,i]^2 + q_ti[t,i]^2+1e-6)^0.5 <=
				T[!, RateSymb][t] + sslack_ti[t,i])
		end
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
	for n = 1:size(N, 1)
		@constraint(m, sum(p_g[g] for g=Gn[n]) - N[!,:Pd][n] - N[!,:Gsh][n]*v_n[n]^2 -
			sum(p_li[Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n])) -
			sum(p_ti[Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n])) ==
			pslackp_n[n] - pslackm_n[n])
		@NLconstraint(m, sum(q_g[g] for g=Gn[n]) - N[!,:Qd][n] -
			(-N[!,:Bsh][n] - sum(b_s[s] for s=SShn[n]))*v_n[n]^2 -
			sum(q_li[Lidxn[n][lix],Lin[n][lix]] for lix=1:length(Lidxn[n])) -
			sum(q_ti[Tidxn[n][tix],Tin[n][tix]] for tix=1:length(Tidxn[n])) ==
			qslackp_n[n] - qslackm_n[n])
	end
	
end

## function to add AGC coupling constraint
function addAGCconstraint!(m::JuMP.Model, o::GoOptions, p::Union{Real, JuMP.VariableRef},
	delta_k::Union{Real, JuMP.VariableRef}, p_k::Union{Real, JuMP.VariableRef},
	Plb::Real, Pub::Real, alpha::Real;
	sclamp::Union{Function, Nothing}=nothing, SetStart::Bool=true)

        if Pub-Plb<1e-8
            #s = @sprintf("Plb=%g and Pub=%g close to each other, not enforcing  AGCconstraint", Plb, Pub)
            #@warn(s)
            return 0
        end
	
	# formulate approximate AGC coupling
	if o.CouplingMode == :approx
		@NLconstraint(m, p_k == Plb + (Pub - Plb)*sclamp((p + alpha*delta_k - Plb)/(Pub - Plb)))
		# note this will fail if sclamp has not been registered! not using try/catch for perfomance
		return nothing
	end
        gen_band = Pub-Plb
        reg = o.SmoothingParamCoupling

	# declare exact AGC coupling auxiliary variables
        # these are unitfree and implicitly <=1
        rhom = @variable(m, lower_bound=0)
        rhop = @variable(m, lower_bound=0)
	if SetStart
            aux = (getstartvalue(p_k) - getstartvalue(p) - alpha*getstartvalue(delta_k))/gen_band

            #p_kv = getstartvalue(p_k)
            rhomv = max(0,  aux)
            rhopv = max(0, -aux)

	    JuMP.set_start_value(rhom, rhomv)
	    JuMP.set_start_value(rhop, rhopv)
	end
	
	# add exact AGC coupling constraints to the model
	@constraint(m, p_k + gen_band*rhop - gen_band*rhom == p + alpha*delta_k)

	if o.CouplingMode == :complementarity
		@NLconstraint(m, rhop*(p_k - Pub) == 0)
		@NLconstraint(m, rhom*(p_k - Plb) == 0)

	elseif o.CouplingMode == :complementarityapprox
            reg = o.SmoothingParamCoupling
            @constraint(m, -reg <= rhop*(p_k/gen_band-Pub/gen_band) <= reg)
            @constraint(m, -reg <= rhom*(p_k/gen_band-Plb/gen_band) <= reg)
            #@constraint(m, rhop * rhom <= reg)
	elseif o.CouplingMode == :fischerburm

            reg = o.SmoothingParamCoupling
	    @NLconstraint(m, ( ((+p_k - Plb)/gen_band)^2 + rhom^2 + reg)^0.5  >= ( p_k - Plb)/gen_band + rhom)
            @NLconstraint(m, ( ((-p_k + Pub)/gen_band)^2 + rhop^2 + reg)^0.5  >= (-p_k + Pub)/gen_band + rhop)

            #reg = max(1e-8, reg^2)
            @constraint(m, rhop * rhom <= reg )

            if false
                viol = ( +getstartvalue(p_k) - Plb)/gen_band + getstartvalue(rhom) - 
                    (((getstartvalue(p_k) - Plb)/gen_band)^2 + getstartvalue(rhom)^2 + reg)^0.5
                if false #abs(viol)>1e-8
                    @printf( "Plb=%11.4e p_k-Plb=%11.4e rhom=%11.4e | Pub=%11.4e Pub-p_k=%11.4e rhop=%11.4e | reg=%11.4e   lhs=-- >= rhs=--  viol=%11.4e\n", 
                             Plb, (+getstartvalue(p_k) - Plb),  getstartvalue(rhom),
                             Pub, (-getstartvalue(p_k) + Pub),  getstartvalue(rhop),                     
                             reg,
                             abs(viol))
                end
            end
	elseif o.CouplingMode == :fischerburmlifted
            tu = @variable(m, lower_bound=0)
            @constraint(m, tu^2 <= (-p_k+Pub)^2 + rhop^2)
            @constraint(m, tu   >=  -p_k+Pub    + rhop  )

            tl = @variable(m, lower_bound=0)
            @constraint(m, tl^2 <= ( p_k-Plb)^2 + rhom^2)
            @constraint(m, tl   >=   p_k+Plb    + rhom  )

            if SetStart
                JuMP.set_start_value(tu, 1)
                JuMP.set_start_value(tl, 1)
            end
	elseif o.CouplingMode == :fischerburmsquare
		@NLconstraint(m, ( ( (+p_k - Plb)^2 + rhom^2 )^0.5 - (+p_k - Plb + rhom) )^2 == 0)
		@NLconstraint(m, ( ( (-p_k + Pub)^2 + rhop^2 )^0.5 - (-p_k + Pub + rhop) )^2 == 0)
	else
	    warn("unrecognized coupling mode")
            return 0
	end
	return 1
	
end

## function to add PV/PQ coupling constraint

function addPVPQconstraint!(m::JuMP.Model, o::GoOptions,
	q_k::Union{Real, JuMP.VariableRef, JuMP.NonlinearExpression},
	v::Union{Real, JuMP.VariableRef}, v_k::Union{Real, JuMP.VariableRef},
	Qlb::Real, Qub::Real;
	sstep::Union{Function, Nothing}=nothing, SetStart::Bool=true)

    gen_band = Qub-Qlb
    if gen_band<1e-8
        #s = @sprintf("Qlb=%g and Qub=%g close to each other, not enforcing  PVPQconstraint", Qlb, Qub); @warn(s)
        return 0
    end

	# formulate approximate PV/PQ coupling
	if o.CouplingMode == :approx
		@NLconstraint(m, q_k == Qub + (Qlb - Qub)*sstep(v_k - v))
		# note this will fail if sstep has not been registered! not using try/catch for perfomance
		return
	end

	# declare exact PV/PQ coupling auxiliary variables
        num = @variable(m, lower_bound=0)#, upper_bound=1.2)
        nup = @variable(m, lower_bound=0)#, upper_bound=1.2)
    
	if SetStart
            aux = (getstartvalue(v_k) - getstartvalue(v))/gen_band
	    JuMP.set_start_value(num, max(0, aux))
	    JuMP.set_start_value(nup, max(0,-aux))
	end

	# add exact PV/PQ coupling constraints to the model
	#@constraint(m, v_k + gen_band*nup - gen_band*num == v)
        @constraint(m, v_k + nup - num == v)
	if o.CouplingMode == :complementarity
		@NLconstraint(m, nup*(q_k - Qub) == 0)
		@NLconstraint(m, num*(q_k - Qlb) == 0)
	elseif o.CouplingMode == :complementarityapprox
            reg = o.SmoothingParamCoupling
            @NLconstraint(m, -reg<=num*(q_k/gen_band-Qlb/gen_band) <=reg)
            @NLconstraint(m, -reg<=nup*(q_k/gen_band-Qub/gen_band) <=reg)
            #@NLconstraint(m, nup * num <=reg)
 	elseif o.CouplingMode == :fischerburm
            reg = o.SmoothingParamCoupling
	    @NLconstraint(m, ( (( q_k - Qlb)/gen_band)^2 + num^2 + reg )^0.5 >= ( q_k - Qlb)/gen_band + num)
	    @NLconstraint(m, ( ((-q_k + Qub)/gen_band)^2 + nup^2 + reg )^0.5 >= (-q_k + Qub)/gen_band + nup)
            
            #reg = max(1e-8, reg^2)
            @constraint(m, nup*num<=reg)
        elseif o.CouplingMode == :fischerburmlifted
            reg = 1e-3*max(1e-2,abs(Qub))
            @constraint(m, nup*num<=reg/10)

            tu = @variable(m)
            @constraint(m, tu^2 <= (-q_k+Qub)^2 + nup^2 + reg)
            @constraint(m, tu   >=  -q_k+Qub    + nup   )
            
            reg = 1e-3*max(1e-2,abs(Qlb))
            tl = @variable(m)
            @constraint(m, tl^2 <= ( q_k-Qlb)^2 + num^2 + reg)
            @constraint(m, tl   >=   q_k-Qlb    + num   )

            if SetStart
                JuMP.set_start_value(tu, 0)
                JuMP.set_start_value(tl, 0)
            end
	elseif o.CouplingMode == :fischerburmsquare
		@NLconstraint(m, ( ((+q_k - Qlb)^2 + num^2 )^0.5 - (+q_k - Qlb + num) )^2 == 0)
		@NLconstraint(m, ( ((-q_k + Qub)^2 + nup^2 )^0.5 - (-q_k + Qub + nup) )^2 == 0)
	else
	    warn("unrecognized PVPQ coupling mode")
            return 0
	end
	return 1
	
end

## returns
## Gk    - indexes of all generators excepting 'outidx' if 'ConType' is :generator
## Gkp   - indexes of participating generators
## Gknop - indexes of non-participating generators
function getAGCParticipation(N, G, ConType, outidx, IndexSets)
    # parse IndexSets used by function
    L_Nidx = IndexSets[1]
    T_Nidx = IndexSets[2]
    G_Nidx = IndexSets[4]
    # determine contingency indices
    Garea = N[!,:Area][G_Nidx]
    Ak = Int[]
    Gk = collect(1:size(G, 1))
    if ConType == :Generator
	deleteat!(Gk, outidx)
	push!(Ak, N[!,:Area][G_Nidx[outidx]])
    elseif ConType == :Line
	append!(Ak, unique(sort(N[!,:Area][L_Nidx[outidx,:]])))
    elseif ConType == :Transformer
	append!(Ak, unique(sort(N[!,:Area][T_Nidx[outidx,:]])))
    else
	warn("unknown contingency type ", ConType)
    end
    Gkareaidx = indexin(Garea[Gk], Ak)
    
    #Gkp = Gk[Gkareaidx .!= nothing]
    #Gknop = Gk[Gkareaidx .== nothing]
    discriminant = .&(Gkareaidx .!= nothing, abs.(G[!,:alpha][Gk]) .> 1E-8)
        # found many generators with alpha = 0 on large instances
    Gkp = Gk[discriminant]
    Gknop = Gk[.!discriminant]
    
    return Gk, Gkp, Gknop
end

## function to add all SC AC OPF coupling constraints for a given contingency

function addcouplingcons!(m::JuMP.Model, o::GoOptions, v_n::AbstractVector,
	p_g::AbstractVector, v_nk::AbstractVector, p_gk::AbstractVector,
	q_gk::AbstractVector, delta_k::Union{Real, JuMP.VariableRef},
	N::DataFrame, G::DataFrame, ConType::Symbol, outidx::Int,
	IndexSets::Tuple; sclamp::Union{Function, Nothing}=nothing,
	sstep::Union{Function, Nothing}=nothing, SetStart::Bool=true)
	
    Gk, Gkp, Gknop = getAGCParticipation(N, G, ConType, outidx, IndexSets)

    # enforce AGC coupling constraints -- non-participant generators
    @constraint(m, [g=Gknop], p_gk[g] == p_g[g])
    
    # enforce AGC coupling constraints -- participating generators
    nAGCCons = 0
    for g = Gkp
	nAGCCons += addAGCconstraint!(m, o, p_g[g], delta_k, p_gk[g], G[!,:Plb][g], G[!,:Pub][g],
			              G[!,:alpha][g], sclamp=sclamp, SetStart=SetStart)
    end
    @printf("AGC:  %d gens participating: added %d constraints\n", size(Gkp,1), nAGCCons)
    
    #G_Nidx = IndexSets[4]
    # enforce PV/PQ coupling constraints
    #nPVPQCons = 0
    #for g = Gk
    #   nPVPQCons += addPVPQconstraint!(m, o, q_gk[g], v_n[G_Nidx[g]], v_nk[G_Nidx[g]], G[!,:Qlb][g],
    #                                   G[!,:Qub][g], sstep=sstep, SetStart=SetStart)
    #end

    G_Nidx = IndexSets[4]
    Gn = IndexSets[10]
    # enforce PV/PQ coupling bus-by-bus (avoids definition of unnecessary aux vars and cons)
    nPVPQGens = 0
    nPVPQCons = 0
    N_PVPQ = unique(sort(G_Nidx[Gk]))
    if ConType == :Generator
        gout = outidx
    else
        gout = nothing
    end
    for n = N_PVPQ
        numout = 0
        numfixed = 0
        gagg = Int[]
        Qagglb = 0
        Qaggub = 0
        for g = Gn[n]
            if g == gout
                numout += 1
                continue
            elseif G[!,:Qub][g]-G[!,:Qlb][g] < 1e-8
                numfixed += 1
                continue
            end
            push!(gagg, g)
            Qagglb += G[!,:Qlb][g]
            Qaggub += G[!,:Qub][g]
        end
	@assert length(gagg) + numout + numfixed == length(Gn[n])  # counting generators
        @assert numout == 0 || length(Gn[n]) > 1		   # this bus should have never been in N_PVPQ
        nPVPQGens += length(gagg) + numfixed
        if length(gagg) == 0
            continue
        end
        qagg = @NLexpression(m, sum(q_gk[g] for g=gagg))
        nPVPQCons += addPVPQconstraint!(m, o, qagg, v_n[n], v_nk[n], Qagglb, Qaggub,
                                        sstep=sstep, SetStart=SetStart)
    end
    @printf("PVPQ: %d gens participating: added %d constraints\n", nPVPQGens, nPVPQCons)
    
    # return to caller
    return nothing
    
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

## function to add piecewise linear penalty for a case (i.e. base or contingency)

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

## function to construct quadratic penalty terms for a case (i.e. base or contingency)

function addquadpenaltyblock!(expr::JuMP.GenericQuadExpr,
	slack::AbstractArray, a::Float64, b::Float64)::Nothing
	# a*sum(x^2) + b*sum(x)
	for i = 1:length(slack)
		add_to_expression!(expr, a, slack[i], slack[i])
	end
	for i = 1:length(slack)
		add_to_expression!(expr, b, slack[i])
	end
	return nothing
end

function quadpenaltyfunctions(
	pslackm_n::AbstractVector, pslackp_n::AbstractVector,
	qslackm_n::AbstractVector, qslackp_n::AbstractVector,
	sslack_li::AbstractArray, sslack_ti::AbstractArray,
	P::DataFrame)
	
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
	
	# return expression containing all penalty terms
	return casepenalty
	
end

## function for solving SC AC OPF
function SolveSCACOPF(o::GoOptions, N::DataFrame, L::DataFrame, T::DataFrame,
	              SSh::DataFrame, G::DataFrame, K::DataFrame, P::DataFrame, NLSolver; 
                      IndexSets::Union{Tuple, Nothing}=nothing,
                      StartingPoint::Union{SCACOPFState,Nothing}=nothing)
	# compute index and sets for performant formulation
        L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn, K_outidx =
            build_or_split_indexsets(IndexSets, N, L, T, SSh, G, K)

	RefBus = G_Nidx[argmax(G[!,:Pub])]	# bus with the biggest generator is reference
	NumK = size(K, 1)

        # (some) output for options used
        println("SolveSCACOPF with ", NumK, " contingencies. \nOptions: CouplingMode=", o.CouplingMode, 
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

            # use the above for all contingencies -> duplicate vectors/arrays for now
            v0_nk, theta0_nk, p0_lik, q0_lik, p0_tik, q0_tik, b0_sk, p0_gk, q0_gk, 
            pslack0m_nk, pslack0p_nk, qslack0m_nk, qslack0p_nk, sslack0_lik, sslack0_tik =
                InitialContingenciesSolutionsFromBase(N, L, T, SSh, G, K,
                                                      v0_n, theta0_n, p0_li, q0_li, p0_ti, q0_ti, b0_s, p0_g, q0_g, 
                                                      pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti,
                                                      IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
		                                                   SShn, Gn, K_outidx))
            delta0_k = zeros(NumK)
        else
            #! this is also temporary code
            @assert num_cont(StartingPoint)==NumK

            ## base case
            v0_n, theta0_n, b0_s, p0_g, q0_g = getBaseStates(StartingPoint)
            
            p0_li, q0_li, p0_ti, q0_ti, pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti = 
                BaseCaseSolution(v0_n, theta0_n, b0_s, p0_g, q0_g,
                                 N, L, T, SSh, G, 
                                 IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
		                              SShn, Gn, K_outidx))
            println("HOT Started")
            #
            # contingencies
            #
            #! (dirty) use this to quickly allocate 
            v0_nk, theta0_nk, p0_lik, q0_lik, p0_tik, q0_tik, b0_sk, p0_gk, q0_gk, 
            pslack0m_nk, pslack0p_nk, qslack0m_nk, qslack0p_nk, sslack0_lik, sslack0_tik =
                InitialContingenciesSolutionsFromBase(N, L, T, SSh, G, K,
                                                      v0_n, theta0_n, p0_li, q0_li, p0_ti, q0_ti, b0_s, p0_g, q0_g, 
                                                      pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti,
                                                      IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
	                                                           SShn, Gn, K_outidx))
            delta0_k = zeros(NumK)
            for k=1:NumK
                delta0_k[k], v0_nk[:,k], theta0_nk[:,k], b0_sk[:,k], p0_gk[:,k], q0_gk[:,k] = 
                    getContingencyStates(StartingPoint, K[!,:Contingency][k])
                
                p0_lik[:,:,k], q0_lik[:,:,k], p0_tik[:,:,k], q0_tik[:,:,k], 
                pslack0m_nk[:,k], pslack0p_nk[:,k], qslack0m_nk[:,k], qslack0p_nk[:,k], sslack0_lik[:,:,k], sslack0_tik[:,:,k] =
                    BaseCaseSolution(v0_nk[:,k], theta0_nk[:,k], b0_sk[:,k], p0_gk[:,k], q0_gk[:,k],
                                     N, L, T, SSh, G, (K[!,:ConType][k], K[!,:IDout][k]),
                                     IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
	                                          SShn, Gn, K_outidx))
                println("reusing starting point passed as argument for conting with index ", K[!,:Contingency][k])
            end
        end
	
	# create model
	m = Model(NLSolver)
	
	# pre-contingency variables
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

	# fix pre-contingency angle at reference bus to zero
	JuMP.fix(theta_n[RefBus], 0.0, force=true)
	
	# add power flow constraints for base case
	addpowerflowcons!(m, v_n, theta_n, p_li, q_li, p_ti, q_ti, b_s, p_g, q_g,
		pslackm_n, pslackp_n, qslackm_n, qslackp_n, sslack_li, sslack_ti,
		N, L, T, (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
		SShn, Gn, K_outidx), o.ThermalLimits)
	
	# generation reserves
	if true
		reservemargins = GenerationReserves(G, K)
		JuMP.set_upper_bound.(p_g, G[!,:Pub] - reservemargins)
	end
	
	# post-contingency variables
        @variable(m, N[!,:EVlb][n] <= v_nk[n=1:size(N, 1), k=1:NumK] <= N[!,:EVub][n], start = v0_nk[n,k])
        @variable(m, theta_nk[n=1:size(N, 1), k=1:NumK], start = theta0_nk[n,k])
	@variable(m, p_lik[l=1:size(L, 1), i=1:2, k=1:NumK], start = p0_lik[l,i,k])
	@variable(m, q_lik[l=1:size(L, 1), i=1:2, k=1:NumK], start = q0_lik[l,i,k])
	@variable(m, p_tik[t=1:size(T, 1), i=1:2, k=1:NumK], start = p0_tik[t,i,k])
	@variable(m, q_tik[t=1:size(T, 1), i=1:2, k=1:NumK], start = q0_tik[t,i,k])
	@variable(m, SSh[!,:Blb][s] <= b_sk[s=1:size(SSh, 1), k=1:NumK] <= SSh[!,:Bub][s], start = b0_sk[s,k])
	@variable(m, G[!,:Plb][g] <= p_gk[g=1:size(G, 1), k=1:NumK] <= G[!,:Pub][g], start = p0_gk[g,k])
	@variable(m, G[!,:Qlb][g] <= q_gk[g=1:size(G, 1), k=1:NumK] <= G[!,:Qub][g], start = q0_gk[g,k])
	@variable(m, delta_k[k=1:NumK], start = delta0_k[k])
	@variable(m, pslackm_nk[n=1:size(N, 1), k=1:NumK] >= 0, start = pslack0m_nk[n,k])
	@variable(m, pslackp_nk[n=1:size(N, 1), k=1:NumK] >= 0, start = pslack0p_nk[n,k])
	@variable(m, qslackm_nk[n=1:size(N, 1), k=1:NumK] >= 0, start = qslack0m_nk[n,k])
	@variable(m, qslackp_nk[n=1:size(N, 1), k=1:NumK] >= 0, start = qslack0p_nk[n,k])
	@variable(m, sslack_lik[l=1:size(L, 1), i=1:2, k=1:NumK] >= 0, start = sslack0_lik[l,i,k])
	@variable(m, sslack_tik[t=1:size(T, 1), i=1:2, k=1:NumK] >= 0, start = sslack0_tik[t,i,k])
	
	# post-contingency constraints
	for k = 1:NumK
	
		# unavailability
		if K[!,:ConType][k] == :Generator
			JuMP.fix(p_gk[K_outidx[k], k], 0.0, force = true)
			JuMP.fix(q_gk[K_outidx[k], k], 0.0, force = true)
			JuMP.set_start_value(p_gk[K_outidx[k], k], 0.0)
			JuMP.set_start_value(q_gk[K_outidx[k], k], 0.0)
		elseif K[!,:ConType][k] == :Line
			for i = 1:2
				JuMP.fix(p_lik[K_outidx[k], i, k], 0.0, force = true)
				JuMP.fix(q_lik[K_outidx[k], i, k], 0.0, force = true)
				JuMP.set_start_value(p_lik[K_outidx[k], i, k], 0.0)
				JuMP.set_start_value(q_lik[K_outidx[k], i, k], 0.0)
			end                        
		elseif K[!,:ConType][k] == :Transformer
			for i = 1:2
				JuMP.fix(p_tik[K_outidx[k], i, k], 0.0, force = true)
				JuMP.fix(q_tik[K_outidx[k], i, k], 0.0, force = true)
				JuMP.set_start_value(p_tik[K_outidx[k], i, k], 0.0)
				JuMP.set_start_value(q_tik[K_outidx[k], i, k], 0.0)
			end
		else
			warn("unrecognized contingency type ", K[!,:ConType][k],
				" found at contingency ", K[!,:Contingency][k])
		end
		
		# fix post-contingency angle at reference to zero
		JuMP.fix(theta_nk[RefBus, k], 0.0, force=true)
	        JuMP.set_start_value(theta_nk[RefBus, k], 0.0)

		# add power flow constraints for contingency
		addpowerflowcons!(m, view(v_nk, :, k), view(theta_nk, :, k),
			view(p_lik, :, :, k), view(q_lik, :, :, k),
			view(p_tik, :, :, k), view(q_tik, :, :, k),
			view(b_sk, :, k), view(p_gk, :, k), view(q_gk, :, k),
			view(pslackm_nk, :, k), view(pslackp_nk, :, k),
			view(qslackm_nk, :, k), view(qslackp_nk, :, k),
			view(sslack_lik, :, :, k), view(sslack_tik, :, :, k),
			N, L, T, (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin,
			Tidxn, Tin, SShn, Gn, K_outidx), o.ThermalLimits,
			SysCond=:Contingency, ConType=K[!,:ConType][k], outidx=K_outidx[k])
		
	end
	
	# register approximation functions
	if NumK > 0 && o.CouplingMode == :approx
		sclamp(x) = smoothclamp(x, o.SmoothingParamCoupling)
		JuMP.register(m, :sclamp, 1, sclamp,
			x -> SmoothApproximations.smoothclampprime(x, o.SmoothingParamCoupling),
			x -> SmoothApproximations.smoothclampprimeprime(x, o.SmoothingParamCoupling))
		sstep(x) = smoothstep(x, o.SmoothingParamCoupling)
		JuMP.register(m, :sstep, 1, sstep,
			x -> SmoothApproximations.smoothstepprime(x, o.SmoothingParamCoupling),
			x -> SmoothApproximations.smoothstepprimeprime(x, o.SmoothingParamCoupling))
	else
		sclamp = nothing
		sstep = nothing
	end
	
	## coupling constraints
	if o.CouplingMode != :ignore
		for k = 1:NumK
			addcouplingcons!(m, o, v_n, p_g, view(v_nk, :, k),
				view(p_gk, :, k), view(q_gk, :, k), delta_k[k],
				N, G, K[!,:ConType][k], K_outidx[k], 
				(L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin,
				Tidxn, Tin, SShn, Gn, K_outidx), sclamp=sclamp,
				sstep=sstep, SetStart=true)
		end
	end
	
	## objective function
	
	# production cost
	prodcost = zero(JuMP.AffExpr)
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
	basepenalty =
		if o.PenaltyType == :pcwslin
			addpenaltyfunctions!(m, pslackm_n, pslackp_n, qslackm_n, qslackp_n,
				sslack_li, sslack_ti, N, L, T, P, true)
		else
			quadpenaltyfunctions(pslackm_n, pslackp_n, qslackm_n, qslackp_n,
				sslack_li, sslack_ti, P)
		end
			

	# contingency penalties
	if o.PenaltyType == :pcwslin
		conpenalty = zero(JuMP.AffExpr)
		for k = 1:NumK
			add_to_expression!(conpenalty, addpenaltyfunctions!(m,
				view(pslackm_nk, :, k), view(pslackp_nk, :, k),
				view(qslackm_nk, :, k), view(qslackp_nk, :, k),
				view(sslack_lik, :, :, k), view(sslack_tik, :, :, k),
				N, L, T, P, true))
		end
	else
		conpenalty = zero(JuMP.QuadExpr)
		for k = 1:NumK
			add_to_expression!(conpenalty, quadpenaltyfunctions(
				view(pslackm_nk, :, k), view(pslackp_nk, :, k),
				view(qslackm_nk, :, k), view(qslackp_nk, :, k),
				view(sslack_lik, :, :, k), view(sslack_tik, :, :, k), P))
		end
	end
	
	# regularization term on delta
        if NumK>0 && o.CouplingMode == :ignore
            delta_pen = @expression(m, sum(delta_k[k]^2 for k=1:NumK))
        else
            delta_pen = zero(JuMP.AffExpr)
        end

	# declare objective
	@objective(m, Min, prodcost + DELTA*basepenalty + (1-DELTA)/max(1, NumK)*conpenalty + delta_pen)
	
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
	prodcostvalue = JuMP.value(prodcost)
	penaltybasevalue = DELTA*JuMP.value(basepenalty)
	if NumK > 0
		penaltyconvalue = (1-DELTA)/max(1, NumK)*JuMP.value(conpenalty)
	else
		penaltyconvalue = 0.0
	end
	objectivevalue = prodcostvalue + penaltybasevalue + penaltyconvalue
	println("Objective:           \$", round(objectivevalue, digits=1), "\n",
		"Generation cost:     \$", round(prodcostvalue, digits=1), "\n",
		"Penalty base:        \$", round(penaltybasevalue, digits=1), "\n",
		"Penalty contingency: \$", round(penaltyconvalue, digits=1))

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
	return SCACOPFState(prodcostvalue, penaltybasevalue, penaltyconvalue,
                            BaseState(JuMP.value.(v_n), 
                                      JuMP.value.(theta_n), 
                                      convert(Vector{Float64}, JuMP.value.(b_s)),
		                      JuMP.value.(p_g), 
                                      JuMP.value.(q_g), 0., 0.),
                            K[!,:Contingency], K[!,:IDout], K[!,:ConType],
                            convert(Array{Float64, 2}, JuMP.value.(v_nk)), 
                            convert(Array{Float64, 2}, JuMP.value.(theta_nk)), 
                            convert(Array{Float64, 2}, JuMP.value.(b_sk)),
		            convert(Array{Float64, 2}, JuMP.value.(p_gk)), 
                            convert(Array{Float64, 2}, JuMP.value.(q_gk)), 
                            convert(Array{Float64, 1}, JuMP.value.(delta_k)))
	
	#return objectivevalue,
	#	JuMP.value.(v_n), JuMP.value.(theta_n), convert(Vector{Float64}, JuMP.value.(b_s)),
	#	JuMP.value.(p_g), JuMP.value.(q_g),
	#	JuMP.value.(v_nk), JuMP.value.(theta_nk), convert(Array{Float64, 2}, JuMP.value.(b_sk)),
	#	JuMP.value.(p_gk), JuMP.value.(q_gk), JuMP.value.(delta_k)
end

## function for solving a contingency given:
##	- base case voltages
##	- base case dispatch
##	- contingency redispatch signal
##	- estimated contingency PV voltages
##
## Arguments notes
##  - K is a triple (:ConType, IDout, ConType)
function SolveContingencyFastBackend(o::GoOptions, N::DataFrame, L::DataFrame, T::DataFrame,
	                             SSh::DataFrame, G::DataFrame, K::Tuple{Symbol, Int, Int}, P::DataFrame,
                                     scopfstate::SCACOPFState, NLSolver;
	                             IndexSets::Union{Tuple, Nothing}=nothing,
                                     DeltaToUse::Union{Float64, Nothing}=nothing)

	println(" -- SolveContingencyFastBackend(base case and w. redispatch):\n",
		" Cont: ID ", K[3], " type: ", K[1], ", failed elem. index: ", K[2], "\n",
        	" Options: ThermalLimits=", o.ThermalLimits, 
                " SmoothingParam=", o.SmoothingParamCoupling)
        
        p_g = scopfstate.base_state.p_g
        v_n = scopfstate.base_state.v_n
        
        if haskey(scopfstate.cont_states, K[3])
            # use starting point given by the scacopf contingency solution
            cont_sol       = scopfstate.cont_states[K[3]]
            vapprox_nk     = cont_sol.v_n
            delta_k        = DeltaToUse==nothing ? cont_sol.delta : DeltaToUse
            thetaapprox_nk = cont_sol.theta_n
            bapprox_sk     = cont_sol.b_s
            qapprox_gk     = cont_sol.q_g
        else
            # use base case from scacopf solution
            vapprox_nk     = scopfstate.base_state.v_n
            @assert(DeltaToUse !== nothing)
            delta_k        = DeltaToUse
            thetaapprox_nk = scopfstate.base_state.theta_n
            bapprox_sk     = scopfstate.base_state.b_s
            qapprox_gk     = scopfstate.base_state.q_g
        end

	# compute index and sets for performant formulation
	if IndexSets == nothing
		L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn,
			K_outidx = indexsets(N, L, T, SSh, G, (K[1], K[2]))
	else
		L_Nidx = IndexSets[1]
		T_Nidx = IndexSets[2]
		SSh_Nidx = IndexSets[3]
		G_Nidx = IndexSets[4]
		Lidxn = IndexSets[5]
		Lin = IndexSets[6]
		Tidxn = IndexSets[7]
		Tin = IndexSets[8]
		SShn = IndexSets[9]
		Gn = IndexSets[10]
		K_outidx = IndexSets[11]
	end
	RefBus = G_Nidx[argmax(G[!,:Pub])]	# bus with the biggest generator is reference
	
	# compute contingency re-dispatch based on the redispatch signal
	p_gk = p_g + G[!,:alpha]*delta_k
	for g = 1:size(G, 1)
		if p_gk[g] < G[!,:Plb][g]
			p_gk[g] = G[!,:Plb][g]
		elseif p_gk[g] > G[!,:Pub][g]
			p_gk[g] = G[!,:Pub][g]
		end
	end
	if K[1] == :Generator
		p_gk[K_outidx] = 0.0
	end
	
	# compute AC OPF starting point
	theta0_nk = copy(thetaapprox_nk)
	theta0_nk .-= theta0_nk[RefBus]
	p0_lik, q0_lik, p0_tik, q0_tik, pslack0m_nk, pslack0p_nk, qslack0m_nk, qslack0p_nk,
		sslack0_lik, sslack0_tik = BaseCaseSolution(vapprox_nk, theta0_nk, bapprox_sk,
		p_gk, qapprox_gk, N, L, T, SSh, G, IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx,
		Lidxn, Lin, Tidxn, Tin, SShn, Gn, K_outidx))
	
	# create contingency model
	m = Model(NLSolver)
	
	# declare free variables
	@variable(m, theta_nk[n=1:size(N, 1)], start=theta0_nk[n])
	@variable(m, p_lik[l=1:size(L, 1), i=1:2], start=p0_lik[l,i])
	@variable(m, q_lik[l=1:size(L, 1), i=1:2], start=q0_lik[l,i])
	@variable(m, p_tik[t=1:size(T, 1), i=1:2], start=p0_tik[t,i])
	@variable(m, q_tik[t=1:size(T, 1), i=1:2], start=q0_tik[t,i])
	@variable(m, SSh[!,:Blb][s] <= b_sk[s=1:size(SSh, 1)] <= SSh[!,:Bub][s], start=bapprox_sk[s])
	@variable(m, pslackm_nk[n=1:size(N, 1)] >= 0, start=pslack0m_nk[n])
	@variable(m, pslackp_nk[n=1:size(N, 1)] >= 0, start=pslack0p_nk[n])
	@variable(m, qslackm_nk[n=1:size(N, 1)] >= 0, start=qslack0m_nk[n])
	@variable(m, qslackp_nk[n=1:size(N, 1)] >= 0, start=qslack0p_nk[n])
	@variable(m, sslack_lik[l=1:size(L, 1), i=1:2] >= 0, start=sslack0_lik[l,i])
	@variable(m, sslack_tik[t=1:size(T, 1), i=1:2] >= 0, start=sslack0_tik[t,i])
	
	# declare variables subject to PV/PQ switching
	v_nk = Vector{Any}(undef, size(N, 1))
	q_gk = Vector{Any}(undef, size(G, 1))
	if K[1] == :Generator
		ngout = G_Nidx[K_outidx]
	else
		ngout = nothing
	end
	for n=1:size(N, 1)
		if length(Gn[n]) == 0
			v_nk[n] = @variable(m, lower_bound=N[!,:EVlb][n],
				upper_bound=N[!,:EVub][n], start=vapprox_nk[n])
		elseif length(Gn[n]) == 1 && n == ngout
			v_nk[n] = @variable(m, lower_bound=N[!,:EVlb][n],
				upper_bound=N[!,:EVub][n], start=vapprox_nk[n])
			q_gk[K_outidx] = 0.0
		else
			if n == ngout
				q_gk[K_outidx] = 0.0
				Gnk = Gn[n][Gn[n] .!= K_outidx]
			else
				Gnk = Gn[n]
			end
			vdev, _ = map2step(vapprox_nk[n] - v_n[n], o.SmoothingParamCoupling)
			if vdev < -EPSILON
				v_nk[n] = @variable(m, lower_bound=N[!,:EVlb][n], upper_bound=v_n[n],
					start=.5*(N[!,:EVlb][n] + v_n[n]))
				for g=Gnk
					q_gk[g] = G[!,:Qub][g]
				end
			elseif vdev > EPSILON
				v_nk[n] = @variable(m, lower_bound=v_n[n], upper_bound=N[!,:EVub][n],
					start=.5*(v_n[n] + N[!,:EVub][n]))
				for g=Gnk
					q_gk[g] = G[!,:Qlb][g]
				end
			else
				v_nk[n] = v_n[n]
				for g=Gnk
					q_gk[g] = @variable(m, lower_bound=G[!,:Qlb][g],
						upper_bound=G[!,:Qub][g], start=.5*(G[!,:Qlb][g] + G[!,:Qub][g]))
				end
			end
		end
	end
	unassigned_voltages = collect(n for n=1:size(N, 1) if !isassigned(v_nk, n))
	unassigned_reactives = collect(g for g=1:size(G, 1) if !isassigned(q_gk, g))
	if length(unassigned_voltages) > 0 || length(unassigned_reactives) > 0
		warn("unassigned voltages at buses ", unassigned_voltages,
			" -- or -- unassigned reactive power at generators ", unassigned_reactives,
			" at buses ", G[!,:Bus][unassigned_reactives], " (G_Nidx[unassigned] = ",
			G_Nidx[unassigned_reactives], ")")
	end

	# branch unavailability
	if K[1] == :Line
		for i = 1:2
			JuMP.fix(p_lik[K_outidx, i], 0.0, force=true)
			JuMP.fix(q_lik[K_outidx, i], 0.0, force=true)
			JuMP.set_start_value(p_lik[K_outidx,i], 0.0)
			JuMP.set_start_value(q_lik[K_outidx,i], 0.0)
			JuMP.set_start_value(sslack_lik[K_outidx,i], 0.0)
		end
	elseif K[1] == :Transformer
		for i = 1:2
			JuMP.fix(p_tik[K_outidx, i], 0.0, force=true)
			JuMP.fix(q_tik[K_outidx, i], 0.0, force=true)
			JuMP.set_start_value(p_tik[K_outidx,i], 0.0)
			JuMP.set_start_value(q_tik[K_outidx,i], 0.0)
			JuMP.set_start_value(sslack_tik[K_outidx,i], 0.0)
		end
	end
	
	# post-contingency angle zero at reference
	JuMP.fix(theta_nk[RefBus], 0.0, force=true)
	
	# add power flow constraints for contingency
	addpowerflowcons!(m, v_nk, theta_nk, p_lik, q_lik, p_tik, q_tik, b_sk, p_gk, q_gk,
		pslackm_nk, pslackp_nk, qslackm_nk, qslackp_nk, sslack_lik, sslack_tik,
		N, L, T, (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
		SShn, Gn, K_outidx), o.ThermalLimits, SysCond=:Contingency,
		ConType=K[1], outidx=K_outidx)
	
	## objective function
	conpenalty =
		if o.PenaltyType == :pcwslin
			addpenaltyfunctions!(m, pslackm_nk, pslackp_nk, qslackm_nk, qslackp_nk,
				sslack_lik, sslack_tik, N, L, T, P, true)
		else
			quadpenaltyfunctions(pslackm_nk, pslackp_nk, qslackm_nk, qslackp_nk,
				sslack_lik, sslack_tik, P)
		end
	@objective(m, Min, conpenalty)
	
	## attempt to solve contingency residual problem
	JuMP.optimize!(m)
        solvestatusprimal = JuMP.primal_status(m)
        solvestatus = JuMP.termination_status(m)
        println("JuMP optimize MOI termination_status: ", solvestatus, " ", solvestatusprimal)
	if (solvestatusprimal != MOI.FEASIBLE_POINT) && (solvestatusprimal != MOI.NEARLY_FEASIBLE_POINT)
		error("solver failed to find a feasible solution")
	end
	
	## return solution
	v_nk_value = Vector{Float64}(undef, size(N, 1))
	q_gk_value = Vector{Float64}(undef, size(G, 1))
	for n=1:size(N, 1)
		if typeof(v_nk[n]) == Float64
			v_nk_value[n] = v_nk[n]
		else
			v_nk_value[n] = JuMP.value(v_nk[n])
		end
	end
	for g=1:size(G, 1)
		if typeof(q_gk[g]) == Float64
			q_gk_value[g] = q_gk[g]
		else
			q_gk_value[g] = JuMP.value(q_gk[g])
		end
	end

        solution = ContingencyState(K[3], K[2], K[1], delta_k, v_nk_value,
                                    JuMP.value.(theta_nk), convert(Vector{Float64}, JuMP.value.(b_sk)), 
                                    p_gk, q_gk_value,
                                    JuMP.value(conpenalty), JuMP.objective_value(m))
        return solution
	#return v_nk_value, JuMP.value.(theta_nk), convert(Vector{Float64}, JuMP.value.(b_sk)),
	#	p_gk, q_gk_value
end


## Arguments notes
##  - K is a triple (:ConType, IDout, ConType)
function SolveContingencyWithFixing(o::GoOptions, N::DataFrame, L::DataFrame, T::DataFrame,
	                            SSh::DataFrame, G::DataFrame, K::Tuple{Symbol, Int, Int}, P::DataFrame,
                                    scopfstate::SCACOPFState, NLSolver;
	                            StartingPoint::ContingencyState,
                                    IndexSets::Union{Tuple, Nothing}=nothing,
                                    fixVolt::Bool=false)

    println("\n -- SolveContingencyWithFixing:\n",
	    " Cont: ID ", K[3], " type: ", K[1], ", failed elem. index: ", K[2], "\n",
            " Options: CouplingMode=", o.CouplingMode, " ThermalLimits=", o.ThermalLimits, 
            " SmoothingParam=", o.SmoothingParamCoupling)
    
    p_g = scopfstate.base_state.p_g
    v_n = scopfstate.base_state.v_n

    # this is the solution with approximate coupling found by SolveContingencyBackend
    cont_sol = StartingPoint

    vapprox_nk     = cont_sol.v_n
    deltaapprox_k  = cont_sol.delta
    thetaapprox_nk = cont_sol.theta_n
    bapprox_sk     = cont_sol.b_s
    qapprox_gk     = cont_sol.q_g
    papprox_gk     = cont_sol.p_g

    # compute index and sets for performant formulation
    if IndexSets == nothing
	L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn,
	K_outidx = indexsets(N, L, T, SSh, G, (K[1], K[2]))
    else
	L_Nidx = IndexSets[1]; T_Nidx = IndexSets[2]; SSh_Nidx = IndexSets[3]; G_Nidx = IndexSets[4]
	Lidxn = IndexSets[5]; Lin = IndexSets[6]; Tidxn = IndexSets[7]; Tin = IndexSets[8]
	SShn = IndexSets[9]; Gn = IndexSets[10]; K_outidx = IndexSets[11]
    end
    RefBus = G_Nidx[argmax(G[!,:Pub])]	# bus with the biggest generator is reference

    @assert size(K_outidx,1)==1


    
    # compute AC OPF starting point
    theta0_nk = copy(thetaapprox_nk)
    theta0_nk .-= theta0_nk[RefBus]
    p0_lik, q0_lik, p0_tik, q0_tik, pslack0m_nk, pslack0p_nk, qslack0m_nk, qslack0p_nk, sslack0_lik, sslack0_tik = 
        BaseCaseSolution(vapprox_nk, theta0_nk, bapprox_sk, papprox_gk, qapprox_gk, 
                         N, L, T, SSh, G, 
                         IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx,
		                      Lidxn, Lin, Tidxn, Tin, SShn, Gn, K_outidx))
    
    # create contingency model
    m = Model(NLSolver)
	
    # declare free variables
    @variable(m, theta_nk[n=1:size(N, 1)], start=theta0_nk[n])
    @variable(m, p_lik[l=1:size(L, 1), i=1:2], start=p0_lik[l,i])
    @variable(m, q_lik[l=1:size(L, 1), i=1:2], start=q0_lik[l,i])
    @variable(m, p_tik[t=1:size(T, 1), i=1:2], start=p0_tik[t,i])
    @variable(m, q_tik[t=1:size(T, 1), i=1:2], start=q0_tik[t,i])
    @variable(m, SSh[!,:Blb][s] <= b_sk[s=1:size(SSh, 1)] <= SSh[!,:Bub][s], start=bapprox_sk[s])
    @variable(m, pslackm_nk[n=1:size(N, 1)] >= 0, start=pslack0m_nk[n])
    @variable(m, pslackp_nk[n=1:size(N, 1)] >= 0, start=pslack0p_nk[n])
    @variable(m, qslackm_nk[n=1:size(N, 1)] >= 0, start=qslack0m_nk[n])
    @variable(m, qslackp_nk[n=1:size(N, 1)] >= 0, start=qslack0p_nk[n])
    @variable(m, sslack_lik[l=1:size(L, 1), i=1:2] >= 0, start=sslack0_lik[l,i])
    @variable(m, sslack_tik[t=1:size(T, 1), i=1:2] >= 0, start=sslack0_tik[t,i])


    # declare variables subject to AGC
    Gk, Gkp, Gknop = getAGCParticipation(N, G,  K[1], K_outidx, IndexSets)

    p_gk = Vector{Any}(undef, size(G, 1))
    for g = Gknop
        p_gk[g] = p_g[g]
    end

    if K[1] == :Generator
        p_gk[K_outidx] = 0.0
    end

    delta_lb = -1e+20
    delta_ub =  1e+20
    rtol  = sqrt(o.SmoothingParamCoupling)
    rtolp = sqrt(o.SmoothingParamCoupling)/10.
    for g = Gkp
        @assert g!=K_outidx || K[1] != :Generator
        gen_band = G[!,:Pub][g] - G[!,:Plb][g]

        #@assert gen_band>=1e-8 "fixed generators [1]"
        if gen_band < 1e-8
            p_gk[g] = 0.5*(G[!,:Pub][g] + G[!,:Plb][g])
            #@printf("p_g=%g Plb=%g Pub=%g p_gk=%g   g=%d\n", p_g[g], G[!,:Plb][g], G[!,:Pub][g], p_gk[g], g)
            if abs(p_gk[g] - p_g[g])>=1e-8
                @warn("fixed p_gk, but p_g does not match, generator ", g)
            end
            continue
        end

        dist_lower = (papprox_gk[g]-G[!,:Plb][g])/gen_band
        dist_upper = (G[!,:Pub][g]-papprox_gk[g])/gen_band
        if dist_lower<0; warn("distlower=", dist_lower); end
        if dist_upper<0; warn("distupper=", dist_upper); end

        distp = (papprox_gk[g] - p_g[g] - G[!,:alpha][g]*deltaapprox_k)/gen_band

        if dist_lower > rtol && dist_upper > rtol
            #p_gk[g] = @variable(m, lower_bound=G[!,:Plb][g], upper_bound=G[!,:Pub][g], start=papprox_gk[g])
            #delta_k = @variable(m)
            #@constraint(m, p_gk[g]   == p_g[g] + G[!,:alpha][g]*delta_k)
        elseif dist_lower <= rtol
            if distp > rtolp
                #strict complementarity, fixing is clear
                #p_gk[g] = G[!,:Plb][g]
                delta_ub = min(delta_ub, (papprox_gk[g] - p_g[g])/G[!,:alpha][g])
            else
                #degenerate complementarity, enforce equality
                #p_gk[g] = @variable(m, lower_bound=G[!,:Plb][g], upper_bound=G[!,:Pub][g], start=papprox_gk[g])
                #delta_k = @variable(m)
                #@constraint(m, p_gk[g]   == p_g[g] + G[!,:alpha][g]*delta_k)
            end
        else #dist_upper <= rtol
            if distp < -rtolp
                #strict complementarity, fixing is clear
                #p_gk[g] = G[!,:Pub][g]
                delta_lb = max(delta_lb, (papprox_gk[g] - p_g[g])/G[!,:alpha][g])
            else
                #degenerate complementarity, enforce equality
                #p_gk[g] = @variable(m, lower_bound=G[!,:Plb][g], upper_bound=G[!,:Pub][g], start=papprox_gk[g])
                #delta_k = @variable(m)
                #@constraint(m, p_gk[g]   == p_g[g] + G[!,:alpha][g]*delta_k)
            end
        end
    end

    deltaIsAVar = true #always enforced to be var #(delta_lb !=-1e+20) || (delta_ub!= 1e+20)
    deltaIsFeas = true
    if deltaIsAVar
        if delta_ub-delta_lb>=1e-16
            @printf("Fixing resulted in feasible bounds for delta: (%12.4e, %12.4e)\n", delta_lb, delta_ub)
        else
            @printf("Fixing resulted in INFEASIBLE bounds for delta: (%12.4e, %12.4e)\n", delta_lb, delta_ub)
            deltaIsFeas = false
        end
    else
        @printf("delta will not be posed as a variable\n")
    end
    
   if deltaIsAVar
       delta_k = @variable(m)
       if deltaIsFeas
           if delta_ub< 1e+20; JuMP.set_upper_bound(delta_k, delta_ub); end
           if delta_lb>-1e+20; JuMP.set_lower_bound(delta_k, delta_lb); end
       else
           #no bounds on delta
       end
   else
       @assert deltaIsFeas
       delta_k = deltaapprox_k
   end

   ## second pass over Gkp
   for g = Gkp
        gen_band = G[!,:Pub][g] - G[!,:Plb][g]

        if gen_band < 1e-8
            p_gk[g] = 0.5*(G[!,:Pub][g] + G[!,:Plb][g])
            if abs(p_gk[g] - p_g[g])>=1e-8
                @warn("fixed p_gk, but p_g does not match, generator ", g)
            end
            continue
        end

        dist_lower = (papprox_gk[g]-G[!,:Plb][g])/gen_band
        dist_upper = (G[!,:Pub][g]-papprox_gk[g])/gen_band

        if dist_lower<0; warn("distlower=", dist_lower); end
        if dist_upper<0; warn("distupper=", dist_upper); end

        distp = (papprox_gk[g] - p_g[g] - G[!,:alpha][g]*deltaapprox_k)/gen_band

        #@printf("AGC:fixing p_gk=%12.4e p_g=%12.4e | p_gk-p_g-a*d=%12.4e  p_gk-plb=%12.4e pub-p_gk=%12.4e | delta_k=%12.4e Plb=%12.4e Pub=%12.4e",
        #        papprox_gk[g], p_g[g], papprox_gk[g] - p_g[g] - G[!,:alpha][g]*deltaapprox_k,
        #        papprox_gk[g]-G[!,:Plb][g], G[!,:Pub][g]-papprox_gk[g], deltaapprox_k,
        #        G[!,:Plb][g], G[!,:Pub][g])

        if dist_lower > rtol && dist_upper > rtol
            p_gk[g] = @variable(m, lower_bound=G[!,:Plb][g], upper_bound=G[!,:Pub][g], start=papprox_gk[g])
            @constraint(m, p_gk[g]   == p_g[g] + G[!,:alpha][g]*delta_k)
            #@printf("  inside, enf. eq\n")
        elseif dist_lower <= rtol
            if distp > rtolp && deltaIsFeas
                #strict complementarity and delta is feasible -> fixing is clear
                p_gk[g] = G[!,:Plb][g]
                #@printf("  fixed p_gk at lower\n")
            else
                #degenerate complementarity or delta is not feasible , enforce equality
                p_gk[g] = @variable(m, lower_bound=G[!,:Plb][g], upper_bound=G[!,:Pub][g], start=papprox_gk[g])
                @constraint(m, p_gk[g]   == p_g[g] + G[!,:alpha][g]*delta_k)
                #@printf("  p_gk at lower but enforcing eq.; deltaIsFeas=%d; degen compl\n",deltaIsFeas)
            end
        else #dist_upper <= rtol
            if distp < -rtolp && deltaIsFeas
                #strict complementarity and delta is feasible, fixing is clear
                p_gk[g] = G[!,:Pub][g]
                #@printf("  fixed p_gk at upper\n")
            else
                #degenerate complementarity or delta is infeasible, enforce equality
                p_gk[g] = @variable(m, lower_bound=G[!,:Plb][g], upper_bound=G[!,:Pub][g], start=papprox_gk[g])
                @constraint(m, p_gk[g]   == p_g[g] + G[!,:alpha][g]*delta_k)
                #@printf("  p_gk at upper but enforcing eq.; deltaIsFeas=%d; degen compl\n",deltaIsFeas)
            end
        end
    end


    #
    ## declare variables subject to PV/PQ switching
    #
    v_nk = Vector{Any}(undef, size(N, 1))
    q_gk = Vector{Any}(undef, size(G, 1))
    if K[1] == :Generator
	ngout = G_Nidx[K_outidx]
    else
	ngout = nothing
    end

    # whether to fix voltage or q in case of degenerate complementarity
    #fixVolt = false

    for n=1:size(N, 1)
	if length(Gn[n]) == 0
	    v_nk[n] = @variable(m, lower_bound=N[!,:EVlb][n], upper_bound=N[!,:EVub][n], start=vapprox_nk[n])
	elseif length(Gn[n]) == 1 && n == ngout
	    v_nk[n] = @variable(m, lower_bound=N[!,:EVlb][n],
				upper_bound=N[!,:EVub][n], start=vapprox_nk[n])
	    q_gk[K_outidx] = 0.0
	else
# #	    if n == ngout
# #		q_gk[K_outidx] = 0.0
# #		Gnk = Gn[n][Gn[n] .!= K_outidx]
# #	    else
# #		Gnk = Gn[n]
# #	    end
	    Gnk = Int[]
	    if n == ngout
		for g = Gn[n]
		    if g == K_outidx
			q_gk[g] = 0.0
		    elseif G[!,:Qub][g]-G[!,:Qlb][g] <= 1e-8
			q_gk[g] = .5*(G[!,:Qub][g]+G[!,:Qlb][g])
		    else
			push!(Gnk, g)
		    end
		end
	    else
		for g = Gn[n]
		    if G[!,:Qub][g]-G[!,:Qlb][g] <= 1e-8
			q_gk[g] = .5*(G[!,:Qub][g]+G[!,:Qlb][g])
		    else
			push!(Gnk, g)
		    end
		end
	    end
	    if length(Gnk) == 0
		v_nk[n] = @variable(m, lower_bound=N[!,:EVlb][n],
				    upper_bound=N[!,:EVub][n], start=vapprox_nk[n])
		continue
	    end


            #@printf("PVPQ: n=%d v_nk=%12.4e v_n=%12.4e reldiff=%12.4e   |   ",
            #        n, vapprox_nk[n], v_n[n], (vapprox_nk[n] - v_n[n]))
            #for g=Gnk
            #    @printf("q_gk=%12.4e q_gk-Qlb=%12.4e Qub-q_gk=%12.4e |  Qlb=%11.4e Qub=%11.4e g=%d",
            #            qapprox_gk[g], -G[!,:Qlb][g]+qapprox_gk[g], G[!,:Qub][g]-qapprox_gk[g], G[!,:Qlb][g], G[!,:Qub][g], g)
            #end

            vdev = (vapprox_nk[n] - v_n[n]) / max(1,abs(v_n[n]))
  	    Qlbn = sum(G[!,:Qlb][Gnk])
	    Qubn = sum(G[!,:Qub][Gnk])
	    gen_band = Qubn - Qlbn
	    qapprox_nk = sum(qapprox_gk[Gnk])
	    dist_lower = (qapprox_nk - Qlbn)/gen_band
	    dist_upper = (Qubn - qapprox_nk)/gen_band
          
            # # ngen = size(Gnk,1)
            # # dist_lower = zeros(ngen); dist_upper = zeros(ngen)
            # # for ig=1:ngen
            # #     g=Gnk[ig]
            # #     gen_band=G[!,:Qub][g]-G[!,:Qlb][g]
            # #     if gen_band<1e-8
            # #         if abs(qapprox_gk[g]-G[!,:Qlb][g])>1e-8 || abs(G[!,:Qub][g]-qapprox_gk[g])>1e-8
            # #             msg = @sprintf("Qlb==Qub is fixed but q_gk in is not close to them. bus=%d gen=%d", n, g) 
            # #             @warn msg
            # #         end
            # #     else
            # #         dist_lower[ig] = max(0, (-G[!,:Qlb][g]+qapprox_gk[g])/gen_band)
            # #         dist_upper[ig] = max(0, ( G[!,:Qub][g]-qapprox_gk[g])/gen_band)
            # #     end
            # # end
            # # min_dist_lower =  minimum(dist_lower); min_dist_upper = minimum(dist_upper)

            ##rtol = 5e-2; rtolv=1e-2
            rtol = 1e-2; rtolv=1e-3
            if dist_lower > rtol && dist_upper > rtol
                #inside -> fix v_nk
                ##@assert abs(vdev) < 1e-2
                if abs(vdev)>=sqrt(o.SmoothingParamCoupling)
                    msg = @sprintf("sum(q_gk)q_gk in is inside the bounds (dist_lower,dist_upper)=(%g,%g), but volt dev is large %g. %d gens at the bus %d.", 
                                   dist_lower, dist_upper, vdev, length(Gnk), n)
                    @warn(msg)
                    @printf("PVPQ: v_nk=%12.4e v_n=%12.4e reldiff=%12.4e   |   ",
                            vapprox_nk[n], v_n[n], (vapprox_nk[n] - v_n[n]))
                    for g=Gnk
                        @printf("q_gk=%12.4e q_gk-Qlb=%12.4e Qub-q_gk=%12.4e |  Qlb=%11.4e Qub=%11.4e ",
                                qapprox_gk[g], -G[!,:Qlb][g]+qapprox_gk[g], G[!,:Qub][g]-qapprox_gk[g], G[!,:Qlb][g], G[!,:Qub][g])
                    end
                    println("")
                end
                v_nk[n] = v_n[n]
                for g=Gnk
		    q_gk[g] = @variable(m, lower_bound=G[!,:Qlb][g], upper_bound=G[!,:Qub][g], start=qapprox_gk[g])
		end
                #@printf("  fixing v_nk to v_n\n")
# #            elseif min_dist_lower==0 && min_dist_upper==0
# #                #at least one of the generators is fixed -> fix voltage
# #                #v_nk[n] = v_n[n]
# #                v_nk[n] = @variable(m, lower_bound=N[!,:EVlb][n], upper_bound=N[!,:EVub][n], start=vapprox_nk[n])
# #               for g=Gnk
# #		    q_gk[g] = @variable(m, lower_bound=G[!,:Qlb][g], upper_bound=G[!,:Qub][g], start=qapprox_gk[g])
# #               end
            elseif dist_lower <= rtol 
                #@assert vdev>=0
                if vdev >= rtolv 
                    ##strict complementarity -> fix q_gk     to Qlb               

                    #@printf("  fixing q_gk to Qlb\n")
                    v_nk[n] = @variable(m, lower_bound=v_n[n], upper_bound=N[!,:EVub][n], start=vapprox_nk[n])
                    for g=Gnk
		        q_gk[g] = G[!,:Qlb][g]
                    end
                else
                    #degenerate complementarity 

                    if fixVolt
                        #@printf("  fixing v_nk to v_n  (degen compl)\n")
                        v_nk[n] = v_n[n]
                        for g=Gnk
		            q_gk[g] = @variable(m, lower_bound=G[!,:Qlb][g], upper_bound=G[!,:Qub][g], start=qapprox_gk[g])
                        end
                        
                    else
                        #@printf("  fixing q_gk to Qlb  (degen compl)\n")
                        v_nk[n] = @variable(m, lower_bound=v_n[n], upper_bound=N[!,:EVub][n], start=vapprox_nk[n])
                        
                        for g=Gnk
		            q_gk[g] = G[!,:Qlb][g]
                        end
                    end
		end
            else #if dist_upper <= rtol
                #@assert vdev<=0
                if vdev <= - rtolv 
                    #strict complementarity -> fix q_gk to Qub
                    #@printf("  fixing q_gk to Qub\n")
                    v_nk[n] = @variable(m, lower_bound=N[!,:EVlb][n], upper_bound=v_n[n], start=vapprox_nk[n])
                    for g=Gnk
		        q_gk[g] = G[!,:Qub][g]
                    end
		else
                    #degenerate complementarity 

                    if fixVolt
                        #@printf("  fixing v_nk to v_n  (degen compl)\n")
                        v_nk[n] = v_n[n]
                        for g=Gnk
		            q_gk[g] = @variable(m, lower_bound=G[!,:Qlb][g], upper_bound=G[!,:Qub][g], start=qapprox_gk[g])
                        end
                    else
                        #@printf("  fixing q_gk to Qub  (degen compl)\n")
                        v_nk[n] = @variable(m, lower_bound=N[!,:EVlb][n], upper_bound=v_n[n], start=vapprox_nk[n])
                        for g=Gnk
		            q_gk[g] = G[!,:Qub][g]
                        end                    
                    end
                end                    
            end

	end
    end
    unassigned_voltages = collect(n for n=1:size(N, 1) if !isassigned(v_nk, n))
    unassigned_reactives = collect(g for g=1:size(G, 1) if !isassigned(q_gk, g))
    if length(unassigned_voltages) > 0 || length(unassigned_reactives) > 0
	warn("unassigned voltages at buses ", unassigned_voltages,
	      " -- or -- unassigned reactive power at generators ", unassigned_reactives,
	      " at buses ", G[!,:Bus][unassigned_reactives], " (G_Nidx[unassigned] = ",
	      G_Nidx[unassigned_reactives], ")")
    end

    # branch unavailability
    if K[1] == :Generator
        @assert typeof(q_gk[K_outidx]) == Float64
        q_gk[K_outidx] = 0.0
    elseif K[1] == :Line
        for i = 1:2
	    JuMP.fix(p_lik[K_outidx, i], 0.0, force=true)
	    JuMP.fix(q_lik[K_outidx, i], 0.0, force=true)
	    JuMP.set_start_value(p_lik[K_outidx,i], 0.0)
	    JuMP.set_start_value(q_lik[K_outidx,i], 0.0)
	    JuMP.set_start_value(sslack_lik[K_outidx,i], 0.0)
        end
    elseif K[1] == :Transformer
	for i = 1:2
	    JuMP.fix(p_tik[K_outidx, i], 0.0, force=true)
	    JuMP.fix(q_tik[K_outidx, i], 0.0, force=true)
	    JuMP.set_start_value(p_tik[K_outidx,i], 0.0)
	    JuMP.set_start_value(q_tik[K_outidx,i], 0.0)
	    JuMP.set_start_value(sslack_tik[K_outidx,i], 0.0)
	end
    end

    # post-contingency angle zero at reference
    JuMP.fix(theta_nk[RefBus], 0.0, force=true)
	
    # add power flow constraints for contingency
    addpowerflowcons!(m, v_nk, theta_nk, p_lik, q_lik, p_tik, q_tik, b_sk, p_gk, q_gk,
 		      pslackm_nk, pslackp_nk, qslackm_nk, qslackp_nk, sslack_lik, sslack_tik,
		      N, L, T, (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
		                SShn, Gn, K_outidx), o.ThermalLimits, SysCond=:Contingency,
		      ConType=K[1], outidx=K_outidx)
	
     ## objective function
     conpenalty =
		if o.PenaltyType == :pcwslin
			addpenaltyfunctions!(m, pslackm_nk, pslackp_nk, qslackm_nk, qslackp_nk,
				sslack_lik, sslack_tik, N, L, T, P, true)
		else
			quadpenaltyfunctions(pslackm_nk, pslackp_nk, qslackm_nk, qslackp_nk,
				sslack_lik, sslack_tik, P)
		end
     @objective(m, Min, conpenalty)
	
     ## attempt to solve contingency residual problem
     JuMP.optimize!(m)
     solvestatusprimal = JuMP.primal_status(m)
     solvestatus = JuMP.termination_status(m)
     println("JuMP optimize MOI termination_status: ", solvestatus, " ", solvestatusprimal)
     if (solvestatusprimal != MOI.FEASIBLE_POINT) && (solvestatusprimal != MOI.NEARLY_FEASIBLE_POINT)
	 error("solver failed to find a feasible solution")
     end
	
	## return solution
	v_nk_value = Vector{Float64}(undef, size(N, 1))
	q_gk_value = Vector{Float64}(undef, size(G, 1))
	p_gk_value = Vector{Float64}(undef, size(G, 1))
	for n=1:size(N, 1)
		if typeof(v_nk[n]) == Float64
			v_nk_value[n] = v_nk[n]
		else
			v_nk_value[n] = JuMP.value(v_nk[n])
		end
	end
	for g=1:size(G, 1)
	    if typeof(q_gk[g]) == Float64
		q_gk_value[g] = q_gk[g]
	    else
		q_gk_value[g] = JuMP.value(q_gk[g])
	    end
            p_gk_value[g] = typeof(p_gk[g]) == Float64 ? p_gk[g] : JuMP.value(p_gk[g])
	end

        delta_k_value = typeof(delta_k) == Float64 ? delta_k : JuMP.value(delta_k)

#for n=1:size(N, 1)
#    if abs(vapprox_nk[n] - v_nk_value[n]) >1e-3
#       @printf("n=%d voltage difference %g\n", n, abs(vapprox_nk[n] - v_nk_value[n]))
#    end
#end
#for g=1:size(G, 1)
#    if abs(papprox_gk[g]-p_gk_value[g])>1e-3 
#        @printf("g=%d pg difference %g\n", g, abs(papprox_gk[g]-p_gk_value[g]))
#    end
#end

#for g=1:size(G, 1)
#    if abs(qapprox_gk[g]-q_gk_value[g])>1e-3 
#        @printf("g=%d Qg difference %g\n", g, abs(qapprox_gk[g]-q_gk_value[g]))
#    end
#end

#println("!!!!1  ", argmax(JuMP.value.(pslackm_nk)), " -> ", maximum(JuMP.value.(pslackm_nk)), " ",  minimum(JuMP.value.(pslackm_nk))),
#println("!!!!2  ", argmax(JuMP.value.(pslackp_nk)), " -> ", maximum(JuMP.value.(pslackp_nk)), " ",  minimum(JuMP.value.(pslackp_nk)));
#println("!!!!3  ", argmax(JuMP.value.(qslackm_nk)), " -> ", maximum(JuMP.value.(qslackm_nk)), " ",  minimum(JuMP.value.(qslackm_nk)))
#println("!!!!4  ", argmax(JuMP.value.(qslackp_nk)), " -> ", maximum(JuMP.value.(qslackp_nk)), " ",  minimum(JuMP.value.(qslackp_nk)))			
#println("!!!!5  ", argmax(JuMP.value.(sslack_lik)), " -> ", maximum(JuMP.value.(sslack_lik)), " ",  minimum(JuMP.value.(sslack_lik)))
#println("!!!!6  ", argmax(JuMP.value.(sslack_tik)), " -> ", maximum(JuMP.value.(sslack_tik)), " ",  minimum(JuMP.value.(sslack_tik)))

        @printf("Fixing: ContingencyID=%d objs: conpenalty=%11.4e jumpobj=%11.4e  new delta=%11.3e\n", 
                K[3], JuMP.value(conpenalty), JuMP.objective_value(m), delta_k_value); flush(stdout)

        solution = ContingencyState(K[3], K[2], K[1], delta_k_value, v_nk_value,
                                    JuMP.value.(theta_nk), convert(Vector{Float64}, JuMP.value.(b_sk)), 
                                    p_gk_value, q_gk_value,
                                    JuMP.value(conpenalty), JuMP.objective_value(m))
        return solution

end


## backend function for solving a contingency using optimization:
## given are
##	- base case voltages
##	- base case dispatch

function SolveContingencyBackend(o::GoOptions, N::DataFrame, L::DataFrame, T::DataFrame,
	                         SSh::DataFrame, G::DataFrame, K::Tuple{Symbol, Int, Int}, P::DataFrame,
	                         scopfstate::SCACOPFState, NLSolver; 
                                 StartingPoint::Union{ContingencyState, Nothing}=nothing,
	                         IndexSets::Union{Tuple, Nothing}=nothing)
	
        println(" -- SolveContingencyBackend: ", 
                "Contingency ID/IDOut/Type=", K[3], "/", K[2], "/", K[1], "\n", 
                "Options: CouplingMode=", o.CouplingMode, " ThermalLimits=", o.ThermalLimits,
                " SmoothingParam=", o.SmoothingParamCoupling)
        #shortcuts
        v_n = scopfstate.base_state.v_n
        p_g = scopfstate.base_state.p_g


	# compute index and sets for performant formulation
	if IndexSets == nothing
		L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn,
			K_outidx = indexsets(N, L, T, SSh, G, (K[1], K[2]))
	else
		L_Nidx = IndexSets[1]
		T_Nidx = IndexSets[2]
		SSh_Nidx = IndexSets[3]
		G_Nidx = IndexSets[4]
		Lidxn = IndexSets[5]
		Lin = IndexSets[6]
		Tidxn = IndexSets[7]
		Tin = IndexSets[8]
		SShn = IndexSets[9]
		Gn = IndexSets[10]
		K_outidx = IndexSets[11]
	end
	RefBus = G_Nidx[argmax(G[!,:Pub])]	# bus with the biggest generator is reference
	
	# compute AC OPF starting point
	if StartingPoint == nothing
            
		#v0_n, theta0_n, p0_li, q0_li, p0_ti, q0_ti, b0_s, p0_g, q0_g, pslack0m_n, pslack0p_n,
		#	qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti = InitialBaseCaseSolution(
		#	N, L, T, SSh, G, IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
		#	SShn, Gn, K_outidx))
		#theta0_n .-= theta0_n[RefBus]
		#delta0 = 0
                
            delta0 = 0
            v0_n     = scopfstate.base_state.v_n     
	    theta0_n = scopfstate.base_state.theta_n 
	    b0_s     = scopfstate.base_state.b_s     
	    p0_g     = scopfstate.base_state.p_g     
	    q0_g     = scopfstate.base_state.q_g     
	    p0_li, q0_li, p0_ti, q0_ti, pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti = 
                BaseCaseSolution(v0_n, theta0_n, b0_s, p0_g, q0_g,
			         N, L, T, SSh, G, (K[1],K[2]), 
                                 IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin,
			                      Tidxn, Tin, SShn, Gn, K_outidx))
	else
		delta0   = StartingPoint.delta   #StartingPoint[1]
		v0_n     = StartingPoint.v_n     #StartingPoint[2]
		theta0_n = StartingPoint.theta_n #StartingPoint[3]
		b0_s     = StartingPoint.b_s     #StartingPoint[4]
		p0_g     = StartingPoint.p_g     #StartingPoint[5]
		q0_g     = StartingPoint.q_g     #StartingPoint[6]
		p0_li, q0_li, p0_ti, q0_ti, pslack0m_n, pslack0p_n, qslack0m_n, qslack0p_n, sslack0_li, sslack0_ti = 
                    BaseCaseSolution(v0_n, theta0_n, b0_s, p0_g, q0_g,
			             N, L, T, SSh, G, (K[1],K[2]), 
                                     IndexSets = (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin,
			                          Tidxn, Tin, SShn, Gn, K_outidx))
	end
	
	# create model
	m = Model(NLSolver)
	
	# contingency OPF variables
	@variable(m, N[!,:Vlb][n] <= v_nk[n=1:size(N, 1)] <= N[!,:Vub][n], start = v0_n[n])
	@variable(m, theta_nk[n=1:size(N, 1)], start = theta0_n[n])
	@variable(m, p_lik[l=1:size(L, 1), i=1:2], start = p0_li[l,i])
	@variable(m, q_lik[l=1:size(L, 1), i=1:2], start = q0_li[l,i])
	@variable(m, p_tik[t=1:size(T, 1), i=1:2], start = p0_ti[t,i])
	@variable(m, q_tik[t=1:size(T, 1), i=1:2], start = q0_ti[t,i])
	@variable(m, SSh[!,:Blb][s] <= b_sk[s=1:size(SSh, 1)] <= SSh[!,:Bub][s], start = b0_s[s])
	@variable(m, G[!,:Plb][g] <= p_gk[g=1:size(G, 1)] <= G[!,:Pub][g], start = p0_g[g])
	@variable(m, G[!,:Qlb][g] <= q_gk[g=1:size(G, 1)] <= G[!,:Qub][g], start = q0_g[g])
	@variable(m, pslackm_nk[n=1:size(N, 1)] >= 0, start = pslack0m_n[n])
	@variable(m, pslackp_nk[n=1:size(N, 1)] >= 0, start = pslack0p_n[n])
	@variable(m, qslackm_nk[n=1:size(N, 1)] >= 0, start = qslack0m_n[n])
	@variable(m, qslackp_nk[n=1:size(N, 1)] >= 0, start = qslack0p_n[n])
	@variable(m, sslack_lik[l=1:size(L, 1), i=1:2] >= 0, start = sslack0_li[l,i])
	@variable(m, sslack_tik[t=1:size(T, 1), i=1:2] >= 0, start = sslack0_ti[t,i])
	@variable(m, delta_k, start=delta0)
	
	# unavailability
	if K[1] == :Generator
		JuMP.fix(p_gk[K_outidx], 0.0, force = true)
		JuMP.fix(q_gk[K_outidx], 0.0, force = true)
		JuMP.set_start_value(p_gk[K_outidx], 0.0)
		JuMP.set_start_value(q_gk[K_outidx], 0.0)
	elseif K[1] == :Line
		for i = 1:2
			JuMP.fix(p_lik[K_outidx, i], 0.0, force = true)
			JuMP.fix(q_lik[K_outidx, i], 0.0, force = true)
			JuMP.set_start_value(p_lik[K_outidx, i], 0.0)
			JuMP.set_start_value(q_lik[K_outidx, i], 0.0)
		end
	elseif K[1] == :Transformer
		for i = 1:2
			JuMP.fix(p_tik[K_outidx, i], 0.0, force = true)
			JuMP.fix(q_tik[K_outidx, i], 0.0, force = true)
			JuMP.set_start_value(p_tik[K_outidx, i], 0.0)
			JuMP.set_start_value(q_tik[K_outidx, i], 0.0)
		end
	else
		warn("unrecognized contingency type ", K[1])
	end
	
	# fix contingency angle at reference bus to zero
	JuMP.fix(theta_nk[RefBus], 0.0, force=true)
	
	# add power flow constraints for contingency
	addpowerflowcons!(m, v_nk, theta_nk, p_lik, q_lik, p_tik, q_tik, b_sk, p_gk, q_gk,
		pslackm_nk, pslackp_nk, qslackm_nk, qslackp_nk, sslack_lik, sslack_tik,
		N, L, T, (L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin,
		SShn, Gn, K_outidx), o.ThermalLimits, SysCond=:Contingency,
		ConType=K[1], outidx=K_outidx)
	
	# register approximation functions if necessary
	if o.CouplingMode == :approx
                SmoothingParam = o.SmoothingParamCoupling
		sclamp(x) = smoothclamp(x, SmoothingParam)
		JuMP.register(m, :sclamp, 1, sclamp,
			x -> SmoothApproximations.smoothclampprime(x, SmoothingParam),
			x -> SmoothApproximations.smoothclampprimeprime(x, SmoothingParam))
		sstep(x) = smoothstep(x, SmoothingParam)
		JuMP.register(m, :sstep, 1, sstep,
			x -> SmoothApproximations.smoothstepprime(x, SmoothingParam),
			x -> SmoothApproximations.smoothstepprimeprime(x, SmoothingParam))
	else
		sclamp = nothing
		sstep = nothing
	end
	
	# add coupling constraints with base case
	if o.CouplingMode != :ignore
		addcouplingcons!(m, o, v_n, p_g, v_nk, p_gk, q_gk, delta_k,
			N, G, K[1], K_outidx, 
			(L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin,
			Tidxn, Tin, SShn, Gn, K_outidx), sclamp=sclamp,
			sstep=sstep, SetStart=true)
	end
	
	## objective function
	conpenalty =
		if o.PenaltyType == :pcwslin
			addpenaltyfunctions!(m, pslackm_nk, pslackp_nk, qslackm_nk, qslackp_nk,
				sslack_lik, sslack_tik, N, L, T, P, true)
		else
			quadpenaltyfunctions(pslackm_nk, pslackp_nk, qslackm_nk, qslackp_nk,
				sslack_lik, sslack_tik, P)
		end
	@objective(m, Min, conpenalty)
	#@objective(m, Min, sum(pslackm_nk[n] + pslackm_nk[n] for n=1:size(N, 1)) +
	#	sum(qslackm_nk[n] + qslackp_nk[n] for n=1:size(N, 1)) +
	#	sum(sslack_lik[l,i] for l=1:size(L,1), i=1:2) +
	#	sum(sslack_tik[t,i] for t=1:size(T,1), i=1:2))
	
	## attempt to solve contingency residual problem
	JuMP.optimize!(m)
        solvestatusprimal = JuMP.primal_status(m)
        solvestatus = JuMP.termination_status(m)
        println("JuMP optimize MOI termination_status: ", solvestatus, " ", solvestatusprimal)
	if (solvestatusprimal != MOI.FEASIBLE_POINT) && (solvestatusprimal != MOI.NEARLY_FEASIBLE_POINT)
		error("solver failed to find a feasible solution")
	end
	
	#@show sum(JuMP.value.(pslackm_nk) + JuMP.value.(pslackp_nk))
	#@show sum(JuMP.value.(qslackm_nk) + JuMP.value.(qslackp_nk))
	#@show sum(JuMP.value.(sslack_lik))
	#@show sum(JuMP.value.(sslack_tik))
	#@show JuMP.value(conpenalty)
	#@show JuMP.objective_value(m)
	#@show conpenalty
	#@show JuMP.objective_function(m, AffExpr)
        if false
            vp_k = JuMP.value.(p_gk)
            vp_k = view(vp_k, :, 1)
            vp   = p_g
            vPlb = G[!,:Plb]
            vPub = G[!,:Pub]
            delta= JuMP.value.(delta_k)
            delta= delta[1]
            alpha = G[!,:alpha]
            for i=1:length(vp_k)
                @printf("g=%4d pg=%12.5e  pg_k=%12.5e diffPlb=%12.5e diffPub=%12.5e  p_k-(p+alpha*delta)=%12.5e fPlb=%12.5e Pub=%12.5e alpha=%12.5e  delta=%12.5e\n", 
                        i, vp[i], vp_k[i], vp_k[i]-vPlb[i], vPub[i]-vp_k[i], vp_k[i]-(vp[i]+alpha[i]*delta), vPlb[i], vPub[i], alpha[i], delta)
                #@printf("alpha[i]=%g  delta=%g\n", alpha[i], delta)
            end
        end

        if false
            vq_gk = JuMP.value.(q_gk)
            #vq_gk = view(vq_gk,:,1)
            #vq_g  = q_g
            vv_nk = JuMP.value.(v_nk)
            #vv_nk = view(vv_nk,:,1)
            vv_n  = v_n
            vQlb = G[!,:Qlb]
            vQub = G[!,:Qub]
            
            for i=1:length(vq_gk)
                @printf("PVPQ gen=%3d q_gk=%12.5e v_n=%12.5e v_nk=%12.5e | v-v_k=%12.5e Qub-q_k=%12.5e q_k-Qlb=%12.5e | Qlb=%12.5e Qub=%12.5e\n",
                        i, vq_gk[i], vv_n[i], vv_nk[i], vv_n[i]- vv_nk[i],  vQub[i]-vq_gk[i],  vq_gk[i]-vQlb[i], vQlb[i], vQub[i])
            end
                    
        end
        @printf("Backend: ContingencyID=%d objs: conpenalty=%11.4e jumpobj=%11.4e\n", 
                K[3], JuMP.value(conpenalty), JuMP.objective_value(m)); flush(stdout)
        solution = ContingencyState(K[3], K[2], K[1], JuMP.value(delta_k), JuMP.value.(v_nk), 
                                    JuMP.value.(theta_nk), convert(Vector{Float64}, JuMP.value.(b_sk)), 
                                    JuMP.value.(p_gk), JuMP.value.(q_gk),
                                    JuMP.value(conpenalty), JuMP.objective_value(m))

        return solution
	#return JuMP.value(delta_k), JuMP.value.(v_nk), JuMP.value.(theta_nk),
	#	convert(Vector{Float64}, JuMP.value.(b_sk)), JuMP.value.(p_gk),
	#	JuMP.value.(q_gk), JuMP.value(conpenalty)
	#
end

# two-step contingency solve
#  - step 1: solve approximatively using one of the smoothing approaches for coupling constraints
#  - step 2: crash the solution from step1, fix "half" of the contingency coupling variables and resolve
function SolveContingencyTwoStep(o_in::GoOptions, N::DataFrame, L::DataFrame, T::DataFrame,
	                         SSh::DataFrame, G::DataFrame, K::Tuple{Symbol, Int, Int}, P::DataFrame,
                                 scacopfsol::SCACOPFState,
	                         NLSolver; 
                                 IndexSets::Union{Tuple, Nothing}=nothing,
                                 StartingPoint::Union{Nothing,ContingencyState}=nothing)
    o = o_in

    println(" -- SolveContingencyTwoStep: ",
            "ContingencyID/IDOut/Type=", K[3], "/", K[2], "/", K[1], "\n",
            "Options: SolveType=", o.PostSolveType, " CouplingMode=", o.CouplingMode,
            " ThermalLimits=", o.ThermalLimits, " SmoothingParam=", o.SmoothingParamCoupling);
   
    
    if haskey(scacopfsol.cont_states, K[3])
        #this contingency was part of the scacopf
        StartingPoint = scacopfsol.cont_states[K[3]]
        println("using conting solution from SC-ACOPF")
    end

    if StartingPoint == nothing
        # this is the first contingency to be solved
        # use options sent by the caller
    else
        o.ipopt.mu_init = 1e-4
        o.ipopt.max_cpu_time = 3600.
        o.ipopt.max_iter = 1000
    end
    NLSolver = with_ipopt_optimizer(o)

    ConSolApprox = StartingPoint
    if ConSolApprox != nothing
        ConSolApprox.cont_id = K[3]
    end

    tries=0
    while tries<=1
        tries +=1
        
        try
            ConSolApprox = SolveContingencyBackend(o, N, L, T, SSh, G, K, P, scacopfsol, NLSolver,
			                      StartingPoint=StartingPoint, IndexSets=IndexSets)
            break
        catch ex
            println(ex)
            println("SolveContingencyBackend[1] exception encountered, relaxing coupling")
            o.SmoothingParamCoupling *= 100
            o.ipopt.mu_init          *= 10
            o.ipopt.tol               = 1e-7
            o.ipopt.max_cpu_time = 3600.
            o.ipopt.max_iter = 600
            NLSolver = with_ipopt_optimizer(o)
        end
    end

    #
    ## Pass with fixing - TRY1
    #
    fixingOK=true #be optimistic
    o.ipopt.mu_init = 1e-1
    o.ipopt.tol     = 1e-8
    ConSol3 = ConSolApprox

    o.ipopt.max_cpu_time = 3600.
    o.ipopt.max_iter = 1000

    tries = 0; maxtries=1
    while tries<maxtries
        tries += 1

        NLSolver = with_ipopt_optimizer(o)
        ConSol3 = ConSolApprox
        try
            ConSol3 = SolveContingencyWithFixing(o, N, L, T, SSh, G, K, P, scacopfsol, NLSolver;
	                                         StartingPoint=ConSolApprox, IndexSets=IndexSets)

            fixingOK=true
        catch ex
            println(ex)
            println("SolveContingencyWithFixing[3] TRY1 exception encountered, relaxing mu (fixVolt=false)")
            fixingOK=false
            o.ipopt.mu_init *=50
            o.ipopt.max_cpu_time = 3600.
            o.ipopt.max_iter = 500
        end

        if fixingOK
            if ConSol3.objective>1e+7 && ConSolApprox.objective<=ConSol3.objective/1000
                #solve smoothing problem with increased smoothing 
                o.SmoothingParamCoupling *= 20
                o.ipopt.max_iter = 1000
                NLSolver = with_ipopt_optimizer(o)
                @printf("TwoStep: ContingencyID=%d detected large fixing obj=%g | backend obj=%g. will relax backend\n",
                        K[3], ConSol3.objective, ConSolApprox.objective)
                try
                    ConSolApprox = SolveContingencyBackend(o, N, L, T, SSh, G, K, P, scacopfsol, NLSolver,
			                                   StartingPoint=StartingPoint, IndexSets=IndexSets)
                    ConSol3 = SolveContingencyWithFixing(o, N, L, T, SSh, G, K, P, scacopfsol, NLSolver;
	                                                 StartingPoint=ConSolApprox, IndexSets=IndexSets)
                    @printf("TwoStep: ContingencyID=%d repaired large new fixing obj=%g | backend obj=%g\n",
                        K[3], ConSol3.objective, ConSolApprox.objective)
                    break
                catch ex
                    println(ex)
                    println("SolveContingencyWithFixing[3] TRY2 exception with relaxed backend and fixing")
                end
            else
                break
            end
        end
    end # of while

    if fixingOK
        return (ConSol3, ConSolApprox)
    end

    # #
    # ## Pass with fixing - TRY2 - this time fix voltages
    # #
    # fixingOK=true #be optimistic
    # o.ipopt.mu_init = 1e-1
    # o.ipopt.tol     = 1e-8
    # o.ipopt.max_cpu_time = 1200.
    # o.ipopt.max_iter = 500
    # ConSol3 = ConSolApprox

    # tries = 0
    # while tries<0
    #     tries += 1
    #     NLSolver = with_ipopt_optimizer(o)
    #     ConSol3 = ConSolApprox
    #     try
    #         ConSol3 = SolveContingencyWithFixing(o, N, L, T, SSh, G, K, P, scacopfsol, NLSolver;
    #                                                  StartingPoint=ConSolApprox, IndexSets=IndexSets, fixVolt=true)
    #         fixingOK=true
    #         break
    #     catch ex
    #         println(ex)
    #         println("SolveContingencyWithFixing TRY2 exception encountered, increasing mu (fixVolt=true)")
    #         fixingOK=false
    #         o.ipopt.mu_init *=50
    #     end
    # end

    # if fixingOK
    #     return (ConSol3, ConSolApprox)
    # end

    #
    ## Pass with fixing - TRY3 - this time relax smoothing parameter and resolve the smooth problem
    #
    #relax smoothing from the get-go
    o.SmoothingParamCoupling *= 50
    fixingOK=true #be optimistic
    mu_init_approx = 1e-4
    o.ipopt.tol     = 1e-8

    o.ipopt.max_iter = 1200
    o.ipopt.max_cpu_time = 3600.
    ConSol3 = ConSolApprox

    tries = 0
    while tries<=2
        tries += 1
        ConSol3 = ConSolApprox

        try
            o.ipopt.mu_init = mu_init_approx            
            NLSolver = with_ipopt_optimizer(o)
            ConSolApprox = SolveContingencyBackend(o, N, L, T, SSh, G, K, P, scacopfsol, NLSolver,
			                      StartingPoint=StartingPoint, IndexSets=IndexSets)

            o.ipopt.mu_init = 1.0
            NLSolver = with_ipopt_optimizer(o)
            ConSol3 = SolveContingencyWithFixing(o, N, L, T, SSh, G, K, P, scacopfsol, NLSolver;
	                                             StartingPoint=ConSolApprox, IndexSets=IndexSets)
            fixingOK = true
            break
        catch ex
            println(ex)
            println("SolveContingencyBackend+WithFixing TRY3 exception encountered, relaxing coupling")
            fixingOK = false
            o.SmoothingParamCoupling *= 10
            mu_init_approx           *= 10
            o.ipopt.max_iter = 800
            o.ipopt.max_cpu_time = 3600.
        end 
    end

    if false==fixingOK
        println("SolveContingencyTwoStep: approximate solution returned for ContingencyID/IDOut/Type=", K[3], "/", K[2], "/", K[1])
    end
    

    #return (sol with exact coupling, sol with smoothing); the last one will be used to warm start the next contingency
    return (ConSol3, ConSolApprox)
end


## function for solving a contingency given:
##	- base case voltages
##	- base case dispatch

function SolveContingencyTwoStep_old(o::GoOptions, N::DataFrame, L::DataFrame, T::DataFrame,
	                         SSh::DataFrame, G::DataFrame, K::Tuple{Symbol, Int, Int}, P::DataFrame,
                                 scacopfsol::SCACOPFState,
	                         #v_n::AbstractVector{Float64}, p_g::AbstractVector{Float64},
	                         NLSolver; IndexSets::Union{Tuple, Nothing}=nothing)

        println(" -- SolveContingencyTwoStep: ",
                "ContingencyID/IDOut/Type=", K[3], "/", K[2], "/", K[1], "\n    ",
                "Options: SolveType=", o.PostSolveType, " CouplingMode=", o.CouplingMode,
                " ThermalLimits=", o.ThermalLimits, " SmoothingParam=", o.SmoothingParamCoupling);

        v_n = scacopfsol.base_state.v_n
        p_g = scacopfsol.base_state.p_g

	# greedy fixing: delta=P0gcon/sum(P0g), vk=v0
	if K[1] == :Generator
		if IndexSets == nothing
			gout = Int(findfirst(G[!,:Generator] .== K[2]))
		else
			gout = IndexSets[11]
		end
		Gkind = [1:(gout-1); (gout+1):size(G, 1)]
		delta_k = find_zero(x -> redispatch(x, view(p_g, Gkind), view(G[!,:alpha], Gkind),
			view(G[!,:Plb], Gkind), view(G[!,:Pub], Gkind)) - p_g[gout], 0.0, Order1())
	else delta_k = 0.0
	end
	vapprox_nk = v_n
	
	if o.PostSolveType == :fast	
		
		# call contingency function
		#return delta_k, SolveContingency(o, N, L, T, SSh, G, K, P, v_n, p_g, vapprox_nk, delta_k,
		#	NLSolver, IndexSets=IndexSets)...
                ContSol = SolveContingencyFastBackend(o, N, L, T, SSh, G, K, P, 
                                                      scacopfsol, NLSolver, 
                                                      IndexSets=IndexSets, DeltaToUse=delta_k)
                return ContSol
	elseif o.PostSolveType == :optimized
		
		# solve contingency using fast approach to obtain a starting point
		#ConSolFast = SolveContingency(o, N, L, T, SSh, G, K, P, v_n, p_g, vapprox_nk, delta_k,
		#	NLSolver, IndexSets=IndexSets)
		ConSolFast = SolveContingencyFastBackend(o, N, L, T, SSh, G, K, P, 
                                                         scacopfsol, NLSolver, 
                                                         IndexSets=IndexSets, DeltaToUse=delta_k)
		
                o.SmoothingParamCoupling *= 1000;
		# compute contingency solution using optimization
		#ConSol = SolveContingencyBackend(o, N, L, T, SSh, G, K, P, v_n, p_g, NLSolver,
		#	StartingPoint=(delta_k, ConSolFast...), IndexSets=IndexSets)
		ConSol = SolveContingencyBackend(o, N, L, T, SSh, G, K, P, scacopfsol, NLSolver,
			StartingPoint=ConSolFast, IndexSets=IndexSets)
                o.SmoothingParamCoupling /= 1000;
		
		# call exact contingency function to enforce coupling exactly if necessary
		
		if o.CouplingMode != :approx
			return ConSol
		else
			#return ConSol[1], SolveContingency(o, N, L, T, SSh, G, K, P, v_n, p_g,
			#	ConSol[2], ConSol[1], NLSolver, ConSol[3], ConSol[4], ConSol[5],
			#	IndexSets=IndexSets)...
		    return SolveContingencyFastBackend(o, N, L, T, SSh, G, K, P, 
                                                       scacopfsol, NLSolver, 
				                       IndexSets=IndexSets, DeltaToUse=ConSol.delta)
		end
		
	else
		warn("unrecognized SolveType value ", o.PostSolveType)
	end
	
end

## function to solve all contingencies given:
##	- base case voltages
##	- base case dispatch

# function SolveAllContingencies(o::GoOptions, N::DataFrame, L::DataFrame, T::DataFrame,
# 	SSh::DataFrame, G::DataFrame, K::DataFrame, P::DataFrame,
# 	v_n::AbstractVector{Float64}, p_g::AbstractVector{Float64}, NLSolver;
# 	IndexSets::Union{Tuple, Nothing}=nothing)

#         println("SolveAllContingencies(base case only): Options: SolveType=", o.PostSolveType,
#                 " CouplingMode=", o.CouplingMode,
#                 " ThermalLimits=", o.ThermalLimits,
#                 " SmoothingParam=", o.SmoothingParamCoupling);
	
# 	# compute index and sets for performant formulation
# 	if IndexSets == nothing
# 		L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn, Lin, Tidxn, Tin, SShn, Gn,
# 			K_outidx = indexsets(N, L, T, SSh, G, K)
# 	else
# 		L_Nidx = IndexSets[1]
# 		T_Nidx = IndexSets[2]
# 		SSh_Nidx = IndexSets[3]
# 		G_Nidx = IndexSets[4]
# 		Lidxn = IndexSets[5]
# 		Lin = IndexSets[6]
# 		Tidxn = IndexSets[7]
# 		Tin = IndexSets[8]
# 		SShn = IndexSets[9]
# 		Gn = IndexSets[10]
# 		K_outidx = IndexSets[11]
# 	end
	
# 	# allocate space to hold the contingency solution
# 	delta_k = zeros(Float64, size(K, 1))
# 	v_nk = Array{Float64, 2}(undef, size(N, 1), size(K, 1))
# 	theta_nk = Array{Float64, 2}(undef, size(N, 1), size(K, 1))
# 	b_sk = Array{Float64, 2}(undef, size(SSh, 1), size(K, 1))
# 	p_gk = Array{Float64, 2}(undef, size(G, 1), size(K, 1))
# 	q_gk =  Array{Float64, 2}(undef, size(G, 1), size(K, 1))
	
# 	# solve contingencies and return
# 	for k = 1:size(K, 1)
# 		print("Solving contingency ", k, " out of ", size(K, 1), "... ")
# 		contime = @elapsed begin
# 			delta_k[k], v_nk[:,k], theta_nk[:,k], b_sk[:,k], p_gk[:,k], q_gk[:,k] =
# 				SolveContingency(o, N, L, T, SSh, G, (K[!,:ConType][k], K[:IDout][k]), P,
# 				                 v_n, p_g, NLSolver,
# 				                 IndexSets=(L_Nidx, T_Nidx, SSh_Nidx, G_Nidx, Lidxn,
# 				                            Lin, Tidxn, Tin, SShn, Gn, K_outidx[k]))
# 		end
# 		print("Solved in ", round(contime, digits=2), " secs.\n")
# 		flush(stdout)
# 	end
# 	return v_nk, theta_nk, b_sk, p_gk, q_gk, delta_k
# end

end
