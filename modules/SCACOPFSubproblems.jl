#__precompile__()

module SCACOPFSubproblems

## elements to be exported

export solve_base_power_flow, solve_basecase, solve_contingency, solve_random_contingency,
       solve_SC_ACOPF, 
       SCACOPFdata, GenericContingency, 
       SubproblemSolution, BasecaseSolution, ContingencySolution, SCACOPFsolution,
       write_solution

## load external modules

using JuMP, Printf, DataFrames, Random, Distributions
const MOI = JuMP.MOI
import JuMP: add_to_expression!

## load internal modules and functions

Base.include(@__MODULE__, "corejugo/CoreJuGO.jl")
using .CoreJuGO

## auxiliary functions and other definitions

include("SCACOPFSubproblems/starting_point.jl")
include("SCACOPFSubproblems/power_flow_constraints.jl")
include("SCACOPFSubproblems/contingency_n-k.jl")
include("SCACOPFSubproblems/huber_loss.jl")

## constants

## module body

# function to solve power flow (assumes x0 gives setpoint)

function solve_base_power_flow(psd::SCACOPFdata, NLSolver)
    
    # get primal starting point
    x0 = get_primal_starting_point(psd)
    
    # create model
    m = Model(NLSolver)

    # base case variables
    @variable(m, psd.N[n,:Vlb] <= v_n[n=1:nrow(psd.N)] <= psd.N[n,:Vub],
              start=x0[:v_n][n])
    @variable(m, theta_n[n=1:nrow(psd.N)], start=x0[:theta_n][n])
    @variable(m, p_li[l=1:nrow(psd.L), i=1:2], start=x0[:p_li][l,i])
    @variable(m, q_li[l=1:nrow(psd.L), i=1:2], start=x0[:q_li][l,i])
    @variable(m, p_ti[t=1:nrow(psd.T), i=1:2], start=x0[:p_ti][t,i])
    @variable(m, q_ti[t=1:nrow(psd.T), i=1:2], start=x0[:q_ti][t,i])
    @variable(m, psd.SSh[s,:Blb] <= b_s[s=1:nrow(psd.SSh)] <=
              psd.SSh[s,:Bub], start=x0[:b_s][s])
    @variable(m, psd.G[g,:Plb] <= p_g[g=1:nrow(psd.G)] <= psd.G[g,:Pub],
              start=x0[:p_g][g])
    @variable(m, psd.G[g,:Qlb] <= q_g[g=1:nrow(psd.G)] <= psd.G[g,:Qub],
              start=x0[:q_g][g])
    @variable(m, p_devm_g[g=1:nrow(psd.G)] >= 0, start=0.0)         # negative deviation w.r.t. setpoint
    @variable(m, p_devp_g[g=1:nrow(psd.G)] >= 0, start=0.0)         # positive deviation w.r.t. setpoint
    @variable(m, pslackm_n[n=1:size(psd.N, 1)] >= 0, start = x0[:pslackm_n][n])
    @variable(m, pslackp_n[n=1:size(psd.N, 1)] >= 0, start = x0[:pslackp_n][n])
    @variable(m, qslackm_n[n=1:size(psd.N, 1)] >= 0, start = x0[:qslackm_n][n])
    @variable(m, qslackp_n[n=1:size(psd.N, 1)] >= 0, start = x0[:qslackp_n][n])
    @variable(m, sslack_li[l=1:nrow(psd.L), i=1:2] >= 0,
              start=x0[:sslack_li][l,i])
    @variable(m, sslack_ti[t=1:nrow(psd.T), i=1:2] >= 0,
              start=x0[:sslack_ti][t,i])
    
    # fix angle at reference bus to zero
    JuMP.fix(theta_n[psd.RefBus], 0.0, force=true)
    
    # add power flow constraints
    addpowerflowcons!(m, v_n, theta_n, p_li, q_li, p_ti, q_ti, b_s, p_g, q_g,
                      pslackm_n, pslackp_n, qslackm_n, qslackp_n,
                      sslack_li, sslack_ti, psd)
    
    # register Huber-like deviation penalty function
    scale_factor = 1.0E-3
    dev = HuberLikePenalty(psd.a[:S] * scale_factor, psd.b[:S] * scale_factor,
                           2 * psd.a[:S] * scale_factor * 0.5 + psd.b[:S] * scale_factor,   
                                                                  # slope at 0.5 pu (e.g., 50MW) is max slope
                           0.05)                                  # 0.05pu (e.g., 5MW) to change curvature
    dev_prime = HuberLikePenaltyPrime(dev)
    dev_prime_prime = HuberLikePenaltyPrimePrime(dev)
    register(m, :dev, 1, x -> dev(x), x -> dev_prime(x), x -> dev_prime_prime(x))
    
    # production deviation w.r.t. set points, except for generators at swing buses
    @constraint(m, [g=1:nrow(psd.G)], p_g[g] == x0[:p_g][g] - p_devm_g[g] + p_devp_g[g])
    p_g_dev_penalty = @NLexpression(m, sum(dev(p_devm_g[g]) for g=1:nrow(psd.G) if 
                                           psd.N[psd.G_Nidx[g],:Type] != :SWING) +
                                       sum(dev(p_devp_g[g]) for g=1:nrow(psd.G) if
                                           psd.N[psd.G_Nidx[g],:Type] != :SWING) )
    
    # register Huber-like violation penalty function
    h = HuberLikePenalty(psd.a[:S], psd.b[:S],
                         2 * psd.a[:S] * 0.1 + psd.b[:S],       # slope at 0.1 pu (e.g., 10MW) is max slope
                         0.05)                                  # 0.05pu (e.g., 5MW) to change curvature
    h_prime = HuberLikePenaltyPrime(h)
    h_prime_prime = HuberLikePenaltyPrimePrime(h)
    register(m, :h, 1, x -> h(x), x -> h_prime(x), x -> h_prime_prime(x))
    
    # collect technical violation penalty terms
    p_bal_penalty = @NLexpression(m, sum(h(pslackm_n[n]) for n=1:nrow(psd.N)) +
                                     sum(h(pslackp_n[n]) for n=1:nrow(psd.N)))
    q_bal_penalty = @NLexpression(m, sum(h(qslackm_n[n]) for n=1:nrow(psd.N)) +
                                     sum(h(qslackp_n[n]) for n=1:nrow(psd.N)))
    lin_overload_penalty = @NLexpression(m, sum(h(sslack_li[l,i]) for l=1:nrow(psd.L), i=1:2))
    trf_overload_penalty = @NLexpression(m, sum(h(sslack_ti[t,i]) for t=1:nrow(psd.T), i=1:2))
    violation_penalty = @NLexpression(m, p_bal_penalty + q_bal_penalty +
                                         lin_overload_penalty + trf_overload_penalty)
    
    # declare objective
    @NLobjective(m, Min, p_g_dev_penalty + violation_penalty)
    
    # attempt to solve SCACOPF
    JuMP.optimize!(m)
    if JuMP.primal_status(m) != MOI.FEASIBLE_POINT && 
       JuMP.primal_status(m) != MOI.NEARLY_FEASIBLE_POINT &&
       JuMP.termination_status(m) != MOI.NEARLY_FEASIBLE_POINT &&
       JuMP.termination_status(m) != MOI.ALMOST_LOCALLY_SOLVED
        @show JuMP.primal_status(m)
        @show JuMP.termination_status(m)
        error("solver failed to find a feasible solution.")
    end
    
    # aggregate deviations to report them
    summary = Dict{Symbol, Float64}()
    summary[:p_deviations] = sum(JuMP.value(p_devm_g[g]) + JuMP.value(p_devp_g[g]) for g=1:nrow(psd.G) if
                                 psd.N[psd.G_Nidx[g],:Type] != :SWING)
    
    # aggregate infeasibilities to report them
    summary[:active_nodal_imbalance] = sum(JuMP.value.(pslackm_n)) + sum(JuMP.value.(pslackp_n))
    summary[:reactive_nodal_imbalance] = sum(JuMP.value.(qslackm_n)) + sum(JuMP.value.(qslackp_n))
    summary[:branch_overloads] = .5 * sum(JuMP.value.(sslack_li)) + .5 * sum(JuMP.value.(sslack_li))
    
    # return solution
    return BasecaseSolution(psd, JuMP.value.(v_n), JuMP.value.(theta_n),
                            convert(Vector{Float64}, JuMP.value.(b_s)),
                            JuMP.value.(p_g), JuMP.value.(q_g),
                            0.0, 0.0),
           summary
    
end

# function to solve base case, possibly with recourse approximations

function solve_basecase(psd::SCACOPFdata, NLSolver;
                       recourse_f::T=nothing,   # recourse function value
                       recourse_g::T=nothing,   # recourse function gradient
                       recourse_H::T=nothing,   # recourse function hessian
                       previous_solution::Union{Nothing,
                                                BasecaseSolution}=nothing
                       )::BasecaseSolution where {T <: Union{Nothing, Function}}
    
    # get primal starting point
    x0 = get_primal_starting_point(psd, previous_solution)
    
    # create model
    m = Model(NLSolver)

    # base case variables
    @variable(m, psd.N[n,:Vlb] <= v_n[n=1:nrow(psd.N)] <= psd.N[n,:Vub],
              start=x0[:v_n][n])
    @variable(m, theta_n[n=1:nrow(psd.N)], start=x0[:theta_n][n])
    @variable(m, p_li[l=1:nrow(psd.L), i=1:2], start=x0[:p_li][l,i])
    @variable(m, q_li[l=1:nrow(psd.L), i=1:2], start=x0[:q_li][l,i])
    @variable(m, p_ti[t=1:nrow(psd.T), i=1:2], start=x0[:p_ti][t,i])
    @variable(m, q_ti[t=1:nrow(psd.T), i=1:2], start=x0[:q_ti][t,i])
    @variable(m, psd.SSh[s,:Blb] <= b_s[s=1:nrow(psd.SSh)] <=
              psd.SSh[s,:Bub], start=x0[:b_s][s])
    @variable(m, psd.G[g,:Plb] <= p_g[g=1:nrow(psd.G)] <= psd.G[g,:Pub],
              start=x0[:p_g][g])
    @variable(m, psd.G[g,:Qlb] <= q_g[g=1:nrow(psd.G)] <= psd.G[g,:Qub],
              start=x0[:q_g][g])
    @variable(m, c_g[g=1:nrow(psd.G)], start=x0[:c_g][g])
    @variable(m, pslackm_n[n=1:size(psd.N, 1)] >= 0, start = x0[:pslackm_n][n])
    @variable(m, pslackp_n[n=1:size(psd.N, 1)] >= 0, start = x0[:pslackp_n][n])
    @variable(m, qslackm_n[n=1:size(psd.N, 1)] >= 0, start = x0[:qslackm_n][n])
    @variable(m, qslackp_n[n=1:size(psd.N, 1)] >= 0, start = x0[:qslackp_n][n])
    @variable(m, sslack_li[l=1:nrow(psd.L), i=1:2] >= 0,
              start=x0[:sslack_li][l,i])
    @variable(m, sslack_ti[t=1:nrow(psd.T), i=1:2] >= 0,
              start=x0[:sslack_ti][t,i])
    
    # fix angle at reference bus to zero
    JuMP.fix(theta_n[psd.RefBus], 0.0, force=true)
    
    # add power flow constraints
    addpowerflowcons!(m, v_n, theta_n, p_li, q_li, p_ti, q_ti, b_s, p_g, q_g,
                      pslackm_n, pslackp_n, qslackm_n, qslackp_n,
                      sslack_li, sslack_ti, psd)
    
    if psd.G.CTYP[1] == 1
        c_g = Vector{JuMP.QuadExpr}(undef, size(psd.G, 1))
        # production cost (continous quadratic formulation)
        for g = 1:size(psd.G, 1)
            c_g[g] = QuadExpr( AffExpr(psd.G.COST[g], p_g[g] => psd.G.COSTLIN[g]   * psd.MVAbase),
                                UnorderedPair(p_g[g], p_g[g]) => psd.G.COSTQUAD[g] * psd.MVAbase * psd.MVAbase, )
        end
    else
        # production cost (epigraph formulation)
        for g = 1:size(psd.G, 1)
            slope_gi = psd.G_epicost_slope[g]
            intercept_gi = psd.G_epicost_intercept[g]
            @constraint(m, [i=1:length(slope_gi)],
                        c_g[g] >= slope_gi[i]*p_g[g] + intercept_gi[i])
        end
    end
    production_cost = @expression(m, sum(c_g[g] for g=1:nrow(psd.G)))
    
    # base case penalty
    basecase_penalty = JuMP.GenericQuadExpr(JuMP.AffExpr(0))
    for n = 1:nrow(psd.N)
        add_to_expression!(basecase_penalty,
                           psd.a[:P], pslackm_n[n], pslackm_n[n])
        add_to_expression!(basecase_penalty, psd.b[:P], pslackm_n[n])
        add_to_expression!(basecase_penalty,
                           psd.a[:P], pslackp_n[n], pslackp_n[n])
        add_to_expression!(basecase_penalty, psd.b[:P], pslackp_n[n])
        add_to_expression!(basecase_penalty,
                           psd.a[:Q], qslackm_n[n], qslackm_n[n])
        add_to_expression!(basecase_penalty, psd.b[:Q], qslackm_n[n])
        add_to_expression!(basecase_penalty,
                           psd.a[:Q], qslackp_n[n], qslackp_n[n])
        add_to_expression!(basecase_penalty, psd.b[:Q], qslackp_n[n])
    end
    for l = 1:nrow(psd.L), i=1:2
        add_to_expression!(basecase_penalty,
                           psd.a[:S], sslack_li[l,i], sslack_li[l,i])
        add_to_expression!(basecase_penalty, psd.b[:S], sslack_li[l,i])
    end
    for t = 1:nrow(psd.T), i=1:2
        add_to_expression!(basecase_penalty,
                           psd.a[:S], sslack_ti[t,i], sslack_ti[t,i])
        add_to_expression!(basecase_penalty, psd.b[:S], sslack_ti[t,i])
    end
    
    # contingency penalty
    if isnothing(recourse_f)
        contingency_penalty = 0.0
    else
        JuMP.register(m, :recourse_f, nrow(psd.G),
                      recourse_f, recourse_g, recourse_H)
        contingency_penalty = @NLexpression(m, recourse_f(p_g...))
    end
    
    # declare objective
    @objective(m, Min, production_cost + psd.delta*basecase_penalty +
               (1-psd.delta)*contingency_penalty)
    
    # attempt to solve SCACOPF
    JuMP.optimize!(m)
    if JuMP.primal_status(m) != MOI.FEASIBLE_POINT &&
       JuMP.primal_status(m) != MOI.NEARLY_FEASIBLE_POINT && 
       JuMP.termination_status(m) != MOI.NEARLY_FEASIBLE_POINT
        error("solver failed to find a feasible solution.")
    end
    
    # objective breakdown
    base_cost = JuMP.value(production_cost) +
                psd.delta*JuMP.value(basecase_penalty)
    recourse_cost = JuMP.objective_value(m) - base_cost
    
    # return solution
    return BasecaseSolution(psd, JuMP.value.(v_n), JuMP.value.(theta_n),
                            convert(Vector{Float64}, JuMP.value.(b_s)),
                            JuMP.value.(p_g), JuMP.value.(q_g),
                            base_cost, recourse_cost)
    
end

function solve_SC_ACOPF(psd::SCACOPFdata, NLSolver;
                       previous_solution::Union{Nothing,
                                                BasecaseSolution}=nothing,
                       quadratic_relaxation_k::Float64=Inf,
                       minutes_since_base::Float64=1.0,
                       use_huber_like_penalty::Bool=true
                       )::SCACOPFsolution 

    # get primal starting point
    x0 = get_primal_starting_point(psd, previous_solution)

    # create model
    m = Model(NLSolver)

    # base case variables
    @variable(m, psd.N[n,:Vlb] <= v_n[n=1:nrow(psd.N)] <= psd.N[n,:Vub],
              start=x0[:v_n][n])
    @variable(m, theta_n[n=1:nrow(psd.N)], start=x0[:theta_n][n])
    @variable(m, p_li[l=1:nrow(psd.L), i=1:2], start=x0[:p_li][l,i])
    @variable(m, q_li[l=1:nrow(psd.L), i=1:2], start=x0[:q_li][l,i])
    @variable(m, p_ti[t=1:nrow(psd.T), i=1:2], start=x0[:p_ti][t,i])
    @variable(m, q_ti[t=1:nrow(psd.T), i=1:2], start=x0[:q_ti][t,i])
    @variable(m, psd.SSh[s,:Blb] <= b_s[s=1:nrow(psd.SSh)] <=
              psd.SSh[s,:Bub], start=x0[:b_s][s])
    @variable(m, psd.G[g,:Plb] <= p_g[g=1:nrow(psd.G)] <= psd.G[g,:Pub],
              start=x0[:p_g][g])
    @variable(m, psd.G[g,:Qlb] <= q_g[g=1:nrow(psd.G)] <= psd.G[g,:Qub],
              start=x0[:q_g][g])
    @variable(m, c_g[g=1:nrow(psd.G)], start=x0[:c_g][g])
    @variable(m, pslackm_n[n=1:size(psd.N, 1)] >= 0, start = x0[:pslackm_n][n])
    @variable(m, pslackp_n[n=1:size(psd.N, 1)] >= 0, start = x0[:pslackp_n][n])
    @variable(m, qslackm_n[n=1:size(psd.N, 1)] >= 0, start = x0[:qslackm_n][n])
    @variable(m, qslackp_n[n=1:size(psd.N, 1)] >= 0, start = x0[:qslackp_n][n])
    @variable(m, sslack_li[l=1:nrow(psd.L), i=1:2] >= 0,
              start=x0[:sslack_li][l,i])
    @variable(m, sslack_ti[t=1:nrow(psd.T), i=1:2] >= 0,
              start=x0[:sslack_ti][t,i])

    # fix angle at reference bus to zero
    JuMP.fix(theta_n[psd.RefBus], 0.0, force=true)

    # add power flow constraints
    addpowerflowcons!(m, v_n, theta_n, p_li, q_li, p_ti, q_ti, b_s, p_g, q_g,
                     pslackm_n, pslackp_n, qslackm_n, qslackp_n,
                     sslack_li, sslack_ti, psd)

    if psd.G.CTYP[1] == 1
        c_g = Vector{JuMP.QuadExpr}(undef, size(psd.G, 1))
        # production cost (continous quadratic formulation)
        for g = 1:size(psd.G, 1)
            c_g[g] = QuadExpr( AffExpr(psd.G.COST[g], p_g[g] => psd.G.COSTLIN[g]   * psd.MVAbase),
                                UnorderedPair(p_g[g], p_g[g]) => psd.G.COSTQUAD[g] * psd.MVAbase * psd.MVAbase, )
        end
    else
        # production cost (epigraph formulation)
        for g = 1:size(psd.G, 1)
            slope_gi = psd.G_epicost_slope[g]
            intercept_gi = psd.G_epicost_intercept[g]
            @constraint(m, [i=1:length(slope_gi)],
                        c_g[g] >= slope_gi[i]*p_g[g] + intercept_gi[i])
        end
    end
    production_cost = @expression(m, sum(c_g[g] for g=1:nrow(psd.G)))

    # base case penalty
    basecase_penalty = JuMP.GenericQuadExpr(JuMP.AffExpr(0))
    for n = 1:nrow(psd.N)
        add_to_expression!(basecase_penalty,
                           psd.a[:P], pslackm_n[n], pslackm_n[n])
        add_to_expression!(basecase_penalty, psd.b[:P], pslackm_n[n])
        add_to_expression!(basecase_penalty,
                           psd.a[:P], pslackp_n[n], pslackp_n[n])
        add_to_expression!(basecase_penalty, psd.b[:P], pslackp_n[n])
        add_to_expression!(basecase_penalty,
                           psd.a[:Q], qslackm_n[n], qslackm_n[n])
        add_to_expression!(basecase_penalty, psd.b[:Q], qslackm_n[n])
        add_to_expression!(basecase_penalty,
                           psd.a[:Q], qslackp_n[n], qslackp_n[n])
        add_to_expression!(basecase_penalty, psd.b[:Q], qslackp_n[n])
    end
    for l = 1:nrow(psd.L), i=1:2
        add_to_expression!(basecase_penalty,
                          psd.a[:S], sslack_li[l,i], sslack_li[l,i])
        add_to_expression!(basecase_penalty, psd.b[:S], sslack_li[l,i])
    end
    for t = 1:nrow(psd.T), i=1:2
        add_to_expression!(basecase_penalty,
                          psd.a[:S], sslack_ti[t,i], sslack_ti[t,i])
        add_to_expression!(basecase_penalty, psd.b[:S], sslack_ti[t,i])
    end

    # Begin contingency case
    if !use_huber_like_penalty
        contingency_penalty = Vector{JuMP.GenericQuadExpr}(undef, nrow(psd.K))
    else
        contingency_penalty = Any[]

        # register Huber-like penalty functions
        hP = HuberLikePenalty(psd.a[:P], psd.b[:P],
                            2 * psd.a[:P] * (10.0/psd.MVAbase) + psd.b[:P], # slope at 10MW is max
                            5.0/psd.MVAbase)                                # 5MW to change curvature
        hP_prime = HuberLikePenaltyPrime(hP)
        hP_prime_prime = HuberLikePenaltyPrimePrime(hP)
        register(m, :hP, 1, x -> hP(x), x -> hP_prime(x), x -> hP_prime_prime(x))
        hQ = HuberLikePenalty(psd.a[:Q], psd.b[:Q],
                            2 * psd.a[:Q] * (10.0/psd.MVAbase) + psd.b[:Q], # slope at 10MVAr is max
                            5.0/psd.MVAbase)                                # 5MVAr to change curvature
        hQ_prime = HuberLikePenaltyPrime(hQ)
        hQ_prime_prime = HuberLikePenaltyPrimePrime(hQ)
        register(m, :hQ, 1, x -> hQ(x), x -> hQ_prime(x), x -> hQ_prime_prime(x))
        hS = HuberLikePenalty(psd.a[:S], psd.b[:S],
                            2 * psd.a[:S] * (10.0/psd.MVAbase) + psd.b[:S], # slope at 10MVA is max
                            5.0/psd.MVAbase)                                # 5MVA to change curvature
        hS_prime = HuberLikePenaltyPrime(hS)
        hS_prime_prime = HuberLikePenaltyPrimePrime(hS)
        register(m, :hS, 1, x -> hS(x), x -> hS_prime(x), x -> hS_prime_prime(x))
    end

    # quadratic relaxation term
    if quadratic_relaxation_k < Inf
        quadratic_relaxation_term = Vector{JuMP.GenericQuadExpr}(undef, nrow(psd.K))
    else
        quadratic_relaxation_term = Any[]
    end

    # contingency variables
    @variable(m, v_nk[n=1:nrow(psd.N), k=1:nrow(psd.K)] )
    @variable(m, theta_nk[n=1:nrow(psd.N), k=1:nrow(psd.K)])
    @variable(m, p_lik[l=1:nrow(psd.L), i=1:2, k=1:nrow(psd.K)])
    @variable(m, q_lik[l=1:nrow(psd.L), i=1:2, k=1:nrow(psd.K)])
    @variable(m, p_tik[t=1:nrow(psd.T), i=1:2, k=1:nrow(psd.K)])
    @variable(m, q_tik[t=1:nrow(psd.T), i=1:2, k=1:nrow(psd.K)])
    @variable(m, b_sk[s=1:nrow(psd.SSh), k=1:nrow(psd.K)])
    @variable(m, p_gk[g=1:nrow(psd.G), k=1:nrow(psd.K)])
    @variable(m, q_gk[g=1:nrow(psd.G), k=1:nrow(psd.K)])
    @variable(m, pslackm_nk[n=1:size(psd.N, 1), k=1:nrow(psd.K)])
    @variable(m, pslackp_nk[n=1:size(psd.N, 1), k=1:nrow(psd.K)])
    @variable(m, qslackm_nk[n=1:size(psd.N, 1), k=1:nrow(psd.K)])
    @variable(m, qslackp_nk[n=1:size(psd.N, 1), k=1:nrow(psd.K)])
    @variable(m, sslack_lik[l=1:nrow(psd.L), i=1:2, k=1:nrow(psd.K)])
    @variable(m, sslack_tik[t=1:nrow(psd.T), i=1:2, k=1:nrow(psd.K)])

    # contingency variables constraint
    @constraint(m, [n=1:nrow(psd.N), k = 1:nrow(psd.K)], psd.N[n,:EVlb] <= v_nk[n, k] <= 
            psd.N[n,:EVub])
    @constraint(m, [s=1:nrow(psd.SSh), k = 1:nrow(psd.K)], psd.SSh[s,:Blb] <= b_sk[s, k] <=
            psd.SSh[s,:Bub])
    @constraint(m, [n=1:size(psd.N, 1), k = 1:nrow(psd.K)], pslackm_nk[n, k] >= 0)
    @constraint(m, [n=1:size(psd.N, 1), k = 1:nrow(psd.K)], pslackp_nk[n, k] >= 0)
    @constraint(m, [n=1:size(psd.N, 1), k = 1:nrow(psd.K)], qslackm_nk[n, k] >= 0)
    @constraint(m, [n=1:size(psd.N, 1), k = 1:nrow(psd.K)], qslackp_nk[n, k] >= 0)
    @constraint(m, [l=1:nrow(psd.L), i=1:2, k = 1:nrow(psd.K)], sslack_lik[l, i, k] >= 0)
    @constraint(m, [t=1:nrow(psd.T), i=1:2, k = 1:nrow(psd.K)], sslack_tik[t, i, k] >= 0)

    # Number of contingencies
    krow = nrow(psd.K)

    for k = 1:nrow(psd.K)
        # create generic contingency object
        con = GenericContingency(psd.K[k,:IDout][findall(psd.K[k,:ConType][:] .== :Generator)], 
                                 psd.K[k,:IDout][findall(psd.K[k,:ConType][:] .== :Line)], 
                                 psd.K[k,:IDout][findall(psd.K[k,:ConType][:] .== :Transformer)])

        x0 = get_primal_starting_point(psd, con)

        # contingency variables starting values
        for n=1:nrow(psd.N)
            set_start_value(v_nk[n, k], x0[:v_nk][n])
            set_start_value(theta_nk[n, k], x0[:theta_nk][n])
            set_start_value(pslackm_nk[n, k], x0[:pslackm_nk][n])
            set_start_value(pslackp_nk[n, k], x0[:pslackp_nk][n])
            set_start_value(qslackm_nk[n, k], x0[:qslackm_nk][n])
            set_start_value(qslackp_nk[n, k], x0[:qslackp_nk][n])
        end
        for l=1:nrow(psd.L), i=1:2
            set_start_value(p_lik[l, i, k], x0[:p_lik][l,i])
            set_start_value(q_lik[l, i, k], x0[:q_lik][l,i])
            set_start_value(sslack_lik[l, i, k], x0[:sslack_lik][l,i])
        end
        for t=1:nrow(psd.T), i=1:2
            set_start_value(p_tik[t, i, k], x0[:p_tik][t,i])
            set_start_value(q_tik[t, i, k], x0[:q_tik][t,i])
            set_start_value(sslack_tik[t, i, k], x0[:sslack_tik][t,i])
        end
        for s=1:nrow(psd.SSh)
            set_start_value(b_sk[s, k], x0[:b_sk][s])
        end
        for g=1:nrow(psd.G)
            set_start_value(p_gk[g, k], x0[:p_gk][g])
            set_start_value(q_gk[g, k], x0[:q_gk][g])
        end
        
        # fix angle at reference bus to zero
        JuMP.fix(theta_nk[psd.RefBus, k], 0.0, force=true)

        # add power flow constraints
        addpowerflowcons!(m, v_nk[:,k], theta_nk[:,k], p_lik[:,:,k], q_lik[:,:,k], p_tik[:,:,k], 
                        q_tik[:,:,k], b_sk[:,k],
                        p_gk[:,k], q_gk[:,k], pslackm_nk[:,k], pslackp_nk[:,k], qslackm_nk[:,k],
                        qslackp_nk[:,k], sslack_lik[:,:,k], sslack_tik[:,:,k], psd, con)

        # ramp rate constraints
        Gonline = if length(con.generators_out)>0 setdiff(1:nrow(psd.G), con.generators_out)
                else 1:nrow(psd.G)
                end
                
        @constraint(m, [g in Gonline], psd.G[g,:Plb] <= p_gk[g, k] <= psd.G[g,:Pub])
        @constraint(m, [g in Gonline], psd.G[g,:Qlb] <= q_gk[g, k] <= psd.G[g,:Qub])
        @constraint(m, [g in Gonline], p_gk[g, k] - p_g[g] <=
                                    psd.G[g, :Pub] * psd.G[g, :RampRate] * minutes_since_base)
        @constraint(m, [g in Gonline], p_g[g] - p_gk[g, k] <=
                                    psd.G[g, :Pub] * psd.G[g, :RampRate] * minutes_since_base)

        # enforce out of service generators
        if length(con.generators_out) > 0
            for congo in con.generators_out
                JuMP.fix.(p_gk[findall(psd.G[!, :Generator] .== congo), k], 0.0, force=true)
                JuMP.fix.(q_gk[findall(psd.G[!, :Generator] .== congo), k], 0.0, force=true)
            end
        end

        # contingency penalty
        if !use_huber_like_penalty
            contingency_penalty[k] = JuMP.GenericQuadExpr(JuMP.AffExpr(0))
            for n = 1:nrow(psd.N)
                add_to_expression!(contingency_penalty[k],
                                psd.a[:P], pslackm_nk[n, k], pslackm_nk[n, k])
                add_to_expression!(contingency_penalty[k], psd.b[:P], pslackm_nk[n, k])
                add_to_expression!(contingency_penalty[k],
                                psd.a[:P], pslackp_nk[n, k], pslackp_nk[n, k])
                add_to_expression!(contingency_penalty[k], psd.b[:P], pslackp_nk[n, k])
                add_to_expression!(contingency_penalty[k],
                                psd.a[:Q], qslackm_nk[n, k], qslackm_nk[n, k])
                add_to_expression!(contingency_penalty[k], psd.b[:Q], qslackm_nk[n, k])
                add_to_expression!(contingency_penalty[k],
                                psd.a[:Q], qslackp_nk[n, k], qslackp_nk[n, k])
                add_to_expression!(contingency_penalty[k], psd.b[:Q], qslackp_nk[n, k])
            end
            for l = 1:nrow(psd.L), i=1:2
                add_to_expression!(contingency_penalty[k],
                                psd.a[:S], sslack_lik[l,i,k], sslack_lik[l,i,k])
                add_to_expression!(contingency_penalty[k], psd.b[:S], sslack_lik[l,i,k])
            end
            for t = 1:nrow(psd.T), i=1:2
                add_to_expression!(contingency_penalty[k],
                                psd.a[:S], sslack_tik[t,i,k], sslack_tik[t,i,k])
                add_to_expression!(contingency_penalty[k], psd.b[:S], sslack_tik[t,i,k])
            end
        else
            # collect penalty terms
            p_bal_penalty = @NLexpression(m, sum(hP(pslackm_nk[n,k]) for n=1:nrow(psd.N)) +
                                            sum(hP(pslackp_nk[n,k]) for n=1:nrow(psd.N)))
            q_bal_penalty = @NLexpression(m, sum(hQ(qslackm_nk[n,k]) for n=1:nrow(psd.N)) +
                                            sum(hQ(qslackp_nk[n,k]) for n=1:nrow(psd.N)))
            lin_overload_penalty = @NLexpression(m, sum(hS(sslack_lik[l,i,k]) for l=1:nrow(psd.L), i=1:2))
            trf_overload_penalty = @NLexpression(m, sum(hS(sslack_tik[t,i,k]) for t=1:nrow(psd.T), i=1:2))
            push!(contingency_penalty, @NLexpression(m, p_bal_penalty + q_bal_penalty +
                                                lin_overload_penalty + trf_overload_penalty))
        end

        # quadratic relaxation term
        if quadratic_relaxation_k < Inf
            @assert length(aux_slack_gk) == nrow(psd.G)
            quadratic_relaxation_term[k] = JuMP.GenericQuadExpr(JuMP.AffExpr(0))
            for g = 1:length(aux_slack_gk)
                add_to_expression!(quadratic_relaxation_term[k],
                                quadratic_relaxation_k,
                                aux_slack_gk[g], aux_slack_gk[g])
            end
        else
            push!(quadratic_relaxation_term, 0.0)
        end
    end
    
    # declare objective
    if !use_huber_like_penalty
        @objective(m, Min, production_cost + psd.delta * basecase_penalty +
                            (1-psd.delta)*(1/krow) * ( sum( cp for cp in contingency_penalty) + 
                            sum( qrt for qrt in quadratic_relaxation_term)))
    else
        @NLobjective(m, Min, production_cost + psd.delta * basecase_penalty +
                              (1-psd.delta)*(1/krow) * ( sum( cp for cp in contingency_penalty) + 
                              sum( qrt for qrt in quadratic_relaxation_term)))
    end

    # attempt to solve SCACOPF
    JuMP.optimize!(m)
    if JuMP.primal_status(m) != MOI.FEASIBLE_POINT &&
        JuMP.primal_status(m) != MOI.NEARLY_FEASIBLE_POINT && 
        JuMP.termination_status(m) != MOI.NEARLY_FEASIBLE_POINT
         error("solver failed to find a feasible solution.")
    end

    # objective breakdown
    base_cost = JuMP.value(production_cost) +
                psd.delta*JuMP.value(basecase_penalty)
    recourse_cost = JuMP.objective_value(m) - base_cost


    println("Production Cost: " *string(JuMP.value.(production_cost))*"")
    println("Basecase Penalty: " *string(JuMP.value.(basecase_penalty))*"")
    println("Contingency Penalty: " *string(JuMP.value.(contingency_penalty))*"")

    # Intial construction of the SCACOPF solution
    solution = SCACOPFsolution(psd, BasecaseSolution(psd, JuMP.value.(v_n), JuMP.value.(theta_n),
                                                    convert(Vector{Float64}, JuMP.value.(b_s)),
                                                    JuMP.value.(p_g), JuMP.value.(q_g),
                                                    base_cost, recourse_cost))

    # Add contingency solutions                                                        
    for k = 1:nrow(psd.K)
        con = GenericContingency(psd.K[k,:IDout][findall(psd.K[k,:ConType][:] .== :Generator)], 
                                 psd.K[k,:IDout][findall(psd.K[k,:ConType][:] .== :Line)], 
                                 psd.K[k,:IDout][findall(psd.K[k,:ConType][:] .== :Transformer)])
        contin_penalty = JuMP.value(contingency_penalty[k])
        contingency = ContingencySolution(psd, con,
                                        JuMP.value.(v_nk[:,k]), JuMP.value.(theta_nk[:,k]),
                                        convert(Vector{Float64}, JuMP.value.(b_sk[:,k])),
                                        JuMP.value.(p_gk[:,k]), JuMP.value.(q_gk[:,k]),
                                        0.0, contin_penalty)
        contingency.cont_id = k
        add_contingency_solution!(solution, contingency)
        println("Contingency "*string(k)*" Penalty: " *string(contin_penalty)*"")   
    end

    # return solution
    return solution

end


# function to solve contingency
    
function solve_contingency(psd::SCACOPFdata, k::Int,
                          basecase_solution::BasecaseSolution, NLSolver;
                          previous_solution::Union{Nothing, T} 
                                             where {T<:SubproblemSolution} =
                                             nothing,
                          quadratic_relaxation_k::Float64=Inf,
                          minutes_since_base::Float64=1.0,
                          use_huber_like_penalty::Bool=true)::ContingencySolution
    
    # check we are given a valid contingency index
    if k <= 0 || k > nrow(psd.K)
        error("instance does not have contingency No ", k, " (N cont: ",
              nrow(psd.K), ").")
    end
    
    # create generic contingency object
    # if length(psd.K[k, :ConType]) == 1
    #     con = GenericContingency(psd, k)
    # else
    con = GenericContingency(psd.K[k,:IDout][findall(psd.K[k,:ConType][:] .== :Generator)], 
                                psd.K[k,:IDout][findall(psd.K[k,:ConType][:] .== :Line)], 
                                psd.K[k,:IDout][findall(psd.K[k,:ConType][:] .== :Transformer)])
    # end
    
    # call function for generic contingencies
    sol = solve_contingency(psd, con, basecase_solution, NLSolver,
                            previous_solution=previous_solution,
                            quadratic_relaxation_k=quadratic_relaxation_k,
                            minutes_since_base=minutes_since_base,
                            use_huber_like_penalty=use_huber_like_penalty)
    
    # override fields for contingencies and return (this is TEMPORARY... only for backwards
    # compatibility)
    sol.cont_Kidx = k
    sol.cont_id = psd.K[k,:Contingency]
    sol.cont_idout = psd.K[k,:IDout]
    sol.cont_type = psd.K[k,:ConType]
    sol.cont_alt = nothing
    return sol
    
end

function solve_contingency(psd::SCACOPFdata, con::GenericContingency,
                          basecase_solution::BasecaseSolution, NLSolver;
                          previous_solution::Union{Nothing, T} 
                                             where {T<:SubproblemSolution} =
                                             nothing,
                          quadratic_relaxation_k::Float64=Inf,
                          minutes_since_base::Float64=1.0,
                          use_huber_like_penalty::Bool=true)::ContingencySolution
    

    quadratic_relaxation_k > 0 || error("quadratic relaxation penalty should be positive")
    
    # get primal starting point
    x0 = get_primal_starting_point(psd, con, basecase_solution, previous_solution)
    
    # create model
    m = Model(NLSolver)

    # contingency variables
    @variable(m, psd.N[n,:EVlb] <= v_nk[n=1:nrow(psd.N)] <= psd.N[n,:EVub],
                 start=x0[:v_nk][n])
    @variable(m, theta_nk[n=1:nrow(psd.N)], start=x0[:theta_nk][n])
    @variable(m, p_lik[l=1:nrow(psd.L), i=1:2], start=x0[:p_lik][l,i])
    @variable(m, q_lik[l=1:nrow(psd.L), i=1:2], start=x0[:q_lik][l,i])
    @variable(m, p_tik[t=1:nrow(psd.T), i=1:2], start=x0[:p_tik][t,i])
    @variable(m, q_tik[t=1:nrow(psd.T), i=1:2], start=x0[:q_tik][t,i])
    @variable(m, psd.SSh[s,:Blb] <= b_sk[s=1:nrow(psd.SSh)] <=
              psd.SSh[s,:Bub], start=x0[:b_sk][s])
    @variable(m, p_gk[g=1:nrow(psd.G)], start = x0[:p_gk][g])
    @variable(m, q_gk[g=1:nrow(psd.G)], start = x0[:q_gk][g])
    @variable(m, psd.G[g,:Plb] <= p_g0[g=1:nrow(psd.G)] <= psd.G[g,:Pub],
              start=x0[:p_gk][g])       # clone of first stage variable
    @variable(m, pslackm_nk[n=1:size(psd.N, 1)] >= 0, start = x0[:pslackm_nk][n])
    @variable(m, pslackp_nk[n=1:size(psd.N, 1)] >= 0, start = x0[:pslackp_nk][n])
    @variable(m, qslackm_nk[n=1:size(psd.N, 1)] >= 0, start = x0[:qslackm_nk][n])
    @variable(m, qslackp_nk[n=1:size(psd.N, 1)] >= 0, start = x0[:qslackp_nk][n])
    @variable(m, sslack_lik[l=1:nrow(psd.L), i=1:2] >= 0,
              start=x0[:sslack_lik][l,i])
    @variable(m, sslack_tik[t=1:nrow(psd.T), i=1:2] >= 0,
              start=x0[:sslack_tik][t,i])
    
    # fix angle at reference bus to zero
    JuMP.fix(theta_nk[psd.RefBus], 0.0, force=true)
    
    # add power flow constraints
    addpowerflowcons!(m, v_nk, theta_nk, p_lik, q_lik, p_tik, q_tik, b_sk,
                      p_gk, q_gk, pslackm_nk, pslackp_nk, qslackm_nk,
                      qslackp_nk, sslack_lik, sslack_tik, psd, con)
    
    # ramp rate constraints
    Gonline = if length(con.generators_out)>0 setdiff(1:nrow(psd.G), con.generators_out)
              else 1:nrow(psd.G)
              end

    @constraint(m, [g in Gonline], psd.G[g,:Plb] <= p_gk[g] <= psd.G[g,:Pub])
    @constraint(m, [g in Gonline], psd.G[g,:Qlb] <= q_gk[g] <= psd.G[g,:Qub])
    @constraint(m, [g in Gonline], p_gk[g] - p_g0[g] <=
                                   psd.G[g, :Pub] * psd.G[g, :RampRate] * minutes_since_base)
    @constraint(m, [g in Gonline], p_g0[g] - p_gk[g] <=
                                   psd.G[g, :Pub] * psd.G[g, :RampRate] * minutes_since_base)
    
    # enforce out of service generators
    if length(con.generators_out) > 0
        for congo in con.generators_out
            JuMP.fix.(p_gk[findall(psd.G[!, :Generator] .== congo)], 0.0, force=true)
            JuMP.fix.(q_gk[findall(psd.G[!, :Generator] .== congo)], 0.0, force=true)
        end
    end
    
    # coupling constraints (relaxed non-anticipativity)
    if quadratic_relaxation_k < Inf
        @variable(m, aux_slack_gk[g=1:nrow(psd.G)])
        @constraint(m, non_anticipativity_con[g=1:nrow(psd.G)],
                    p_g0[g] - basecase_solution.p_g[g] == aux_slack_gk[g])
    else
        @constraint(m, non_anticipativity_con[g=1:nrow(psd.G)], p_g0[g] == basecase_solution.p_g[g])
    end
    
    # contingency penalty
    contingency_penalty = nothing
    if !use_huber_like_penalty
        contingency_penalty = JuMP.GenericQuadExpr(JuMP.AffExpr(0))
        for n = 1:nrow(psd.N)
            add_to_expression!(contingency_penalty,
                               psd.a[:P], pslackm_nk[n], pslackm_nk[n])
            add_to_expression!(contingency_penalty, psd.b[:P], pslackm_nk[n])
            add_to_expression!(contingency_penalty,
                               psd.a[:P], pslackp_nk[n], pslackp_nk[n])
            add_to_expression!(contingency_penalty, psd.b[:P], pslackp_nk[n])
            add_to_expression!(contingency_penalty,
                               psd.a[:Q], qslackm_nk[n], qslackm_nk[n])
            add_to_expression!(contingency_penalty, psd.b[:Q], qslackm_nk[n])
            add_to_expression!(contingency_penalty,
                               psd.a[:Q], qslackp_nk[n], qslackp_nk[n])
            add_to_expression!(contingency_penalty, psd.b[:Q], qslackp_nk[n])
        end
        for l = 1:nrow(psd.L), i=1:2
            add_to_expression!(contingency_penalty,
                               psd.a[:S], sslack_lik[l,i], sslack_lik[l,i])
            add_to_expression!(contingency_penalty, psd.b[:S], sslack_lik[l,i])
        end
        for t = 1:nrow(psd.T), i=1:2
            add_to_expression!(contingency_penalty,
                               psd.a[:S], sslack_tik[t,i], sslack_tik[t,i])
            add_to_expression!(contingency_penalty, psd.b[:S], sslack_tik[t,i])
        end
    else
        # register Huber-like penalty functions
        hP = HuberLikePenalty(psd.a[:P], psd.b[:P],
                              2 * psd.a[:P] * (10.0/psd.MVAbase) + psd.b[:P], # slope at 10MW is max
                              5.0/psd.MVAbase)                                # 5MW to change curvature
        hP_prime = HuberLikePenaltyPrime(hP)
        hP_prime_prime = HuberLikePenaltyPrimePrime(hP)
        register(m, :hP, 1, x -> hP(x), x -> hP_prime(x), x -> hP_prime_prime(x))
        hQ = HuberLikePenalty(psd.a[:Q], psd.b[:Q],
                              2 * psd.a[:Q] * (10.0/psd.MVAbase) + psd.b[:Q], # slope at 10MVAr is max
                              5.0/psd.MVAbase)                                # 5MVAr to change curvature
        hQ_prime = HuberLikePenaltyPrime(hQ)
        hQ_prime_prime = HuberLikePenaltyPrimePrime(hQ)
        register(m, :hQ, 1, x -> hQ(x), x -> hQ_prime(x), x -> hQ_prime_prime(x))
        hS = HuberLikePenalty(psd.a[:S], psd.b[:S],
                              2 * psd.a[:S] * (10.0/psd.MVAbase) + psd.b[:S], # slope at 10MVA is max
                              5.0/psd.MVAbase)                                # 5MVA to change curvature
        hS_prime = HuberLikePenaltyPrime(hS)
        hS_prime_prime = HuberLikePenaltyPrimePrime(hS)
        register(m, :hS, 1, x -> hS(x), x -> hS_prime(x), x -> hS_prime_prime(x))
        # collect penalty terms
        p_bal_penalty = @NLexpression(m, sum(hP(pslackm_nk[n]) for n=1:nrow(psd.N)) +
                                         sum(hP(pslackp_nk[n]) for n=1:nrow(psd.N)))
        q_bal_penalty = @NLexpression(m, sum(hQ(qslackm_nk[n]) for n=1:nrow(psd.N)) +
                                         sum(hQ(qslackp_nk[n]) for n=1:nrow(psd.N)))
        lin_overload_penalty = @NLexpression(m, sum(hS(sslack_lik[l,i]) for l=1:nrow(psd.L), i=1:2))
        trf_overload_penalty = @NLexpression(m, sum(hS(sslack_tik[t,i]) for t=1:nrow(psd.T), i=1:2))
        contingency_penalty = @NLexpression(m, p_bal_penalty + q_bal_penalty +
                                               lin_overload_penalty + trf_overload_penalty)
    end
    
    # quadratic relaxation term
    if quadratic_relaxation_k < Inf
        @assert length(aux_slack_gk) == nrow(psd.G)
        quadratic_relaxation_term = JuMP.GenericQuadExpr(JuMP.AffExpr(0))
        for g = 1:length(aux_slack_gk)
            add_to_expression!(quadratic_relaxation_term,
                               quadratic_relaxation_k,
                               aux_slack_gk[g], aux_slack_gk[g])
        end
    else
        quadratic_relaxation_term = 0.0
    end
    
    # declare objective
    if !use_huber_like_penalty
        @objective(m, Min, contingency_penalty + quadratic_relaxation_term)
    else
        @NLobjective(m, Min, contingency_penalty + quadratic_relaxation_term)
    end
    
    # attempt to solve contingency subproblem
    JuMP.optimize!(m)
    if JuMP.primal_status(m) != MOI.FEASIBLE_POINT &&
       JuMP.primal_status(m) != MOI.NEARLY_FEASIBLE_POINT && 
       JuMP.termination_status(m) != MOI.NEARLY_FEASIBLE_POINT
        error("solver failed to find a feasible solution.")
    end
    
    # compute gradient/subgradient
    if quadratic_relaxation_k < Inf
        obj_grad = 2 * quadratic_relaxation_k * (basecase_solution.p_g[g] - JuMP.value.(p_g0))
    else
        if has_duals(m)
            obj_grad = JuMP.dual.(non_anticipativity_con)
                       # check sign consistency before using these duals; JuMP does not use the same
                       # duality conventions as most OR modeling software
        else
            obj_grad = nothing
        end
    end
    
    # return solution
    return ContingencySolution(psd, con,
                            JuMP.value.(v_nk), JuMP.value.(theta_nk),
                            convert(Vector{Float64}, JuMP.value.(b_sk)),
                            JuMP.value.(p_gk), JuMP.value.(q_gk),
                            0.0, objective_value(m), obj_grad)
    
end

# function to solve a random contingency of up to n_failures elements of the grid

function solve_random_contingency(psd::SCACOPFdata, n_failures::Int, 
                                  basecase_solution::BasecaseSolution, NLSolver;
                                  failure_probability_single_element::Float64=0.01,
                                  rng=Random.GLOBAL_RNG,
                                  previous_solution::Union{Nothing, T} 
                                                     where {T<:SubproblemSolution} =
                                                     nothing,
                                  quadratic_relaxation_k::Float64=Inf,
                                  minutes_since_base::Float64=1.0,
                                  use_huber_like_penalty::Bool=true)::ContingencySolution
    
    # generate random contingency
    con = random_up_to_n_elem_contingency(psd, n_failures,
                                          failure_probability_single_element, rng)
    # NOTE: this function draws uniformly from all failures with less than n_failures,
    # assuming that each element has independent failure bernoulli distribution; it does this by
    # sampling from a conditional binomial distribution.
    
    # solve and return solution (solution struct contains information about the generated solution
    return solve_contingency(psd, con, basecase_solution, NLSolver,
                             previous_solution=previous_solution,
                             quadratic_relaxation_k=quadratic_relaxation_k,
                             minutes_since_base=minutes_since_base,
                             use_huber_like_penalty=use_huber_like_penalty)
    
end

end # this is the end of the module
