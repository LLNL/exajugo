# function (method) to compute full base case solution

function get_full_solution(psd::SCACOPFdata, sol::BasecaseSolution)
    if hash(psd) != sol.psd_hash
        error("base case solution does not correspond to power system data.")
    end
    flows_and_slacks = get_power_flows_and_slacks(psd, sol)
    c_g = get_production_cost(psd, sol)
    return flows_and_slacks..., c_g
end

# function (method) to compute full contingency solution

function get_full_solution(psd::SCACOPFdata, sol::ContingencySolution)
    if hash(psd) != sol.psd_hash
        error("contingency solution does not correspond to power system data.")
    end
    return get_power_flows_and_slacks(psd, sol)
end

# function to compute full initial solution for base case

function get_full_initial_solution(psd::SCACOPFdata)
    sol = BasecaseSolution(psd, psd.N[!,:v0], psd.N[!,:theta0],
                           psd.SSh[!,:b0], psd.G[!,:p0], psd.G[!,:q0],
                           0.0, 0.0)
    return sol.v_n, sol.theta_n, sol.b_s, sol.p_g, sol.q_g,
            get_full_solution(psd, sol)...
end

# function to compute full initial solution for contingency

function get_full_initial_solution(psd::SCACOPFdata, con::GenericContingency)
    p0_g = copy(psd.G[!,:p0])
    q0_g = copy(psd.G[!,:q0])
    if length(con.generators_out) > 0
        for congo in con.generators_out
            p0_g[findall(psd.G[!, :Generator] .== congo)] .= 0.0
            q0_g[findall(psd.G[!, :Generator] .== congo)] .= 0.0
        end
    end
    sol = ContingencySolution(psd, con, psd.N[!,:v0], psd.N[!,:theta0],
                              psd.SSh[!,:b0], p0_g, q0_g, 0.0, 0.0)
    return sol.v_n, sol.theta_n, sol.b_s, sol.p_g, sol.q_g,
           get_full_solution(psd, sol)...
end

# function to compute full initial solution for contingency using a base case solution

function get_full_initial_solution(psd::SCACOPFdata, con::GenericContingency,
                                   sol::BasecaseSolution)
    if hash(psd) != sol.psd_hash
        error("base case solution does not correspond to power system data.")
    end
    p_g = copy(sol.p_g)
    q_g = copy(sol.q_g)
    if length(con.generators_out) > 0
        for congo in con.generators_out
            p_g[findall(psd.G[!, :Generator] .== congo)] .= 0.0
            q_g[findall(psd.G[!, :Generator] .== congo)] .= 0.0
        end
    end
    sol = ContingencySolution(psd, con, sol.v_n, sol.theta_n, sol.b_s,
                              p_g, q_g, 0.0, 0.0)
    return sol.v_n, sol.theta_n, sol.b_s, p_g, q_g,
           get_full_solution(psd, sol)...
end

## auxiliary functions

# power flow on branch

function powerflow(vi::Real, vj::Real, deltaij::Real,
                   A::Real, B::Real, C::Real, Theta::Real=0)::Real
    return A*vi^2 + B*vi*vj*cos(deltaij + Theta) + C*vi*vj*sin(deltaij + Theta)
end

# function to assign flows to short circuit branches

function assign_short_circuit_flows(psd::SCACOPFdata,
                                    v_n::Vector{Float64}, pslack_n::Vector{Float64}, qslack_n::Vector{Float64},
                                    L_short_circuit::Vector{Int}, T_short_circuit::Vector{Int}, RateSymb::Symbol)
    
    # build short circuit (sc) network indices
    Nidx_sc = sort(unique([vec(psd.L_Nidx[L_short_circuit, :]);
                           vec(psd.T_Nidx[T_short_circuit, :])]))
    L_Nidx_sc = Array{Int, 2}(undef, length(L_short_circuit), 2)
    L_Nidx_sc[:, 1] .= indexin(psd.L_Nidx[L_short_circuit, 1], Nidx_sc)
    L_Nidx_sc[:, 2] .= indexin(psd.L_Nidx[L_short_circuit, 2], Nidx_sc)
    T_Nidx_sc = Array{Int, 2}(undef, length(T_short_circuit), 2)
    T_Nidx_sc[:, 1] .= indexin(psd.T_Nidx[T_short_circuit, 1], Nidx_sc)
    T_Nidx_sc[:, 2] .= indexin(psd.T_Nidx[T_short_circuit, 2], Nidx_sc)
    
    # backward mappings
    Lidxn_sc = Vector{Vector{Int}}(undef, length(Nidx_sc))
    Lsgnn_sc = Vector{Vector{Int}}(undef, length(Nidx_sc))
    Tidxn_sc = Vector{Vector{Int}}(undef, length(Nidx_sc))
    Tsgnn_sc = Vector{Vector{Int}}(undef, length(Nidx_sc))
    for i = 1:length(Nidx_sc)
        Lidxn_sc[i] = Int[]
        Lsgnn_sc[i] = Int[]
        Tidxn_sc[i] = Int[]
        Tsgnn_sc[i] = Int[]
    end
    for i = 1:length(L_short_circuit)
        push!(Lidxn_sc[L_Nidx_sc[i,1]], i)
        push!(Lsgnn_sc[L_Nidx_sc[i,1]], 1)
        push!(Lidxn_sc[L_Nidx_sc[i,2]], i)
        push!(Lsgnn_sc[L_Nidx_sc[i,2]], -1)
    end
    for i = 1:length(T_short_circuit)
        push!(Tidxn_sc[T_Nidx_sc[i,1]], i)
        push!(Tsgnn_sc[T_Nidx_sc[i,1]], 1)
        push!(Tidxn_sc[T_Nidx_sc[i,2]], i)
        push!(Tsgnn_sc[T_Nidx_sc[i,2]], -1)
    end
    
    # create mathematical programming model to minimize total slacks
	opt = optimizer_with_attributes(Ipopt.Optimizer, "sb" => "yes")
    m = Model(opt)
    
    # declare variables: slacks and short circuit flows
    @variable(m, pslackm_n_sc[1:length(Nidx_sc)] >= 0)
    @variable(m, pslackp_n_sc[1:length(Nidx_sc)] >= 0)
    @variable(m, qslackm_n_sc[1:length(Nidx_sc)] >= 0)
    @variable(m, qslackp_n_sc[1:length(Nidx_sc)] >= 0)
    @variable(m, p_l_sc[1:length(L_short_circuit)])
    @variable(m, q_l_sc[1:length(L_short_circuit)])
    @variable(m, sslack_l_sc[1:length(L_short_circuit)] >= 0)
    @variable(m, p_t_sc[1:length(T_short_circuit)])
    @variable(m, q_t_sc[1:length(T_short_circuit)])
    @variable(m, sslack_t_sc[1:length(T_short_circuit)] >= 0)
    
    # short circuit line constraints
    vmin_l_sc = [minimum(v_n[vec(psd.L_Nidx[i,:])]) for i in L_short_circuit]
    for i=1:length(L_short_circuit)
        if isinf(psd.L[L_short_circuit[i],RateSymb])
            continue
        end
        @constraint(m, p_l_sc[i]^2 + q_l_sc[i]^2 <= 
                    (psd.L[L_short_circuit[i],RateSymb]*vmin_l_sc[i] + sslack_l_sc[i])^2)
    end
    
    # short circuit transformer constraints
    for i=1:length(T_short_circuit)
        if isinf(psd.T[T_short_circuit[i],RateSymb])
            continue
        end
        @constraint(m, p_t_sc[i]^2 + q_t_sc[i]^2 <= 
                       (psd.T[T_short_circuit[i],RateSymb] + sslack_t_sc[i])^2)
    end
    
    # slack balance constraints
    for i = 1:length(Nidx_sc)
        @constraint(m, pslack_n[Nidx_sc[i]] -
                       sum(Lsgnn_sc[i][j] * p_l_sc[Lidxn_sc[i][j]] for j=1:length(Lidxn_sc[i])) -
                       sum(Tsgnn_sc[i][j] * p_t_sc[Tidxn_sc[i][j]] for j=1:length(Tidxn_sc[i])) ==
                       pslackp_n_sc[i] - pslackm_n_sc[i])
        @constraint(m, qslack_n[Nidx_sc[i]] -
                       sum(Lsgnn_sc[i][j] * q_l_sc[Lidxn_sc[i][j]] for j=1:length(Lidxn_sc[i])) -
                       sum(Tsgnn_sc[i][j] * q_t_sc[Tidxn_sc[i][j]] for j=1:length(Tidxn_sc[i])) ==
                       qslackp_n_sc[i] - qslackm_n_sc[i])
    end
    
    # objective: minimize slacks
    @objective(m, Min, sum(pslackm_n_sc[i] for i=1:length(Nidx_sc)) + 
                       sum(pslackp_n_sc[i] for i=1:length(Nidx_sc)) +
                       sum(qslackm_n_sc[i] for i=1:length(Nidx_sc)) +
                       sum(qslackp_n_sc[i] for i=1:length(Nidx_sc)) +
                       2.0 * sum(sslack_l_sc[i] for i=1:length(L_short_circuit)) +
                       2.0 * sum(sslack_t_sc[i] for i=1:length(T_short_circuit)))
    
    # solve model
    JuMP.optimize!(m)
    if JuMP.primal_status(m) != MOI.FEASIBLE_POINT && 
       JuMP.primal_status(m) != MOI.NEARLY_FEASIBLE_POINT
        error("solver failed to find a feasible solution.")
    end
    
    # return solution: from bus flow at each short circuit branch
    return JuMP.value.(p_l_sc), JuMP.value.(q_l_sc),
           JuMP.value.(p_t_sc), JuMP.value.(q_t_sc)
    
end

# computing power flows and slacks

function get_power_flows_and_slacks(psd::SCACOPFdata, sol::SubproblemSolution)
    
    # determine thermal rate type
    if typeof(sol) <: BasecaseSolution
        RateSymb = :RateBase
    elseif typeof(sol) <: ContingencySolution
        RateSymb = :RateEmer
    else
        @warn "undetermined thermal rate type for solution type " *
              string(typeof(sol)) * ". Will assume RateBase."
        RateSymb = :RateBase
    end
       
    # compute line flows and slacks
    p_li = Array{Float64}(undef, nrow(psd.L), 2)
    q_li = Array{Float64}(undef, nrow(psd.L), 2)
    sslack_li = Array{Float64}(undef, nrow(psd.L), 2)
    L_short_circuit = Int[]
    for l=1:nrow(psd.L), i=1:2
        y_l = psd.L[l,:G] + im * psd.L[l,:B]
        if isinf(abs(y_l))
            push!(L_short_circuit, l)
            p_li[l, i] = 0.0
            q_li[l, i] = 0.0
            sslack_li[l,i] = 0.0
            continue
        end
        vi = sol.v_n[psd.L_Nidx[l,i]]
        vj = sol.v_n[psd.L_Nidx[l,3-i]]
        deltaij = sol.theta_n[psd.L_Nidx[l,i]] - sol.theta_n[psd.L_Nidx[l,3-i]]
        p_li[l,i] = powerflow(vi, vj, deltaij,
                              psd.L[l,:G], -psd.L[l,:G], -psd.L[l,:B])
        q_li[l,i] = powerflow(vi, vj, deltaij,
                              -psd.L[l,:B] - psd.L[l,:Bch]/2,
                              psd.L[l,:B], -psd.L[l,:G])
        sslack_li[l,i] = max(0.0, sqrt(p_li[l,i]^2 + q_li[l,i]^2) -
                             psd.L[l,RateSymb]*vi)
    end
    unique!(L_short_circuit)
    
    # compute transformer flows and slacks
    p_ti = Array{Float64}(undef, nrow(psd.T), 2)
    q_ti = Array{Float64}(undef, nrow(psd.T), 2)
    sslack_ti = Array{Float64}(undef, nrow(psd.T), 2)
    T_short_circuit = Int[]
    for t=1:nrow(psd.T)
        y_t = psd.T[t, :G] + im * psd.T[t, :B]
        if isinf(abs(y_t))
            push!(T_short_circuit, t)
            p_ti[t,:] .= 0.0
            q_ti[t,:] .= 0.0
            sslack_ti[t,:] .= 0.0
            continue
        end
        v1 = sol.v_n[psd.T_Nidx[t,1]]
        v2 = sol.v_n[psd.T_Nidx[t,2]]
        delta12 = sol.theta_n[psd.T_Nidx[t,1]] - sol.theta_n[psd.T_Nidx[t,2]]
        # i == 1
        p_ti[t,1] = powerflow(v1, v2, delta12,
                              psd.T[t,:G]/psd.T[t,:Tau]^2 + psd.T[t,:Gm],
                              -psd.T[t,:G]/psd.T[t,:Tau],
                              -psd.T[t,:B]/psd.T[t,:Tau], -psd.T[t,:Theta])
        q_ti[t,1] = powerflow(v1, v2, delta12,
                              -psd.T[t,:B]/psd.T[t,:Tau]^2 - psd.T[t,:Bm],
                              psd.T[t,:B]/psd.T[t,:Tau],
                              -psd.T[t,:G]/psd.T[t,:Tau], -psd.T[t,:Theta])
        # i == 2
        p_ti[t,2] = powerflow(v2, v1, -delta12,
                              psd.T[t,:G], -psd.T[t,:G]/psd.T[t,:Tau],
                              -psd.T[t,:B]/psd.T[t,:Tau], psd.T[t,:Theta])
        q_ti[t,2] = powerflow(v2, v1, -delta12,
                              -psd.T[t,:B], psd.T[t,:B]/psd.T[t,:Tau],
                              -psd.T[t,:G]/psd.T[t,:Tau], psd.T[t,:Theta])
        for i=1:2
            sslack_ti[t,i] = max(0.0, sqrt(p_ti[t,i]^2 + q_ti[t,i]^2) -
                                 psd.T[t,RateSymb])
        end
    end
    
    # set power flows and slacks to zero for failed branches
    if typeof(sol) <: ContingencySolution
        @assert typeof(sol.cont_alt) <: GenericContingency
        if length(sol.cont_alt.lines_out) > 0
            for conalo in sol.cont_alt.lines_out
                p_li[findall(psd.L[!, :Line] .== conalo), :] .= 0.0
                q_li[findall(psd.L[!, :Line] .== conalo), :] .= 0.0
                sslack_li[findall(psd.L[!, :Line] .== conalo), :] .= 0.0
            end
        end
        if length(sol.cont_alt.transformers_out) > 0
            for conato in sol.cont_alt.transformers_out
                p_ti[findall(psd.T[!, :Transformer] .== conato), :] .= 0.0
                q_ti[findall(psd.T[!, :Transformer] .== conato), :] .= 0.0
                sslack_ti[findall(psd.T[!, :Transformer] .== conato), :] .= 0.0
            end
        end
    end
    
    # substract contingency lines and transformers from short circuit records
    if typeof(sol) <: ContingencySolution
        setdiff!(L_short_circuit, sol.cont_alt.lines_out)
        setdiff!(T_short_circuit, sol.cont_alt.transformers_out)
    end
    
    # pre compute power balance slacks
    pslackm_n = Vector{Float64}(undef, nrow(psd.N))
    pslackp_n = Vector{Float64}(undef, nrow(psd.N))
    qslackm_n = Vector{Float64}(undef, nrow(psd.N))
    qslackp_n = Vector{Float64}(undef, nrow(psd.N))
    for n=1:nrow(psd.N)
        pslack = 0.0
        qslack = 0.0
        if length(psd.Gn[n]) > 0
            pslack += sum(sol.p_g[g] for g=psd.Gn[n])
            qslack += sum(sol.q_g[g] for g=psd.Gn[n])
        end
        pslack -= psd.N[n,:Pd] + psd.N[n,:Gsh]*sol.v_n[n]^2
        qslack -= psd.N[n,:Qd] - psd.N[n,:Bsh]*sol.v_n[n]^2
        if length(psd.SShn[n]) > 0
            qslack -= sum(-sol.b_s[s] for s=psd.SShn[n])*sol.v_n[n]^2
        end
        if length(psd.Lidxn[n]) > 0
            pslack -= sum(p_li[psd.Lidxn[n][lix], psd.Lin[n][lix]]
                          for lix=1:length(psd.Lidxn[n]))
            qslack -= sum(q_li[psd.Lidxn[n][lix], psd.Lin[n][lix]]
                          for lix=1:length(psd.Lidxn[n]))
        end
        if length(psd.Tidxn[n]) > 0
            pslack -= sum(p_ti[psd.Tidxn[n][tix], psd.Tin[n][tix]]
                          for tix=1:length(psd.Tidxn[n]))
            qslack -= sum(q_ti[psd.Tidxn[n][tix], psd.Tin[n][tix]]
                          for tix=1:length(psd.Tidxn[n]))
        end
        pslackm_n[n] = max(0.0, -pslack)
        pslackp_n[n] = max(0.0, pslack)
        qslackm_n[n] = max(0.0, -qslack)
        qslackp_n[n] = max(0.0, qslack)
    end
    
    # attemp to reduce slacks using short circuit lines
    if length(L_short_circuit) + length(T_short_circuit) > 0
        
        # compute flows that would minimize slacks
        p_l_sc, q_l_sc, p_t_sc, q_t_sc = 
            assign_short_circuit_flows(psd, sol.v_n, pslackp_n - pslackm_n, qslackp_n - qslackm_n,
                                       L_short_circuit, T_short_circuit, RateSymb)
        
        # set flows for short circuit branches
        for k = 1:length(L_short_circuit)
            l = L_short_circuit[k]
            p_li[l, 1] = p_l_sc[k]
            p_li[l, 2] = -1.0 * p_l_sc[k]
            q_li[l, 1] = q_l_sc[k]
            q_li[l, 2] = -1.0 * q_l_sc[k]
            for i = 1:2
                vi = sol.v_n[psd.L_Nidx[l,i]]
                sslack_li[l,i] = max(0.0, sqrt(p_li[l,i]^2 + q_li[l,i]^2) -
                                     psd.L[l,RateSymb]*vi)
            end
        end
        for k = 1:length(T_short_circuit)
            t = T_short_circuit[k]
            p_ti[t, 1] = p_t_sc[k]
            p_ti[t, 2] = -1.0 * p_t_sc[k]
            q_ti[t, 1] = q_t_sc[k]
            q_ti[t, 2] = -1.0 * q_t_sc[k]
            for i=1:2
                sslack_ti[t,i] = max(0.0, sqrt(p_ti[t,i]^2 + q_ti[t,i]^2) -
                                     psd.T[t,RateSymb])
            end
        end
        
        # update imbalance slacks based on these flows
        for k = 1:length(L_short_circuit), i=1:2
            l = L_short_circuit[k]
            n = psd.L_Nidx[l, i]
            pslack_n = pslackp_n[n] - pslackm_n[n] - p_li[l, i]
            pslackm_n[n] = max(0, -pslack_n) 
            pslackp_n[n] = max(0, pslack_n)
            qslack_n = qslackp_n[n] - qslackm_n[n] - q_li[l, i]
            qslackm_n[n] = max(0, -qslack_n) 
            qslackp_n[n] = max(0, qslack_n)
        end
        for k = 1:length(T_short_circuit), i=1:2
            t = T_short_circuit[k]
            n = psd.T_Nidx[t, i]
            pslack_n = pslackp_n[n] - pslackm_n[n] - p_ti[t, i]
            pslackm_n[n] = max(0, -pslack_n) 
            pslackp_n[n] = max(0, pslack_n)
            qslack_n = qslackp_n[n] - qslackm_n[n] - q_ti[t, i]
            qslackm_n[n] = max(0, -qslack_n) 
            qslackp_n[n] = max(0, qslack_n)
        end
        
    end
    
    # return power flows and slacks
    return p_li, q_li, p_ti, q_ti,
           pslackm_n, pslackp_n, qslackm_n, qslackp_n,
           sslack_li, sslack_ti
    
end

# function to compute production cost

function get_production_cost(psd::SCACOPFdata, sol::BasecaseSolution)
    c_g = Vector{Float64}(undef, nrow(psd.G))
    if psd.G.CTYP[1] == 1
        for g = 1:nrow(psd.G)
            c_g[g]  = psd.G.COSTQUAD[g]    * sol.p_g[g]^2 * psd.MVAbase * psd.MVAbase
                        + psd.G.COSTLIN[g] * sol.p_g[g]   * psd.MVAbase
                        + psd.G.COST[g]
        end
    else
        for g = 1:nrow(psd.G)
            c_g[g] = maximum(sol.p_g[g]*psd.G_epicost_slope[g] +
                            psd.G_epicost_intercept[g])
        end
    end
    return c_g
end